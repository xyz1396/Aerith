#include "lib/peptide.h"
#include "lib/ms2scan.h"
#include "lib/initSIP.h"
#include "lib/PSMpeakAnnotator.h"
#include <Rcpp.h>
using namespace Rcpp;

class ms2scanWithNewScoreFunction : public MS2Scan
{
public:
    double scoreIntensitySimple(const double realIntensity, const double expectedIntensity)
    {
        return (1 - std::erf(std::abs(realIntensity - expectedIntensity) /
                             std::sqrt((realIntensity * realIntensity + expectedIntensity * expectedIntensity) / 2)));
    };

    // get pearson coefficient, A B should be the same size
    double pearson(const vector<double> &A, const vector<double> &B)
    {

        double sumA(0.0), sumB(0.0), aveA(0.0), aveB(0.0);
        size_t Length = A.size();

        // sum
        sumA = std::accumulate(A.begin(), A.end(), 0.0);
        sumB = std::accumulate(B.begin(), B.end(), 0.0);

        // average
        aveA = sumA / double(Length);
        aveB = sumB / double(Length);

        // pearson coefficient
        double R1(0), R2(0), R3(0);
        for (size_t i = 0; i < Length; i++)
        {
            R1 += (A[i] - aveA) * (B[i] - aveB);
            R2 += std::pow((A[i] - aveA), 2);
            R3 += std::pow((B[i] - aveB), 2);
        }
        // avoid divide by 0
        double R4 = R2 * R3;
        // for test
        // Rcout << R1 << "\t" << R1 / std::sqrt(R4) << endl;
        if (R4 == 0)
            return 0;
        else
            return (R1 / std::sqrt(R4));
    }

    double cosine_similarity(const vector<double> &A, const vector<double> &B)
    {
        double dot = 0.0, denom_a = 0.0, denom_b = 0.0;
        size_t Length = A.size();
        for (size_t i = 0; i < Length; ++i)
        {
            dot += A[i] * B[i];
            denom_a += A[i] * A[i];
            denom_b += B[i] * B[i];
        }
        double R = std::sqrt(denom_a) * std::sqrt(denom_b);
        if (R == 0)
            return 0;
        else
            return dot / R;
    }

    double crossEntropy(const vector<double> &A, const vector<double> &B)
    {
        double loss = 0, sumA = 0, sumB = 0;
        sumA = std::accumulate(A.begin(), A.end(), 0.0);
        sumB = std::accumulate(B.begin(), B.end(), 0.0);
        double Arelative = 0, Brelative = 0;
        // put B behind because B has 0
        for (size_t i = 0; i < A.size(); ++i)
        {
            Brelative = B[i] / sumB;
            Arelative = A[i] / sumA;
            loss += -Brelative * std::log(Arelative) - (1.0 - Brelative) * std::log(1.0 - Arelative);
        }
        return loss / (double)A.size();
    }

    // this method is best for score isotopic peaks envelope
    double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity)
    {
        if (observed)
            return 0.5 * (1 - std::erf(std::abs(realIntensity - expectedIntensity) /
                                       std::sqrt((realIntensity * realIntensity + expectedIntensity * expectedIntensity) / 2)));
        else
        {
            // deductionCoefficient is -0.55 when 13Cpct=0 , -0.05 when 13Cpct=0.5
            // Rcout << ProNovoConfig::getDeductionCoefficient() << endl;
            return ProNovoConfig::getDeductionCoefficient() * expectedIntensity;
        }
    };

    bool findProductIonSIP(const vector<double> &vdIonMass,
                           const vector<double> &vdIonProb,
                           const int &iCharge,
                           double &dScoreWeight,
                           double &dAverageMZError,
                           double &dMostAbundantObservedMZ,
                           int &iMostAbundantPeakIndex)
    {
        int iIndex4MaxInt = getMaxValueIndex(vdIonProb);
        double dMaxIntExpectedMZ = (vdIonMass[iIndex4MaxInt] / (double)iCharge) + dProtonMass;
        int iIndex4SelectedFound = 0;
        dScoreWeight = 0;
        dAverageMZError = dMassTolerance;
        // search for the most abundant peak
        if (!searchMZ2D(dMaxIntExpectedMZ, iIndex4SelectedFound))
        {
            return false;
        }
        iMostAbundantPeakIndex = iIndex4SelectedFound;
        dMostAbundantObservedMZ = vdpreprocessedMZ[iIndex4SelectedFound];
        double dMostAbundantObservedIntensity = vdpreprocessedIntensity[iIndex4SelectedFound];
        double dMostAbundantMZError = dMostAbundantObservedMZ - dMaxIntExpectedMZ;
        // compute expected MZ and intensity for this product ion
        vector<bool> vbObserved(vdIonProb.size(), false);
        vector<double> vdObservedMZ(vdIonProb.size(), 0);
        vector<double> vdObservedRelativeInt(vdIonProb.size(), 0);
        vector<double> vdMZError(vdIonProb.size(), dMassTolerance);
        vector<double> vdExpectedMZ(vdIonProb.size(), 0);
        vector<double> vdExpectedRelativeInt(vdIonProb.size(), 0);
        int i;
        for (i = 0; i < (int)vdIonProb.size(); ++i)
        {
            vdExpectedMZ[i] = vdIonMass[i] / (double)iCharge + dProtonMass;
            vdExpectedRelativeInt[i] = vdIonProb[i] / vdIonProb[iIndex4MaxInt];
        }

        // set max int value
        vbObserved[iIndex4MaxInt] = true;
        vdObservedMZ[iIndex4MaxInt] = dMostAbundantObservedMZ;
        vdMZError[iIndex4MaxInt] = dMostAbundantMZError;
        vdObservedRelativeInt[iIndex4MaxInt] = 1.0;

        if (vbObserved.size() == 1)
        {
            // there is only one expected peak in the isotopic distribution
            dScoreWeight = 1.0;
            dAverageMZError = dMostAbundantMZError;
            return true;
        }

        // search for  ions on the right of the most abundant peak
        vector<int> viIndex4Found;
        double dShiftedExpectedMZ;
        unsigned int j;
        double minIntensityError = 0, currentIntensityError = 0;
        double dCurrentRelativeInt = 0;
        // for test
        // cout << "current isotopic peak" << endl;
        for (i = iIndex4MaxInt + 1; i < (int)vdExpectedMZ.size(); ++i)
        {
            // shift expected MZ by the error of the most abundant ion and search for it
            dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
            if (binarySearch(dShiftedExpectedMZ, vdpreprocessedMZ, dMassTolerance / 2, viIndex4Found))
            {
                // relative intensity error should be less than 1.5 fold of vdExpectedRelativeInt
                minIntensityError = 1.5 * vdExpectedRelativeInt[i];
                for (j = 0; j < viIndex4Found.size(); ++j)
                {
                    iIndex4SelectedFound = viIndex4Found[j];
                    dCurrentRelativeInt = vdpreprocessedIntensity[iIndex4SelectedFound] / dMostAbundantObservedIntensity;
                    currentIntensityError = std::abs(dCurrentRelativeInt - vdExpectedRelativeInt[i]);
                    // for test
                    // cout << i << j << " " << iCharge << " " << vdpreprocessedMZ[iIndex4SelectedFound] << " " << vdExpectedMZ[i] << endl;
                    if (currentIntensityError < minIntensityError)
                    {
                        vbObserved[i] = true;
                        vdObservedMZ[i] = vdpreprocessedMZ[iIndex4SelectedFound];
                        vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
                        vdObservedRelativeInt[i] = dCurrentRelativeInt;
                        minIntensityError = currentIntensityError;
                        // for test
                        // Rcout << vdExpectedRelativeInt[i] << " " << dCurrentRelativeInt << " " << minIntensityError << endl;
                    }
                }
            }
            else
                break;
        }
        // search for ions on the left of the most abundant peak
        // identical to the above function except the index
        for (i = iIndex4MaxInt - 1; i >= 0; --i)
        {
            // shift expected MZ by the error of the most abundant ion and search for it
            dShiftedExpectedMZ = vdExpectedMZ[i] + dMostAbundantMZError;
            if (binarySearch(dShiftedExpectedMZ, vdpreprocessedMZ, dMassTolerance / 2, viIndex4Found))
            {
                minIntensityError = 1.5 * vdExpectedRelativeInt[i];
                for (j = 0; j < viIndex4Found.size(); ++j)
                {
                    iIndex4SelectedFound = viIndex4Found[j];
                    dCurrentRelativeInt = vdpreprocessedIntensity[iIndex4SelectedFound] / dMostAbundantObservedIntensity;
                    currentIntensityError = std::abs(dCurrentRelativeInt - vdExpectedRelativeInt[i]);
                    if (currentIntensityError < minIntensityError)
                    {
                        vbObserved[i] = true;
                        vdObservedMZ[i] = vdpreprocessedMZ[iIndex4SelectedFound];
                        vdMZError[i] = vdObservedMZ[i] - vdExpectedMZ[i];
                        vdObservedRelativeInt[i] = dCurrentRelativeInt;
                        minIntensityError = currentIntensityError;
                    }
                }
            }
            else
                break;
        }
        // calculate score weight for this product ion
        dScoreWeight = 1.0;
        // vector<double> vdTempExpectedRelativeInt = vdExpectedRelativeInt;
        // vdTempExpectedRelativeInt[iIndex4MaxInt] = 0;
        // int iIndex4SecondHighestInt = getMaxValueIndex(vdTempExpectedRelativeInt);
        // // if the second highest peak is found
        // // if (vbObserved[iIndex4SecondHighestInt])
        // //     dScoreWeight = 2.0;
        // // else
        // //     dScoreWeight = 1.0;
        // add score of isotopic peak intensity
        for (size_t i = 0; i < vbObserved.size(); i++)
        {
            if (i != (size_t)iIndex4MaxInt)
                dScoreWeight += scoreIntensity(vbObserved[i], vdObservedRelativeInt[i], vdExpectedRelativeInt[i]);
        }
        // avoid negtive
        if (dScoreWeight <= 0)
        {
            dScoreWeight = 0;
            return false;
        }
        // test whether the iCharge is consistant with viZinput
        // if not, lower the dScoreWeight
        if (vipreprocessedCharge[iMostAbundantPeakIndex] != 0)
        {
            if (vipreprocessedCharge[iMostAbundantPeakIndex] != iCharge)
                dScoreWeight = dScoreWeight / 2;
        }
        // calculate average mass error for this product ion
        double dTotalRelativeIntensity = 0;
        double dTotalMZError = 0;
        for (i = 0; i < (int)vdMZError.size(); ++i)
        {
            if (vbObserved[i])
            {
                dTotalMZError += vdMZError[i] * vdExpectedRelativeInt[i];
                dTotalRelativeIntensity += vdExpectedRelativeInt[i];
            }
        }
        dAverageMZError = dTotalMZError / dTotalRelativeIntensity;
        //      dAverageMZError = dMostAbundantMZError;
        return true;
    };

    void scoreWeightSumHighMS2(Peptide *currentPeptide)
    {
        double dScore = 0;
        int iPeptideLength = currentPeptide->getPeptideLength(), i, j, iMostAbundantPeakIndex = 0;
        int n; // Ion number starting from one
        int z; // charge state
        vector<ProductIon> vFoundIons;
        double dScoreWeight = 0, dMZError = 1, dMostAbundantObservedMZ = 0, dAverageMZError = 0,
               dBonus4ComplementaryFragmentObserved = 1.0;
        string sPeptide = currentPeptide->getPeptideSeq();
        //    cout<<currentPeptide->vvdYionMass.size()<<endl;
        for (n = 0; n < (int)currentPeptide->vvdYionMass.size(); ++n)
            for (z = 1; z <= iParentChargeState; ++z)
            {
                ProductIon currentIon;
                currentIon.setProductIon('y', n + 1, z);
                if (ProNovoConfig::getSearchType() == "SIP")
                {
                    if (findProductIonSIP(currentPeptide->vvdYionMass[n], currentPeptide->vvdYionProb[n], z,
                                          dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
                    {
                        currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
                        vFoundIons.push_back(currentIon);
                    }
                }
                else
                {
                    if (findProductIon(currentPeptide->vvdYionMass[n], currentPeptide->vvdYionProb[n], z,
                                       dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
                    {
                        currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
                        vFoundIons.push_back(currentIon);
                    }
                }
            }

        for (n = 0; n < (int)currentPeptide->vvdBionMass.size(); ++n)
            for (z = 1; z <= iParentChargeState; ++z)
            {
                ProductIon currentIon;
                currentIon.setProductIon('b', n + 1, z);
                if (ProNovoConfig::getSearchType() == "SIP")
                {
                    if (findProductIonSIP(currentPeptide->vvdBionMass[n], currentPeptide->vvdBionProb[n], z,
                                          dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
                    {
                        currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
                        vFoundIons.push_back(currentIon);
                    }
                }
                else
                {
                    if (findProductIon(currentPeptide->vvdBionMass[n], currentPeptide->vvdBionProb[n], z,
                                       dScoreWeight, dMZError, dMostAbundantObservedMZ, iMostAbundantPeakIndex))
                    {
                        currentIon.setObservedInfo(dMZError, dScoreWeight, dMostAbundantObservedMZ, iMostAbundantPeakIndex);
                        vFoundIons.push_back(currentIon);
                    }
                }
            }

        for (i = 0; i < (int)vFoundIons.size(); ++i)
            vFoundIons[i].setComplementaryFragmentObserved(false);

        for (i = 0; i < (int)vFoundIons.size(); ++i)
            for (j = i + 1; j < (int)vFoundIons.size(); ++j)
                if (vFoundIons[i].getIonNumber() + vFoundIons[j].getIonNumber() == iPeptideLength)
                    if ((vFoundIons[i].getIonType() == 'y' && vFoundIons[j].getIonType() == 'b') || (vFoundIons[i].getIonType() == 'b' && vFoundIons[j].getIonType() == 'y'))
                    {
                        vFoundIons[i].setComplementaryFragmentObserved(true);
                        vFoundIons[j].setComplementaryFragmentObserved(true);
                    }
        for (i = 0; i < (int)vFoundIons.size(); ++i)
        {
            dAverageMZError += vFoundIons[i].getMZError();
            // cout<<vFoundIons[i].getMZError()<<endl;
        }
        dAverageMZError = dAverageMZError / (double)vFoundIons.size();

        //   cout<<vFoundIons.size()<<endl;

        for (i = 0; i < (int)vFoundIons.size(); ++i)
        {
            if (vFoundIons[i].getComplementaryFragmentObserved())
                dBonus4ComplementaryFragmentObserved = 2.0;
            else
                dBonus4ComplementaryFragmentObserved = 1.0;
            if (ProNovoConfig::getSearchType() == "SIP")
                dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError() -
                                                         dAverageMZError)) *
                          vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;
            else
                // no mass error calibration
                dScore += ProNovoConfig::scoreError(fabs(vFoundIons[i].getMZError())) * vFoundIons[i].getScoreWeight() * dBonus4ComplementaryFragmentObserved;

            // cout<<dScore<<endl;
        }

        saveScore(dScore, currentPeptide, vpWeightSumTopPeptides, vdWeightSumAllScores);
    };
};

//' scoreIntensity
//' @param observed this peak is observed or not
//' @param realIntensity real intensity in MS2 scan
//' @param expectedIntensity expected intensity
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
//' @return a score of this intensity match
//' @examples
//' scoreIntensity(TRUE, 1200.0, 1180.0, "C13", 0.02)
//' @export
// [[Rcpp::export]]
double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity,
                      const String &Atom, double Prob)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    // compute residue mass and prob again
    computeResidueMassIntensityAgain(Atom, Prob);
    ms2scanWithNewScoreFunction myScan;
    return myScan.scoreIntensity(observed, realIntensity, expectedIntensity);
}

//' scoreIntensityByCrossEntropy
//' @param expectedIntensity expected intensityreal
//' @param observedIntensity observed intensity in MS2 scan
//' @return numeric, a score of this intensity match
//' @examples
//' scoreIntensityByCE(c(10.0, 20.0, 30.0), c(9.5, 21.0, 28.0))
//' @export
// [[Rcpp::export]]
double scoreIntensityByCE(const NumericVector &expectedIntensity, const NumericVector &observedIntensity)
{
    ms2scanWithNewScoreFunction myScan;
    return myScan.crossEntropy(as<vector<double>>(expectedIntensity), as<vector<double>>(observedIntensity));
}

//' scorePSM
//' @param realMZ mz vector in MS2 scan
//' @param realIntensity intensity vector in MS2 scan
//' @param realCharge charge vector in MS2 scan
//' @param parentCharge int parent charge of MS2 scan
//' @param pepSeq a string of peptide
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
//' @return a score of this PSM
//' @examples
//' demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
//' scan1 <- readOneScanMS2(ftFile = demo_file, 107728)
//' score <- scorePSM(scan1$peaks$mz,
//'         scan1$peaks$intensity, scan1$peaks$charge, 2,
//'         "[HSQVFSTAEDNQSAVTIHVLQGER]", "C13", 0.0107)
//' @export
// [[Rcpp::export]]
double scorePSM(const NumericVector &realMZ, const NumericVector &realIntensity,
                const NumericVector &realCharge, int parentCharge,
                const String &pepSeq, const String &Atom, double Prob)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    // compute residue mass and prob again
    computeResidueMassIntensityAgain(Atom, Prob);
    ProNovoConfig::setDeductionCoefficient();
    // for test
    // Rcout << ProNovoConfig::getSetMinValue() << "\t" << ProNovoConfig::getSetFold() << endl;
    Peptide myPep;
    string sOriginalPeptide, sProteinName;
    int ibeginPos;
    double dPeptideMass;
    char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
    myPep.setPeptide(pepSeq, sOriginalPeptide, sProteinName,
                     ibeginPos, dPeptideMass, cIdentifyPrefix,
                     cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
    map<char, double> mapResidueMass;
    // mapResidueMass is not used in this function
    myPep.preprocessing(true, mapResidueMass);
    ms2scanWithNewScoreFunction myScan;
    myScan.isMS2HighRes = true;
    myScan.vdIntensity = as<vector<double>>(realIntensity);
    myScan.vdMZ = as<vector<double>>(realMZ);
    myScan.viCharge = as<vector<int>>(realCharge);
    myScan.iParentChargeState = parentCharge;
    myScan.preprocess();
    myScan.scoreWeightSumHighMS2(&myPep);
    return myScan.vpWeightSumTopPeptides[0]->dScore;
}

std::vector<std::string> enumsToStrings(std::vector<PSMpeakAnnotator::ionKind> ionKinds)
{
    std::vector<std::string> ionKindStrs(ionKinds.size());
    for (size_t i = 0; i < ionKinds.size(); i++)
    {
        switch (ionKinds[i])
        {
        case PSMpeakAnnotator::ionKind::B:
            ionKindStrs[i] = "B";
            break;
        case PSMpeakAnnotator::ionKind::Y:
            ionKindStrs[i] = "Y";
            break;
        case PSMpeakAnnotator::ionKind::BisotopicPeak:
            ionKindStrs[i] = "BisotopicPeak";
            break;
        case PSMpeakAnnotator::ionKind::YisotopicPeak:
            ionKindStrs[i] = "YisotopicPeak";
            break;
        default:
            ionKindStrs[i] = "UNKNOWN";
            break;
        }
    }
    return ionKindStrs;
}

//' annotatePSM
//' @param realMZ mz vector in MS2 scan
//' @param realIntensity intensity vector in MS2 scan
//' @param realCharge charge vector in MS2 scan
//' @param pepSeq a string of peptide
//' @param charges charges of product ions in consideration
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
//' @param isoCenter isolation window center, set it 0 as default if not remove peaks in isolation window
//' @param isoWidth isolation window width, set it 0 as default if not remove peaks in isolation window
//' @param calScores, FALSE as default, calculate WDP MVH Xcor scores or not
//' @return a List about matched peaks information of this PSM
//' @examples
//' demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
//' scan1 <- readOneScanMS2(ftFile = demo_file, 107728)
//' anno <- annotatePSM(
//'   scan1$peaks$mz, scan1$peaks$intensity,
//'   scan1$peaks$charge,
//'   "HSQVFSTAEDNQSAVTIHVLQGER", 1:2, "C13",
//'   0.0107, 886.65, 4.0, TRUE
//' )
//' @export
// [[Rcpp::export]]
List annotatePSM(const NumericVector &realMZ, const NumericVector &realIntensity,
                 const NumericVector &realCharge, const String &pepSeq, const NumericVector charges,
                 const String &Atom, double Prob,
                 const double isoCenter = 0, const double isoWidth = 0, const bool calScores = false)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    // compute residue mass and prob again
    computeResidueMassIntensityAgain(Atom, Prob);
    Scan mScan;
    mScan.mz = as<vector<double>>(realMZ);
    mScan.intensity = as<vector<double>>(realIntensity);
    mScan.charge = as<vector<int>>(realCharge);
    // set tolerance in ppm
    PSMpeakAnnotator mAnnotator(10);
    if (mScan.mz.size() > 10 && mScan.intensity.size() > 10)
        mAnnotator.analyzePSM(pepSeq, &mScan, as<vector<int>>(charges), isoCenter, isoWidth, calScores);
    else
    {
        Rcerr << "Too less peaks or empty Scan" << endl;
        return List();
    }
    std::vector<std::string> ionKindStrs = enumsToStrings(mAnnotator.getIonKinds());
    DataFrame ExpectedBYions = DataFrame::create(Named("mz") = mAnnotator.getExpectedMZs(),
                                                 _("intensity") = mAnnotator.getExpectedIntensities(),
                                                 _("charge") = mAnnotator.getExpectedCharges(),
                                                 _("ionkind") = ionKindStrs,
                                                 _("residuePositions") = mAnnotator.getResiduePositions(),
                                                 _("matchedIndices") = mAnnotator.getMatchedIndices(),
                                                 _("SIPabundances") = mAnnotator.getSIPabundances());
    DataFrame realPeaks = DataFrame::create(Named("mz") = mScan.mz,
                                            _("intensity") = mScan.intensity,
                                            _("charge") = mScan.charge);
    List re = List::create(Named("ExpectedBYions") = ExpectedBYions,
                           _("RealPeaks") = realPeaks,
                           _("MatchedSpectraEntropyScore") = mAnnotator.getMatchedSpectraEntropyScore(),
                           _("MVHscore") = mAnnotator.getMVHscore(),
                           _("WDPscore") = mAnnotator.getWDPscore(),
                           _("XcorrScore") = mAnnotator.getXcorrScore(),
                           _("Peptide") = pepSeq);
    return re;
}

//' scorePSMsimple Score a PSM without isotopic envelope shape modeling
//' @param realMZ mz vector in MS2 scan
//' @param realIntensity intensity vector in MS2 scan
//' @param realCharge charge vector in MS2 scan
//' @param parentCharge precursor charge, 2 for example
//' @param pepSeq a string of peptide
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
//' @return a score of this PSM
//' @examples
//' demo_file <- system.file("extdata", "107728.FT2", package = "Aerith")
//' scan1 <- readOneScanMS2(ftFile = demo_file, scanNumber = 107728)
//' score <- scorePSMsimple(
//'   scan1$peaks$mz,
//'   scan1$peaks$intensity,
//'   scan1$peaks$charge,
//'   2,
//'   "[HSQVFSTAEDNQSAVTIHVLQGER]",
//'   "C13",
//'   0.0107
//' )
//' @export
// [[Rcpp::export]]
double scorePSMsimple(const NumericVector &realMZ, const NumericVector &realIntensity,
                   const NumericVector &realCharge, int parentCharge, const String &pepSeq, const String &Atom, double Prob)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    // compute residue mass and prob again
    computeResidueMassIntensityAgain(Atom, Prob);
    Peptide myPep;
    string sOriginalPeptide, sProteinName;
    int ibeginPos;
    double dPeptideMass;
    char cIdentifyPrefix, cIdentifySuffix, cOriginalPrefix, cOriginalSuffix;
    myPep.setPeptide(pepSeq, sOriginalPeptide, sProteinName,
                     ibeginPos, dPeptideMass, cIdentifyPrefix,
                     cIdentifySuffix, cOriginalPrefix, cOriginalSuffix);
    map<char, double> mapResidueMass;
    // mapResidueMass is not used in this function
    myPep.preprocessing(true, mapResidueMass);
    MS2Scan myScan;
    myScan.isMS2HighRes = true;
    myScan.vdIntensity = as<vector<double>>(realIntensity);
    myScan.vdMZ = as<vector<double>>(realMZ);
    myScan.viCharge = as<vector<int>>(realCharge);
    myScan.iParentChargeState = parentCharge;
    myScan.preprocess();
    myScan.scoreWeightSumHighMS2(&myPep);
    return myScan.vpWeightSumTopPeptides[0]->dScore;
}
