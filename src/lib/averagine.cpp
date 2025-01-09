#include "averagine.h"

averagine::averagine(const int minPepLen, const int maxPepLen)
    : minPepLen(minPepLen), maxPepLen(maxPepLen), pepLenRange(maxPepLen - minPepLen + 1)
{
    // get averagine's Atom content and mass distribution
    averagineAtomCounts = ProNovoConfig::configIsotopologue.mResidueAtomicComposition[averagineResidue];
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(averagineAtomCounts,
                                                                  averagineSIPdistribution);
    diffAtomCounts.resize(averagineAtomCounts.size());
    C12Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vMass[0];
    C13Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vMass[1];
    C13Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[1];
    H1Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vMass[0];
    H2Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vMass[1];
    H2Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vProb[1];
    O16Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[0];
    O17Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[1];
    O17Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vProb[1];
    O18Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[2];
    O18Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vProb[2];
    N14Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vMass[0];
    N15Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vMass[1];
    N15Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[1];
    PfakeLowMass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vMass[0];
    PfakeMass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vMass[1];
    PfakeAbundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vProb[1];
    S32Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[0];
    S33Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[1];
    S33Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[1];
    S34Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[2];
    S34Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[2];
    S36Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[3];
    S36Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[3];
    adjustEstimatePrecursorMassbyNP();
}

averagine::averagine()
{
    diffAtomCounts.resize(6);
    C12Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vMass[0];
    C13Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vMass[1];
    C13Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[1];
    H1Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vMass[0];
    H2Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vMass[1];
    H2Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[1].vProb[1];
    O16Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[0];
    O17Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[1];
    O17Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vProb[1];
    O18Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vMass[2];
    O18Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[2].vProb[2];
    N14Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vMass[0];
    N15Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vMass[1];
    N15Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[1];
    PfakeLowMass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vMass[0];
    PfakeMass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vMass[1];
    PfakeAbundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[4].vProb[1];
    S32Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[0];
    S33Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[1];
    S33Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[1];
    S34Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[2];
    S34Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[2];
    S36Mass = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vMass[3];
    S36Abundance = &ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[5].vProb[3];
    adjustEstimatePrecursorMassbyNP();
}

averagine::~averagine()
{
}

void averagine::adjustEstimatePrecursorMassbyNP()
{
    // for SIP atom C13
    SIPatomIX = 0;
    estimatePrecursorMassbyNP = [&](double baseMass, int atomCount, double neutronMass)
    { return (baseMass + std::round(atomCount * (*C13Abundance)) * neutronMass); };
    // for SIP atom other than C13
    // SIP element shif + C13 element shif because C13 is abundant in nature (1%)
    if (ProNovoConfig::getSetSIPelement() == "H")
    {
        SIPatomIX = 1;
        estimatePrecursorMassbyNP = [&](double baseMass, int atomCount, double neutronMass)
        { return (baseMass +
                  std::round(atomCount * (*H2Abundance) + pepAtomCounts[0] * (*C13Abundance)) * neutronMass); };
    }
    else if (ProNovoConfig::getSetSIPelement() == "O")
    {
        SIPatomIX = 2;
        estimatePrecursorMassbyNP = [&](double baseMass, int atomCount, double neutronMass)
        { return (baseMass +
                  std::round(2 * atomCount * (*O18Abundance) + pepAtomCounts[0] * (*C13Abundance)) * neutronMass); };
    }
    else if (ProNovoConfig::getSetSIPelement() == "N")
    {
        SIPatomIX = 3;
        estimatePrecursorMassbyNP = [&](double baseMass, int atomCount, double neutronMass)
        { return (baseMass +
                  std::round(atomCount * (*N15Abundance) + pepAtomCounts[0] * (*C13Abundance)) * neutronMass); };
    }
    else if (ProNovoConfig::getSetSIPelement() == "S")
    {
        SIPatomIX = 5;
        estimatePrecursorMassbyNP = [&](double baseMass, int atomCount, double neutronMass)
        { return (baseMass +
                  std::round(2 * atomCount * (*S34Abundance) + pepAtomCounts[0] * (*C13Abundance)) * neutronMass); };
    }
}

void averagine::changeAtomSIPabundance(const char SIPatom, const double pct)
{
    bool SIPatomFound = false;
    for (size_t i = 0; i < SIPatoms.size(); i++)
    {
        if (SIPatoms[i] == SIPatom)
        {
            SIPatomIX = i;
            SIPatomFound = true;
            break;
        }
    }
    if (SIPatomFound)
    {
        ProNovoConfig::getSetSIPelement() = SIPatom;
        adjustEstimatePrecursorMassbyNP();
        ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[SIPatomIX].vProb[0] =
            1.0 - pct;
        // for O18 and S34
        if (SIPatom == 'O' || SIPatom == 'S')
        {            
            ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[SIPatomIX].vProb[2] =
            pct;
        }
        // for C13, H2, N15
        else
        {
            ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[SIPatomIX].vProb[1] =
            pct;
        }
        ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
            ProNovoConfig::configIsotopologue.mResidueAtomicComposition[averagineResidue],
            averagineSIPdistribution);
    }
    else
    {
        cout << SIPatom << "element is not supported!" << endl;
    }
}

void averagine::calAveraginePepAtomCounts()
{
    vector<int> averaginePepAtomCounts;
    averaginePepAtomCountss.reserve(pepLenRange);
    averaginePepAtomCounts.reserve(averagineAtomCounts.size());
    for (int pepLen = minPepLen; pepLen <= maxPepLen; pepLen++)
    {
        averaginePepAtomCounts.clear();
        for (size_t j = 0; j < averagineAtomCounts.size(); j++)
            averaginePepAtomCounts.push_back(averagineAtomCounts[j] * pepLen);
        averaginePepAtomCountss.push_back(averaginePepAtomCounts);
    }
}

vector<int> *averagine::getAveraginePepAtomCounts(const int pepLen)
{
    return &averaginePepAtomCountss[pepLen - minPepLen];
}

void averagine::calAveraginePepSIPdistributions()
{
    IsotopeDistribution NtermSIPdistribution =
        ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution["Nterm"];
    IsotopeDistribution CtermSIPdistribution =
        ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution["Cterm"];
    averaginePepSIPdistributions.reserve(pepLenRange);
    IsotopeDistribution tempSIPdistribution = NtermSIPdistribution;
    ProNovoConfig::configIsotopologue.sum(CtermSIPdistribution, tempSIPdistribution);
    for (int i = 0; i < minPepLen; i++)
    {
        ProNovoConfig::configIsotopologue.sum(averagineSIPdistribution, tempSIPdistribution);
    }
    for (int i = 0; i <= pepLenRange; i++)
    {
        ProNovoConfig::configIsotopologue.sum(averagineSIPdistribution, tempSIPdistribution);
        averaginePepSIPdistributions.push_back(tempSIPdistribution);
    }
}

IsotopeDistribution *averagine::getAveraginePepSIPdistribution(const int pepLen)
{
    return &averaginePepSIPdistributions[pepLen - minPepLen];
}

void averagine::calPepAtomCounts(const string &pepSeq)
{
    std::fill(pepAtomCounts.begin(), pepAtomCounts.end(), 0);
    vector<int> *residueAtomCounts;
    residueAtomCounts = &ProNovoConfig::configIsotopologue
                             .mResidueAtomicComposition["Nterm"];
    for (size_t j = 0; j < residueAtomCounts->size(); j++)
        pepAtomCounts[j] += residueAtomCounts->operator[](j);
    residueAtomCounts = &ProNovoConfig::configIsotopologue
                             .mResidueAtomicComposition["Cterm"];
    for (size_t j = 0; j < residueAtomCounts->size(); j++)
        pepAtomCounts[j] += residueAtomCounts->operator[](j);
    for (size_t i = 0; i < pepSeq.size(); i++)
    {
        if (ProNovoConfig::configIsotopologue
                .mResidueAtomicComposition.find(pepSeq.substr(i, 1)) == ProNovoConfig::configIsotopologue
                                                                            .mResidueAtomicComposition.end())
        {
            cerr << "ERROR: cannot find " << pepSeq.substr(i, 1) << " residue or PTM in the config file." << endl;
            return;
        }
        residueAtomCounts = &ProNovoConfig::configIsotopologue
                                 .mResidueAtomicComposition[pepSeq.substr(i, 1)];
        for (size_t j = 0; j < residueAtomCounts->size(); j++)
            pepAtomCounts[j] += residueAtomCounts->operator[](j);
    }
}

void averagine::calBYionsAtomCounts(const string &pepSeq)
{
    int peptideLength = pepSeq.size();
    BionsAtomCounts.clear();
    YionsAtomCounts.clear();
    BionsAtomCounts.reserve(peptideLength);
    YionsAtomCounts.reserve(peptideLength);
    // one hydrogen atom shifts from b ion to its conuterpart when cleavage
    std::array<int, 6> bIonCounts = {0, 0, 0, 0, 0, 0};
    // one hydrogen atom shifts to y ion from its conuterpart when cleavage
    std::array<int, 6> yIonCounts = {0, 2, 1, 0, 0, 0};
    vector<int> *residueAtomCounts;
    // Calculate atom counts for B ions
    for (int i = 0; i < peptideLength; i++)
    {
        if (ProNovoConfig::configIsotopologue
                .mResidueAtomicComposition.find(pepSeq.substr(i, 1)) == ProNovoConfig::configIsotopologue
                                                                            .mResidueAtomicComposition.end())
        {
            cerr << "ERROR: cannot find " << pepSeq.substr(i, 1) << " residue or PTM in the config file." << endl;
            return;
        }
        if (!isalpha(pepSeq[i]))
        {
            string ptmSymbol = pepSeq.substr(i, 1);
            residueAtomCounts =
                &ProNovoConfig::configIsotopologue
                     .mResidueAtomicComposition[ptmSymbol];
            for (size_t k = 0; k < residueAtomCounts->size(); ++k)
                bIonCounts[k] += residueAtomCounts->operator[](k);
            // if the PTM is not at Nterm
            if (BionsAtomCounts.size() > 0)
                BionsAtomCounts.back() = bIonCounts;
            continue; // Move to the next character
        }
        const string residue = pepSeq.substr(i, 1);
        residueAtomCounts =
            &ProNovoConfig::configIsotopologue.mResidueAtomicComposition[residue];
        for (size_t k = 0; k < residueAtomCounts->size(); ++k)
            bIonCounts[k] += residueAtomCounts->operator[](k);
        BionsAtomCounts.push_back(bIonCounts);
    }
    // remove the last element, which is the full peptide without cterm
    BionsAtomCounts.pop_back();

    // Calculate atom counts for Y ions
    for (int i = peptideLength - 1; i >= 0; i--)
    {
        if (!isalpha(pepSeq[i]))
        {
            string ptmSymbol = pepSeq.substr(i, 1);
            residueAtomCounts =
                &ProNovoConfig::configIsotopologue
                     .mResidueAtomicComposition[ptmSymbol];
            for (size_t k = 0; k < residueAtomCounts->size(); ++k)
                yIonCounts[k] += residueAtomCounts->operator[](k);
            continue; // Move to the next character
        }
        const string residue = pepSeq.substr(i, 1);
        residueAtomCounts =
            &ProNovoConfig::configIsotopologue.mResidueAtomicComposition[residue];
        // Sum the counts
        for (size_t k = 0; k < residueAtomCounts->size(); ++k)
            yIonCounts[k] += residueAtomCounts->operator[](k);
        YionsAtomCounts.push_back(yIonCounts);
    }
    // remove the last element, which is the full peptide without nterm
    YionsAtomCounts.pop_back();
}

double averagine::weighted_mean(const std::vector<double> &values, const std::vector<double> &weights)
{
    double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
    double sum_product = std::inner_product(values.begin(), values.end(), weights.begin(), 0.0);
    return (sum_product / sum_weights);
}

double averagine::calNetronMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    std::vector<double> weights = {*C13Abundance, *H2Abundance, *O17Abundance + *O18Abundance,
                                   *N15Abundance, *PfakeAbundance, *S33Abundance + *S34Abundance + *S36Abundance};
    for (size_t i = 0; i < pepAtomCounts.size(); i++)
    {
        weights[i] *= pepAtomCounts[i];
    }
    double OdeltaMass = weighted_mean({*O17Mass - *O16Mass, *O18Mass - *O17Mass}, {*O17Abundance, *O18Abundance});
    double SdeltaMass = weighted_mean({*S33Mass - *S32Mass, *S34Mass - *S33Mass, (*S36Mass - *S34Mass) / 2.0},
                                      {*S33Abundance, *S34Abundance, *S36Abundance * 2});
    double neutronMass = weighted_mean({*C13Mass - *C12Mass, *H2Mass - *H1Mass, OdeltaMass,
                                        *N15Mass - *N14Mass, *PfakeMass - *PfakeLowMass, SdeltaMass},
                                       weights);
    return neutronMass;
}

double averagine::calPrecursorBaseMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    double baseMass = pepAtomCounts[0] * (*C12Mass) + pepAtomCounts[1] * (*H1Mass) + pepAtomCounts[2] * (*O16Mass) +
                      pepAtomCounts[3] * (*N14Mass) + pepAtomCounts[4] * (*PfakeLowMass) +
                      pepAtomCounts[5] * (*S32Mass);
    return baseMass;
}

void averagine::calBYionBaseMasses(const string &pepSeq)
{
    calBYionsAtomCounts(pepSeq);
    BionsBaseMasses.clear();
    YionsBaseMasses.clear();
    BionsBaseMasses.reserve(BionsAtomCounts.size());
    YionsBaseMasses.reserve(YionsAtomCounts.size());
    for (size_t i = 0; i < BionsAtomCounts.size(); i++)
    {
        BionsBaseMasses.push_back(
            BionsAtomCounts[i][0] * (*C12Mass) + BionsAtomCounts[i][1] * (*H1Mass) +
            BionsAtomCounts[i][2] * (*O16Mass) + BionsAtomCounts[i][3] * (*N14Mass) +
            BionsAtomCounts[i][4] * (*PfakeLowMass) +
            BionsAtomCounts[i][5] * (*S32Mass));
    }
    for (size_t i = 0; i < YionsAtomCounts.size(); i++)
    {
        YionsBaseMasses.push_back(
            YionsAtomCounts[i][0] * (*C12Mass) + YionsAtomCounts[i][1] * (*H1Mass) +
            YionsAtomCounts[i][2] * (*O16Mass) + YionsAtomCounts[i][3] * (*N14Mass) +
            YionsAtomCounts[i][4] * (*PfakeLowMass) +
            YionsAtomCounts[i][5] * (*S32Mass));
    }
}

double averagine::calPrecursorMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    double neutronMass = calNetronMass(pepSeq);
    double baseMass = pepAtomCounts[0] * (*C12Mass) + pepAtomCounts[1] * (*H1Mass) + pepAtomCounts[2] * (*O16Mass) +
                      pepAtomCounts[3] * (*N14Mass) + pepAtomCounts[4] * (*PfakeLowMass) +
                      pepAtomCounts[5] * (*S32Mass);
    double precursorMass = estimatePrecursorMassbyNP(baseMass, pepAtomCounts[SIPatomIX], neutronMass);
    return (precursorMass);
}

void averagine::calDiffAtomCounts(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    diffAtomCounts = vector<int>(pepAtomCounts.begin(), pepAtomCounts.end());
    vector<int> *averaginePepAtomCountPtr = getAveraginePepAtomCounts(pepSeq.length());
    for (size_t k = 0; k < diffAtomCounts.size(); k++)
        diffAtomCounts[k] -= (*averaginePepAtomCountPtr)[k];
}

void averagine::calPrecursorIsotopeDistribution(const string &pepSeq, IsotopeDistribution &tempSIPdistribution)
{
    calDiffAtomCounts(pepSeq);
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(diffAtomCounts,
                                                                  tempSIPdistribution);
    ProNovoConfig::configIsotopologue.sum(*(getAveraginePepSIPdistribution(pepSeq.length())),
                                          tempSIPdistribution);
}