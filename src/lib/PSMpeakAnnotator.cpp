#include "PSMpeakAnnotator.h"

PSMpeakAnnotator::PSMpeakAnnotator(const double tolerancePPM) : tolerancePPM(tolerancePPM) {}

PSMpeakAnnotator::~PSMpeakAnnotator() {}

void PSMpeakAnnotator::removePeaksInIsolationWindow(Scan *mScan, const double isoCenter, const double isoWidth)
{
    auto lower = std::lower_bound(mScan->mz.begin(), mScan->mz.end(), isoCenter + isoWidth / 2.0);
    auto upper = std::upper_bound(mScan->mz.begin(), mScan->mz.end(), isoCenter - isoWidth / 2.0);
    std::vector<size_t> positions;
    for (auto it = upper; it != lower; ++it)
    {
        positions.push_back(std::distance(mScan->mz.begin(), it));
    }
    // Remove elements from vec and vec2 in reverse order to avoid shifting issues
    for (auto it = positions.rbegin(); it != positions.rend(); ++it)
    {
        mScan->mz.erase(mScan->mz.begin() + *it);
        mScan->intensity.erase(mScan->intensity.begin() + *it);
        mScan->charge.erase(mScan->charge.begin() + *it);
    }
}

void PSMpeakAnnotator::generateTheoreticalSpectra(const std::string &peptide)
{
    std::string pepStr = "[" + peptide + "]";
    ProNovoConfig::configIsotopologue.computeProductIon(pepStr, vvdYionMass,
                                                        vvdYionProb, vvdBionMass, vvdBionProb);
}

size_t PSMpeakAnnotator::binarySearchPeak(const Scan *mScan, double Mz)
{
    // init peakIX to invalid value
    size_t peakIX = std::numeric_limits<size_t>::max();
    size_t low = 0;
    size_t high = mScan->mz.size() - 1;
    size_t mid = 0;
    double diff = 0;
    const double mzTolerance = Mz * tolerancePPM / 1e6;
    double currentIntensity = 0;
    while (low <= high)
    {
        mid = (low + high) / 2;
        diff = std::abs(mScan->mz[mid] - Mz);
        if (diff <= mzTolerance)
        {
            // find the peak with highest intensity in tolerance range
            if (mScan->intensity[mid] > currentIntensity)
            {
                peakIX = mid;
                currentIntensity = mScan->intensity[mid];
            }
        }
        if (mScan->mz[mid] < Mz) // Search the right half
        {
            low = mid + 1;
        }
        else // Search the left half
        {
            // in case searched the first peak in mScan
            if (mid == 0)
                break;
            high = mid - 1;
        }
    }
    return peakIX;
}

double PSMpeakAnnotator::calSIPabundancesOfBYion(
    const double baseMass,
    const std::vector<int> &matchedIXs,
    const Scan *mScan, const int SIPelementCount, const int charge)
{
    double baseMZ = baseMass / charge + ProNovoConfig::getProtonMass();
    double sumOfIntensities = 0;
    size_t i = 0;
    int firstDeltaNeutron = 0;
    vector<double> matchedIntensities;
    matchedIntensities.reserve(matchedIXs.size());
    while (i < matchedIXs.size())
    {
        if (matchedIXs[i] != -1)
        {
            sumOfIntensities += mScan->intensity[matchedIXs[i]];
            firstDeltaNeutron = static_cast<int>(std::round((mScan->mz[matchedIXs[i]] - baseMZ) /
                                                            ProNovoConfig::getNeutronMass() * charge));
            matchedIntensities.push_back(mScan->intensity[matchedIXs[i]]);
            break;
        }
        i++;
    }
    i++;
    while (i < matchedIXs.size())
    {
        if (matchedIXs[i] != -1)
        {
            sumOfIntensities += mScan->intensity[matchedIXs[i]];
            matchedIntensities.push_back(mScan->intensity[matchedIXs[i]]);
        }
        i++;
    }
    if (matchedIntensities.size() == 0)
    {
        return -1.0;
    }
    for (auto &intensity : matchedIntensities)
    {
        intensity /= sumOfIntensities;
    }
    double pct = 0.0;
    for (size_t i = 0; i < matchedIntensities.size(); i++)
    {
        pct += matchedIntensities[i] * (i + firstDeltaNeutron);
    }
    pct /= SIPelementCount;
    pct *= 100.;
    return pct;
}

void PSMpeakAnnotator::
    findIsotopicPeaks(
        const std::vector<double> &ionMasses,
        const std::vector<double> &ionIntensities,
        const Scan *mScan,
        const int residuePosition,
        const int charge,
        const PSMpeakAnnotator::ionKind BYkind)
{
    size_t highestExpectedPeakIX = std::distance(ionIntensities.begin(),
                                                 std::max_element(ionIntensities.begin(),
                                                                  ionIntensities.end()));
    std::vector<double> ionMZs(ionMasses.size());
    for (size_t i = 0; i < ionMasses.size(); i++)
    {
        ionMZs[i] = ionMasses[i] / (double)charge + ProNovoConfig::getProtonMass();
    }
    double highestPeakMZ = ionMZs[highestExpectedPeakIX];
    size_t highestObservedPeakIX = binarySearchPeak(mScan, highestPeakMZ);
    std::vector<double> ionCharges(ionMasses.size(), charge);
    std::vector<ionKind> mIonKinds(ionMasses.size(),
                                   (BYkind == B) ? BisotopicPeak : YisotopicPeak);
    mIonKinds[highestExpectedPeakIX] = BYkind;
    std::vector<int> mPositions(ionMasses.size(), residuePosition);
    std::vector<int> mMatchedIndices(ionMasses.size(), -1);

    size_t currentIX = 0;
    size_t foundIX = 0;
    double currentMZ = 0;
    double expectedMZ = 0;
    double maxIntensity = 0;
    double mzTolerance = highestPeakMZ * tolerancePPM / 1e6;
    bool foundIsotopicPeak = false;

    if (highestObservedPeakIX != std::numeric_limits<size_t>::max())
    {
        mMatchedIndices[highestExpectedPeakIX] = (int)highestObservedPeakIX;
        // Go left and right side
        for (int direction : {-1, 1})
        {
            // currentIX is the index of observed isotopic peak
            currentIX = highestObservedPeakIX + direction;
            currentMZ = mScan->mz[currentIX];
            // iso is the index of expected isotopic peak
            for (size_t iso = highestExpectedPeakIX + direction;
                 iso >= 0 && iso < ionMasses.size(); iso += direction)
            {
                foundIsotopicPeak = false;
                maxIntensity = 0;
                expectedMZ = ionMZs[iso];
                while (currentIX >= 0 && currentIX < mScan->mz.size() &&
                       (direction * (currentMZ - expectedMZ) < mzTolerance))
                {
                    // find the matched isotopic peak with max intensity
                    if (std::abs(expectedMZ - currentMZ) < mzTolerance &&
                        mScan->intensity[currentIX] > maxIntensity)
                    {
                        foundIsotopicPeak = true;
                        foundIX = currentIX;
                        maxIntensity = mScan->intensity[currentIX];
                    }
                    currentIX += direction;
                    currentMZ = mScan->mz[currentIX];
                }
                if (foundIsotopicPeak)
                    mMatchedIndices[iso] = foundIX;
                else
                    break;
            }
        }
    }

    double pct = 0;
    if (BYkind == B)
    {
        pct = calSIPabundancesOfBYion(mAveragine.BionsBaseMasses[residuePosition - 1], mMatchedIndices, mScan,
                                      mAveragine.BionsAtomCounts[residuePosition - 1][mAveragine.SIPatomIX], charge);
    }
    else
    {
        pct = calSIPabundancesOfBYion(mAveragine.YionsBaseMasses[residuePosition - 1], mMatchedIndices, mScan,
                                      mAveragine.YionsAtomCounts[residuePosition - 1][mAveragine.SIPatomIX], charge);
    }
    std::vector<double> pcts(mMatchedIndices.size(), -1);
    for (size_t i = 0; i < mMatchedIndices.size(); i++)
    {
        if (mMatchedIndices[i] != -1)
            pcts[i] = pct;
    }

    expectedMZs.insert(expectedMZs.end(), ionMZs.begin(), ionMZs.end());
    expectedIntensities.insert(expectedIntensities.end(), ionIntensities.begin(), ionIntensities.end());
    expectedCharges.insert(expectedCharges.end(), ionCharges.begin(), ionCharges.end());
    ionKinds.insert(ionKinds.end(), mIonKinds.begin(), mIonKinds.end());
    residuePositions.insert(residuePositions.end(), mPositions.begin(), mPositions.end());
    matchedIndices.insert(matchedIndices.end(), mMatchedIndices.begin(), mMatchedIndices.end());
    SIPabundances.insert(SIPabundances.end(), pcts.begin(), pcts.end());
}

void PSMpeakAnnotator::matchIsotopicEnvelopes(Scan *mRealScan, const int charge)
{
    realScan = mRealScan;
    for (size_t i = 0; i < vvdBionMass.size(); i++)
    {
        findIsotopicPeaks(vvdBionMass[i], vvdBionProb[i], mRealScan, i + 1, charge, B);
    }
    for (size_t i = 0; i < vvdYionMass.size(); i++)
    {
        findIsotopicPeaks(vvdYionMass[i], vvdYionProb[i], mRealScan, i + 1, charge, Y);
    }
}

void PSMpeakAnnotator::calMatchedSpectraEntropy()
{
    double totalIntensity = 0;
    double entropy = 0;
    for (size_t i = 0; i < matchedIndices.size(); i++)
    {
        if (matchedIndices[i] != -1)
        {
            totalIntensity += realScan->intensity[matchedIndices[i]];
        }
    }
    for (size_t i = 0; i < matchedIndices.size(); i++)
    {
        if (matchedIndices[i] != -1)
        {
            double prob = realScan->intensity[matchedIndices[i]] / totalIntensity;
            entropy += -prob * std::log(prob);
        }
    }
    matchedSpectraEntropy = entropy;
}

// Compute log(binomial(n, k)) using lgamma for numerical stability.
double PSMpeakAnnotator::logBinom(int n, int k) {
    if (k < 0 || k > n)
        return -INFINITY;
    return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

// Compute one term of the hypergeometric probability:
// term = (n_theo choose j) * ((N - n_theo) choose (n_obs - j)) / (N choose n_obs)
double PSMpeakAnnotator::hypergeomProbability(int n_theo, int n_obs, int j, int N) {
    if (j > n_theo || (n_obs - j) > (N - n_theo))
        return 0.0;
    double logTerm = logBinom(n_theo, j)
                   + logBinom(N - n_theo, n_obs - j)
                   - logBinom(N, n_obs);
    return exp(logTerm);
}

void PSMpeakAnnotator::calMVHscore()
{
    double maxMz = *std::max_element(realScan->mz.begin(), realScan->mz.end());
    double minMz = *std::min_element(realScan->mz.begin(), realScan->mz.end());
    // total number of tolerance window
    int N = (maxMz - minMz) / ((maxMz + minMz) / 2 * tolerancePPM / 1e6 * 2);
    int n_theo = 0, n_match = 0, n_obs = realScan->mz.size();
    // consider top 3 peaks of isotopic envelope
    const int topN = 3;
    size_t ix = 0, nISO = 0;
    int lastResiduePosition = 1;
    while (ix <- expectedMZs.size())
    {
        if (residuePositions[ix] != lastResiduePosition)
        {
            lastResiduePosition = residuePositions[ix];
            nISO = 0;
        }
        if (expectedMZs[ix] >= minMz && expectedMZs[ix] <= maxMz)
        {
            if (ionKinds[ix] == B || ionKinds[ix] == Y)
                n_theo = n_theo + topN;
            if (matchedIndices[ix] != -1 && nISO <= topN)
            {
                n_match ++;                   
                nISO ++; 
            }
        }
        ix ++;
    }    
    double pValue = 0.0;
    for (int j = n_match; j <= n_theo; ++j) {
        pValue += hypergeomProbability(n_theo, n_obs, j, N);
    }
    // Prevent taking log of zero.
    if (pValue < 1e-300)
        pValue = 1e-300;
    MVHscore = -std::log10(pValue);
}

void PSMpeakAnnotator::calWDPscore()
{
}

void PSMpeakAnnotator::calXcorrScore()
{
}

void PSMpeakAnnotator::scorePSM()
{
}


void PSMpeakAnnotator::analyzePSM(const std::string &peptide, Scan *realScan,
                                  const std::vector<int> &charges,
                                  const double isoCenter, const double isoWidth)
{
    if (isoCenter != 0 && isoWidth != 0)
        removePeaksInIsolationWindow(realScan, isoCenter, isoWidth);
    generateTheoreticalSpectra(peptide);
    mAveragine.calBYionBaseMasses(peptide);
    for (int charge : charges)
    {
        matchIsotopicEnvelopes(realScan, charge);
    }
    calMatchedSpectraEntropy();
}

double PSMpeakAnnotator::getScore()
{
    return Score;
}

double PSMpeakAnnotator::getMatchedSpectraEntropy()
{
    return matchedSpectraEntropy;
}

std::vector<int> PSMpeakAnnotator::getMatchedIndices()
{
    return matchedIndices;
}

std::vector<PSMpeakAnnotator::ionKind> PSMpeakAnnotator::getIonKinds()
{
    return ionKinds;
}

std::vector<double> PSMpeakAnnotator::getSIPabundances()
{
    return SIPabundances;
}

std::vector<int> PSMpeakAnnotator::getResiduePositions()
{
    return residuePositions;
};

std::vector<double> PSMpeakAnnotator::getExpectedMZs()
{
    return expectedMZs;
}

std::vector<double> PSMpeakAnnotator::getExpectedIntensities()
{
    return expectedIntensities;
}

std::vector<int> PSMpeakAnnotator::getExpectedCharges()
{
    return expectedCharges;
}