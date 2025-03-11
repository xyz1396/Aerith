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

void PSMpeakAnnotator::calMatchedSpectraEntropyScore()
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
    matchedSpectraEntropyScore = entropy;
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
    while (ix < expectedMZs.size())
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

SparseBins PSMpeakAnnotator::sparseBinSpectrum(const std::vector<double> &mz,
    const vector<double> &intensity, double binWidth,
    double minMz, double maxMz) {
        SparseBins bins;
        bins.reserve(mz.size());
        for (size_t i = 0; i < mz.size(); ++i) {
            double m = mz[i];
            double inten = intensity[i];
            if (m < minMz || m > maxMz)
                continue;
            int binIndex = static_cast<int>(std::floor((m - minMz) / binWidth));
            // Update bin to hold the max intensity.
            auto it = bins.find(binIndex);
            if (it == bins.end())
                bins[binIndex] = inten;
            else
                it->second = max(it->second, inten);
        }
        return bins;
    }

// Compute correlation at a given shift using sparse bins.
double PSMpeakAnnotator::sparseCorrelationAtShift(const SparseBins &theoBins, const SparseBins &obsBins, int shift) {
    double sum = 0.0;
    for (const auto &p : theoBins) {
        int binIndex = p.first;
        int shiftIndex = binIndex + shift;
        auto it = obsBins.find(shiftIndex);
        if (it != obsBins.end()) {
            sum += p.second * it->second;
        }
    }
    return sum;
}

void PSMpeakAnnotator::selectTopPeaks(std::vector<double> &theo_mz,
    std::vector<double> &theo_intensity, int topN, double intensityThreshold)
{
    auto processIons = [&](const std::vector<std::vector<double>> &ionMasses,
                             const std::vector<std::vector<double>> &ionProbs)
    {
        for (size_t i = 0; i < ionMasses.size(); i++)
        {
            if (ionMasses[i].size() != ionProbs[i].size())
            {   
                std::cerr << "Error: ionMasses and ionProbs have different sizes." << std::endl;
                continue;
            }
            double maxInten = *std::max_element(ionProbs[i].begin(), ionProbs[i].end());
            // Make a copy of the probability vector so we can modify it.
            std::vector<double> probCopy = ionProbs[i];
            for (int n = 0; n < topN; n++)
            {
                auto maxIt = std::max_element(probCopy.begin(), probCopy.end());
                if (maxIt == probCopy.end() || *maxIt < 0)
                    break; // No valid peak remaining.
                size_t index = std::distance(probCopy.begin(), maxIt);
                double normalizedInten = ionProbs[i][index] / maxInten;
                if (normalizedInten < intensityThreshold)
                    break; // No peak with intensity above threshold.
                theo_intensity.push_back(normalizedInten);
                theo_mz.push_back(ionMasses[i][index]);
                // Mark this peak so it won't be chosen again.
                *maxIt = -std::numeric_limits<double>::infinity();
            }
        }
    };
    processIons(vvdBionMass, vvdBionProb);
    processIons(vvdYionMass, vvdYionProb);
}

void PSMpeakAnnotator::calXcorrScore()
{
    std::vector<double> obs_mz = realScan->mz;
    // convert mz to single charge and remove proton mass
    for (size_t i = 0; i < obs_mz.size(); i++)
    {
        if (realScan->charge[i] == 0)
            obs_mz[i] = obs_mz[i] - ProNovoConfig::getProtonMass();
        else
            obs_mz[i] = obs_mz[i] * realScan->charge[i] - realScan->charge[i] * ProNovoConfig::getProtonMass();
    }
    std::vector<double> obs_intensity = realScan->intensity;
    double maxInten = *std::max_element(obs_intensity.begin(), obs_intensity.end());
    for (size_t i = 0; i < obs_intensity.size(); i++)
    {
        obs_intensity[i] /= maxInten;
    }

    std::vector<double> theo_mz, theo_intensity;
    int topN = 3;
    theo_mz.reserve(vvdBionMass.size() * topN + vvdYionMass.size() * topN);
    theo_intensity.reserve(vvdBionMass.size() * topN + vvdYionMass.size() * topN);
    double intensityThreshold = 0.01;
    selectTopPeaks(theo_mz, theo_intensity, topN, intensityThreshold);

    auto obsRange = std::minmax_element(obs_mz.begin(), obs_mz.end());
    double minMz = *obsRange.first;
    double maxMz = *obsRange.second;
    double binWidth = tolerancePPM / 1e6 * (minMz + maxMz) / 2;
    int shiftWindow = 75;
    // Add margins to account for shifts.
    minMz -= shiftWindow * binWidth;
    minMz = std::max(minMz, 0.0);
    maxMz += shiftWindow * binWidth;
  
    // Create sparse bins for both the theoretical and observed spectra.
    SparseBins theoBins = sparseBinSpectrum(theo_mz, theo_intensity, binWidth, minMz, maxMz);
    SparseBins obsBins  = sparseBinSpectrum(obs_mz, obs_intensity, binWidth, minMz, maxMz);
  
    // Correlation at zero shift.
    double cc0 = sparseCorrelationAtShift(theoBins, obsBins, 0);
  
    // Compute the average correlation over shifted windows (excluding zero).
    double ccSum = 0.0;
    int cnt = 0;
    for (int shift = -shiftWindow; shift <= shiftWindow; ++shift) {
        if (shift == 0)
            continue;
        ccSum += sparseCorrelationAtShift(theoBins, obsBins, shift);
        cnt++;
    }
    double avgShift = (cnt > 0) ? ccSum / cnt : 0.0;
  
    // XCorr is defined as the zero-shift correlation minus the average shifted correlation.
    XcorrScore = cc0 - avgShift;
    if (XcorrScore < 0)
        XcorrScore = 0;
}

void PSMpeakAnnotator::scorePSM()
{
}


void PSMpeakAnnotator::analyzePSM(const std::string &peptide, Scan *realScan,
                                  const std::vector<int> &charges,
                                  const double isoCenter, const double isoWidth, const bool calScores)
{
    if (isoCenter != 0 && isoWidth != 0)
        removePeaksInIsolationWindow(realScan, isoCenter, isoWidth);
    generateTheoreticalSpectra(peptide);
    mAveragine.calBYionBaseMasses(peptide);
    for (int charge : charges)
    {
        matchIsotopicEnvelopes(realScan, charge);
    }
    if (calScores)
    {
        calMVHscore();
        calWDPscore();
        calXcorrScore();
        calMatchedSpectraEntropyScore();
    }
}

double PSMpeakAnnotator::getScore()
{
    return Score;
}

double PSMpeakAnnotator::getWDPscore()
{
    return WDPscore;
}

double PSMpeakAnnotator::getMVHscore()
{
    return MVHscore;
}

double PSMpeakAnnotator::getXcorrScore()
{
    return XcorrScore;
}

double PSMpeakAnnotator::getMatchedSpectraEntropyScore()
{
    return matchedSpectraEntropyScore;
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