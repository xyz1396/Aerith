#include "PSMfeatureExtractor.h"

PSMfeatureExtractor::PSMfeatureExtractor()
{
}

size_t PSMfeatureExtractor::binarySearchPeak(const Scan *mScan, double Mz, int charge)
{
    // init peakIX to invalid value
    size_t peakIX = std::numeric_limits<size_t>::max();
    size_t low = 0;
    size_t high = mScan->mz.size() - 1;
    size_t mid = 0;
    double diff = 0;
    const double mzTolerance = 0.01 / charge;
    double currentIntensity = 0;
    while (low <= high)
    {
        mid = (low + high) / 2;
        diff = std::abs(mScan->mz[mid] - Mz);
        if (diff <= mzTolerance)
        {
            // find the peak with highest intensity in tolerance range
            if (mScan->charge[mid] == charge && mScan->intensity[mid] > currentIntensity)
            {
                peakIX = mid;
                currentIntensity = mScan->intensity[mid];
                // in case matched the first peak in mScan
                if (mid == 0)
                    break;
            }
        }
        if (mScan->mz[mid] < Mz) // Search the right half
        {
            low = mid + 1;
        }
        else // Search the left half
        {
            // in case matched the first peak in mScan
            if (mid == 0)
                high = 0;
            else
                high = mid - 1;
        }
    }
    return peakIX;
}

void PSMfeatureExtractor::filterIsotopicPeaks(std::vector<isotopicPeak> &isotopicPeaks,
                                              const double calculatedPrecursorMZ)
{
    double intensityDiff = 0, lastIntensityDiff = 1.0;
    vertexIXs.clear();
    for (size_t i = 1; i < isotopicPeaks.size(); i++)
    {
        intensityDiff = isotopicPeaks[i].intensity - isotopicPeaks[i - 1].intensity;
        if (lastIntensityDiff > 0 && intensityDiff < 0)
            vertexIXs.push_back(i - 1);
        lastIntensityDiff = intensityDiff;
    }
    if (lastIntensityDiff > 0)
        vertexIXs.push_back(isotopicPeaks.size() - 1);
    // find closest vertexIX to calculatedPrecursorMZ
    int closestVertexIX = 0;
    double closestVertexMZdiff = std::numeric_limits<double>::max(), MZdiff = 0;
    for (size_t i = 0; i < vertexIXs.size(); i++)
    {
        MZdiff = std::abs(isotopicPeaks[vertexIXs[i]].mz - calculatedPrecursorMZ);
        if (MZdiff < closestVertexMZdiff)
        {
            closestVertexMZdiff = MZdiff;
            closestVertexIX = vertexIXs[i];
        }
    }
    int currentPeakIX = closestVertexIX - 1;
    while (currentPeakIX >= 0)
    {
        intensityDiff = isotopicPeaks[currentPeakIX].intensity - isotopicPeaks[currentPeakIX + 1].intensity;
        if (intensityDiff > 0)
        {
            isotopicPeaks.erase(isotopicPeaks.begin(), isotopicPeaks.begin() + currentPeakIX + 1);
            // shift closestVertexIX because of erasing
            closestVertexIX = closestVertexIX - currentPeakIX - 1;
            break;
        }
        currentPeakIX--;
    }
    currentPeakIX = closestVertexIX + 1;
    while (currentPeakIX < (int)isotopicPeaks.size())
    {
        intensityDiff = isotopicPeaks[currentPeakIX].intensity - isotopicPeaks[currentPeakIX - 1].intensity;
        if (intensityDiff > 0)
        {
            isotopicPeaks.erase(isotopicPeaks.begin() + currentPeakIX, isotopicPeaks.end());
            break;
        }
        currentPeakIX++;
    }
}

std::vector<isotopicPeak> PSMfeatureExtractor::
    findIsotopicPeaks(const size_t MS1ScanNumber,
                      const int precursorCharge,
                      const double observedPrecursorMass,
                      const double calculatedPrecursorMass)
{
    std::vector<isotopicPeak> isotopicPeaks = {};
    double observedPrecursorMZ = observedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    double calculatedPrecursorMZ = calculatedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    Scan *MS1Scan = scanNumerFT1ScanMap[MS1ScanNumber];
    size_t peakIX = binarySearchPeak(MS1Scan, observedPrecursorMZ, precursorCharge);
    size_t currentIX = 0;
    size_t foundIX = 0;
    double currentMass = 0;
    double expectedMass = 0;
    double maxIntensity = 0;
    double massTolerance = 0.01;
    bool foundIsotopicPeak = false;
    if (peakIX != std::numeric_limits<size_t>::max())
    {
        // try NisotopicPeak isotopic peaks each side
        int tryISO = NisotopicPeak;
        isotopicPeaks.reserve(tryISO);
        isotopicPeaks.push_back({MS1Scan->mz[peakIX],
                                 precursorCharge, MS1Scan->intensity[peakIX]});
        // Go left and right side
        for (int direction : {-1, 1})
        {
            currentIX = peakIX + direction;
            currentMass = MS1Scan->mz[currentIX] * precursorCharge;
            for (int iso = 1; iso <= tryISO; iso++)
            {
                foundIsotopicPeak = false;
                maxIntensity = 0;
                expectedMass = MS1Scan->mz[peakIX] * precursorCharge +
                               direction * iso * ProNovoConfig::getNeutronMass();
                while (currentIX >= 0 && currentIX < MS1Scan->mz.size() &&
                       (direction * (currentMass - expectedMass) < massTolerance))
                {
                    // find the matched isotopic peak with max intensity
                    if (std::abs(expectedMass - currentMass) < massTolerance &&
                        // MS1Scan->charge[currentIX] == precursorCharge &&
                        MS1Scan->intensity[currentIX] > maxIntensity)
                    {
                        foundIsotopicPeak = true;
                        foundIX = currentIX;
                        maxIntensity = MS1Scan->intensity[currentIX];
                    }
                    currentIX += direction;
                    currentMass = MS1Scan->mz[currentIX] * precursorCharge;
                }
                if (foundIsotopicPeak)
                    isotopicPeaks.push_back({MS1Scan->mz[foundIX],
                                             //  precursorCharge,
                                             MS1Scan->charge[foundIX],
                                             MS1Scan->intensity[foundIX]});
                else
                    break;
            }
        }
        std::sort(isotopicPeaks.begin(), isotopicPeaks.end(), [](const isotopicPeak &a, const isotopicPeak &b)
                  { return a.mz < b.mz; });
    }
    if (isotopicPeaks.size() > 2)
    {
        filterIsotopicPeaks(isotopicPeaks, calculatedPrecursorMZ);
    }
    return isotopicPeaks;
}

double PSMfeatureExtractor::getSIPelementAbundanceFromMS1(const std::string &peptideSeq,
                                                          const std::vector<isotopicPeak> &isotopicPeaks,
                                                          const int precursorCharge)
{
    std::size_t start = peptideSeq.find_first_of('[');
    std::size_t end = peptideSeq.find_last_of(']');
    std::string seq = peptideSeq.substr(start + 1, end - start - 1);
    double baseMass = mAveragine.calPrecursorBaseMass(seq);
    int charge = precursorCharge;
    double baseMZ = baseMass / charge + ProNovoConfig::getProtonMass();
    double MZthreshold = baseMZ - 0.5 / charge;
    std::vector<double> usefulIsotopicPeakIntensity;
    usefulIsotopicPeakIntensity.reserve(NisotopicPeak);
    for (const auto &peak : isotopicPeaks)
    {
        if (peak.mz > MZthreshold)
        {
            usefulIsotopicPeakIntensity.push_back(peak.intensity);
        }
    }
    int firstUsefulIsotopicPeakIX = isotopicPeaks.size() - usefulIsotopicPeakIntensity.size();
    double firstUsefulIsotopicPeakMZ = isotopicPeaks[firstUsefulIsotopicPeakIX].mz;
    int firstDeltaNeutron = static_cast<int>(std::round((firstUsefulIsotopicPeakMZ - baseMZ) /
                                                        ProNovoConfig::getNeutronMass() * charge));
    double sumOfIntensities = std::accumulate(usefulIsotopicPeakIntensity.begin(),
                                              usefulIsotopicPeakIntensity.end(), 0.0);
    for (auto &intensity : usefulIsotopicPeakIntensity)
    {
        intensity /= sumOfIntensities;
    }
    double pct = 0.0;
    for (size_t i = 0; i < usefulIsotopicPeakIntensity.size(); i++)
    {
        pct += usefulIsotopicPeakIntensity[i] * (i + firstDeltaNeutron);
    }
    double atomCnumber = mAveragine.pepAtomCounts[0];
    pct /= atomCnumber;
    return pct;
}

void PSMfeatureExtractor::loadFT1(const std::string &FTfileBasePath)
{
    ftFileReader FT1FileReader(FTfileBasePath + ".FT1");
    FT1FileReader.readAllScan();
    FT1Scans = std::move(FT1FileReader.Scans);
    FT1FileReader.Scans = {};
    scanNumerFT1ScanMap.reserve(FT1Scans.size());
    for (size_t i = 0; i < FT1Scans.size(); i++)
    {
        scanNumerFT1ScanMap.insert({FT1Scans[i].scanNumber, FT1Scans.data() + i});
    }
}

void PSMfeatureExtractor::loadFT1FT2fileParallel(const std::string &FTfileBasePath)
{
#pragma omp parallel sections
    {
#pragma omp section
        {
            ftFileReader FT2FileReader(FTfileBasePath + ".FT2");
            FT2FileReader.readAllScan();
            FT2Scans = std::move(FT2FileReader.Scans);
            FT2FileReader.Scans = {};
            scanNumerFT2ScanMap.reserve(FT2Scans.size());
            for (size_t i = 0; i < FT2Scans.size(); i++)
            {
                scanNumerFT2ScanMap.insert({FT2Scans[i].scanNumber, FT2Scans.data() + i});
            }
        }
#pragma omp section
        {
            ftFileReader FT1FileReader(FTfileBasePath + ".FT1");
            FT1FileReader.readAllScan();
            FT1Scans = std::move(FT1FileReader.Scans);
            FT1FileReader.Scans = {};
            scanNumerFT1ScanMap.reserve(FT1Scans.size());
            for (size_t i = 0; i < FT1Scans.size(); i++)
            {
                scanNumerFT1ScanMap.insert({FT1Scans[i].scanNumber, FT1Scans.data() + i});
            }
        }
    }
}

void ompForDemo()
{
    int numThreads = omp_get_max_threads();
    int numIterations = 10000;
    double sum = 0.0;

    // #pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < numIterations; i++)
    {
        double result = std::sqrt(std::pow(i, 20) + std::pow(i + 1, 20));
        sum += result;
    }

    printf("Number of threads: %d\n", numThreads);
    printf("Number of iterations: %d\n", numIterations);
    printf("Sum: %f\n", sum);
}

void PSMfeatureExtractor::extractFeaturesOfEachPSM()
{
    for (size_t i = 0; i < mSipPSM->isotopicPeakss.size(); i++)
    {
        mSipPSM->isotopicPeakss[i] = findIsotopicPeaks(mSipPSM->precursorScanNumbers[i],
                                                       mSipPSM->parentCharges[i],
                                                       mSipPSM->measuredParentMasses[i],
                                                       mSipPSM->calculatedParentMasses[i]);
        mSipPSM->istopicPeakNumbers[i] = mSipPSM->isotopicPeakss[i].size();
        mSipPSM->MS1IsotopicAbundances[i] = getSIPelementAbundanceFromMS1(mSipPSM->identifiedPeptides[i],
                                                                          mSipPSM->isotopicPeakss[i],
                                                                          mSipPSM->parentCharges[i]);
    }
}

void PSMfeatureExtractor::extractPSMfeature(const std::string &Spe2PepFilePath, const int topN,
                                            const std::string &ftFilepath)
{
    mSpe2PepFileReader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(Spe2PepFilePath, topN);
    for (size_t i = 0; i < mSpe2PepFileReader.FT2s.size(); i++)
    {
        // mExtractor.loadFT1FT2fileParallel(ftFilepath + "/" + mSpe2PepFileReader.FT2s[i]);
        loadFT1(ftFilepath + "/" + mSpe2PepFileReader.FT2s[i]);
        mSipPSM = &mSpe2PepFileReader.sipPSMs[i];
        // init vectors for features to be extracted;
        mSipPSM->isotopicPeakss = std::vector<std::vector<isotopicPeak>>(mSipPSM->scanNumbers.size());
        mSipPSM->istopicPeakNumbers = std::vector<int>(mSipPSM->scanNumbers.size());
        extractFeaturesOfEachPSM();
    }
}

void PSMfeatureExtractor::extractPSMfeatureParallel(const std::string &Spe2PepFilePath, const int topN,
                                                    const std::string &ftFilepath, const int threadNumber)
{
    mSpe2PepFileReader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(Spe2PepFilePath, topN);
    int num_cores = omp_get_num_procs();
    // Set the number of threads to the lesser of num_cores and 10
    int num_threads = std::min(num_cores, 10);
    num_threads = std::min(num_threads, threadNumber);
    omp_set_num_threads(num_threads);
    // Enable nested parallelism
    omp_set_nested(1);
#pragma omp parallel for
    for (size_t i = 0; i < mSpe2PepFileReader.FT2s.size(); i++)
    {
        PSMfeatureExtractor mExtractor;
        // mExtractor.loadFT1FT2fileParallel(ftFilepath + "/" + mSpe2PepFileReader.FT2s[i]);
        mExtractor.loadFT1(ftFilepath + "/" + mSpe2PepFileReader.FT2s[i]);
        mExtractor.mSipPSM = &mSpe2PepFileReader.sipPSMs[i];
        // init vectors for features to be extracted;
        mExtractor.mSipPSM->isotopicPeakss = std::vector<std::vector<isotopicPeak>>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->istopicPeakNumbers = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->MS1IsotopicAbundances = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.extractFeaturesOfEachPSM();
    }
}

void PSMfeatureExtractor::writeTSV(const std::string &fileName)
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    std::ofstream file(fileName);
    if (!file)
    {
        std::cerr << "Unable to open file for writing.\n";
    }
    file << std::fixed << std::setprecision(6);

    const size_t chunkSize = 10000;
    std::stringstream ss;
    ss << "fileNames"
       << "\t"
       << "scanNumbers"
       << "\t"
       << "precursorScanNumbers"
       << "\t"
       << "parentCharges"
       << "\t"
       << "isolationWindowCenterMZs"
       << "\t"
       << "measuredParentMasses"
       << "\t"
       << "calculatedParentMasses"
       << "\t"
       << "searchNames"
       << "\t"
       << "retentionTimes"
       << "\t"
       << "MVHscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "WDPscores"
       << "\t"
       << "ranks"
       << "\t"
       << "identifiedPeptides"
       << "\t"
       << "originalPeptides"
       << "\t"
       << "proteinNames"
       << "\n";
    std::vector<sipPSM> &sipPSMs = mSpe2PepFileReader.sipPSMs;
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        for (size_t j = 0; j < sipPSMs[i].fileNames.size(); j += chunkSize)
        {
            for (size_t k = j; k < std::min(k + chunkSize, sipPSMs[i].fileNames.size()); ++k)
            {
                ss << sipPSMs[i].fileNames[k] << "\t" << sipPSMs[i].scanNumbers[k] << "\t"
                   << sipPSMs[i].precursorScanNumbers[k] << "\t"
                   << sipPSMs[i].parentCharges[k] << "\t"
                   << sipPSMs[i].isolationWindowCenterMZs[k] << "\t"
                   << sipPSMs[i].measuredParentMasses[k] << "\t" << sipPSMs[i].calculatedParentMasses[k] << "\t"
                   << sipPSMs[i].searchNames[k] << "\t" << sipPSMs[i].retentionTimes[k] << "\t" << sipPSMs[i].MVHscores[k] << "\t"
                   << sipPSMs[i].XcorrScores[k] << "\t" << sipPSMs[i].WDPscores[k] << "\t"
                   << sipPSMs[i].ranks[k] << "\t" << sipPSMs[i].identifiedPeptides[k] << "\t"
                   << sipPSMs[i].originalPeptides[k] << "\t" << sipPSMs[i].proteinNames[k] << "\n";
            }
            file << ss.str();
            ss.str(std::string()); // Clear the stringstream
        }
    }
    file.close();
}