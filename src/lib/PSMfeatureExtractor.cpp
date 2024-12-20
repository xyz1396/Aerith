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
            // if (mScan->charge[mid] == charge && mScan->intensity[mid] > currentIntensity)
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
    // go left slope of isotopic envelope
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
    // go right slope of isotopic envelope
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
    findIsotopicPeaks(int &MS1ScanNumber,
                      const int precursorCharge,
                      const double observedPrecursorMass,
                      const double calculatedPrecursorMass)
{
    std::vector<isotopicPeak> isotopicPeaks = {};
    double observedPrecursorMZ = observedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    double calculatedPrecursorMZ = calculatedPrecursorMass / precursorCharge + ProNovoConfig::getProtonMass();
    Scan *MS1Scan = scanNumerFT1ScanMap[MS1ScanNumber];
    size_t peakIX = std::numeric_limits<size_t>::max();
    // if first search failed, search 2 scans before it
    for (int i = 0; i < 3; i++)
    {
        peakIX = binarySearchPeak(MS1Scan, observedPrecursorMZ, precursorCharge);
        if (peakIX != std::numeric_limits<size_t>::max())
        {
            MS1ScanNumber = MS1Scan->scanNumber;
            break;
        }
        else if (MS1Scan != FT1Scans.data())
            MS1Scan--;
        else
            break;
    }
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

double PSMfeatureExtractor::getSIPelementAbundanceFromMS1(const std::string &nakePeptide,
                                                          const std::vector<isotopicPeak> &isotopicPeaks,
                                                          const int precursorCharge)
{
    // in case of no isotopic peaks
    if (isotopicPeaks.size() == 0)
        return 0.0;
    double baseMass = mAveragine.calPrecursorBaseMass(nakePeptide);
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
    double atomCnumber = mAveragine.pepAtomCounts[mAveragine.SIPatomIX];
    pct /= atomCnumber;
    pct *= 100.;
    return pct;
}

std::pair<int, int> PSMfeatureExtractor::getSeqLengthAndMissCleavageSiteNumber(const std::string &peptideSeq)
{
    std::size_t start = peptideSeq.find_first_of('[');
    std::size_t end = peptideSeq.find_last_of(']');
    int seqLength = end - start - 1;
    // find miss cleavage site not at margin
    std::string seq = peptideSeq.substr(start + 2, end - start - 3);
    int count = 0;
    for (char &A : cleavageSites)
    {
        for (char &S : seq)
        {
            if (S == A)
                count++;
        }
    }
    return {seqLength, count};
}

int PSMfeatureExtractor::getPTMnumber(const std::string &peptideSeq)
{
    int count = 0;
    for (char c : peptideSeq)
    {
        if (!isalpha(c))
        {
            count++;
        }
    }
    return count;
}

std::pair<int, double> PSMfeatureExtractor::getMassWindowShiftAndError(const double observedPrecursorMass,
                                                                       const double calculatedPrecursorMass)
{
    int massWindowShift = static_cast<int>(round(std::abs(observedPrecursorMass - calculatedPrecursorMass) /
                                                 ProNovoConfig::getNeutronMass()));
    double massError = std::fmod(std::abs(observedPrecursorMass - calculatedPrecursorMass),
                                 ProNovoConfig::getNeutronMass());
    if (massError > ProNovoConfig::getNeutronMass() / 2)
    {
        massError = ProNovoConfig::getNeutronMass() - massError;
    }
    // convert it to ppm
    massError = massError / calculatedPrecursorMass * 1000000;
    return {massWindowShift, massError};
}

double PSMfeatureExtractor::getMS2IsotopicAbundance(const std::string &searchName)
{
    if (searchName == "SE")
        return 1.07;
    std::string pct = searchName;
    std::string delimiter = "_";
    size_t pos = 0;
    while ((pos = pct.find(delimiter)) != std::string::npos)
    {
        pct.erase(0, pos + delimiter.length());
    }
    pct = pct.substr(0, pct.size() - 3);
    double pct_num = std::stod(pct);
    pct_num = pct_num;
    return pct_num;
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

void PSMfeatureExtractor::loadFT2(const std::string &FTfileBasePath)
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

void PSMfeatureExtractor::loadFT1FT2fileParallel(const std::string &FTfileBasePath)
{
#pragma omp parallel sections
    {
#pragma omp section
        loadFT1(FTfileBasePath);
#pragma omp section
        loadFT2(FTfileBasePath);
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
    float topScore = 0;
    for (size_t i = 0; i < mSipPSM->isotopicPeakss.size(); i++)
    {
        mSipPSM->isotopicPeakss[i] = findIsotopicPeaks(mSipPSM->precursorScanNumbers[i],
                                                       mSipPSM->parentCharges[i],
                                                       mSipPSM->measuredParentMasses[i],
                                                       mSipPSM->calculatedParentMasses[i]);
        mSipPSM->isotopicPeakNumbers[i] = mSipPSM->isotopicPeakss[i].size();
        mSipPSM->MS1IsotopicAbundances[i] = getSIPelementAbundanceFromMS1(mSipPSM->nakePeptides[i],
                                                                          mSipPSM->isotopicPeakss[i],
                                                                          mSipPSM->parentCharges[i]);
        std::tie(mSipPSM->peptideLengths[i], mSipPSM->missCleavageSiteNumbers[i]) =
            getSeqLengthAndMissCleavageSiteNumber(mSipPSM->originalPeptides[i]);
        mSipPSM->PTMnumbers[i] = getPTMnumber(mSipPSM->nakePeptides[i]);

        if (mSipPSM->ranks[i] == 1)
        {
            // topScore = mSipPSM->MVHscores[i];
            topScore = mSipPSM->WDPscores[i];
        }
        // mSipPSM->MVHdiffScores[i] = topMVHscore - mSipPSM->MVHscores[i];
        mSipPSM->diffScores[i] = topScore - mSipPSM->WDPscores[i];

        mSipPSM->mzShiftFromisolationWindowCenters[i] = std::abs(
            mSipPSM->isolationWindowCenterMZs[i] -
            mSipPSM->measuredParentMasses[i] / mSipPSM->parentCharges[i] - ProNovoConfig::getProtonMass());
        std::tie(mSipPSM->isotopicMassWindowShifts[i], mSipPSM->massErrors[i]) = getMassWindowShiftAndError(
            mSipPSM->measuredParentMasses[i], mSipPSM->calculatedParentMasses[i]);

        mSipPSM->precursorIntensities[i] = 0;
        for (auto &peak : mSipPSM->isotopicPeakss[i])
        {
            mSipPSM->precursorIntensities[i] += peak.intensity;
        }

        mSipPSM->MS2IsotopicAbundances[i] = getMS2IsotopicAbundance(mSipPSM->searchNames[i]);
        mSipPSM->isotopicAbundanceDiffs[i] = mSipPSM->MS1IsotopicAbundances[i] - mSipPSM->MS2IsotopicAbundances[i];
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
        mSipPSM->isotopicPeakNumbers = std::vector<int>(mSipPSM->scanNumbers.size());
        mSipPSM->MS1IsotopicAbundances = std::vector<double>(mSipPSM->scanNumbers.size());
        mSipPSM->MS2IsotopicAbundances = std::vector<double>(mSipPSM->scanNumbers.size());
        mSipPSM->isotopicAbundanceDiffs = std::vector<double>(mSipPSM->scanNumbers.size());
        mSipPSM->peptideLengths = std::vector<int>(mSipPSM->scanNumbers.size());
        mSipPSM->missCleavageSiteNumbers = std::vector<int>(mSipPSM->scanNumbers.size());
        mSipPSM->diffScores = std::vector<float>(mSipPSM->scanNumbers.size());
        mSipPSM->mzShiftFromisolationWindowCenters = std::vector<double>(mSipPSM->scanNumbers.size());
        mSipPSM->isotopicMassWindowShifts = std::vector<int>(mSipPSM->scanNumbers.size());
        mSipPSM->massErrors = std::vector<double>(mSipPSM->scanNumbers.size());
        mSipPSM->precursorIntensities = std::vector<double>(mSipPSM->scanNumbers.size(), 0);
        extractFeaturesOfEachPSM();
    }
}

void PSMfeatureExtractor::extractPSMfeatureParallel(
    const std::string &ftFilepath, const int threadNumber)
{
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
        mExtractor.mSipPSM->isotopicPeakNumbers = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->MS1IsotopicAbundances = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->MS2IsotopicAbundances = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->isotopicAbundanceDiffs = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->peptideLengths = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->missCleavageSiteNumbers = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->PTMnumbers = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->diffScores = std::vector<float>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->mzShiftFromisolationWindowCenters = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->isotopicMassWindowShifts = std::vector<int>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->massErrors = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size());
        mExtractor.mSipPSM->precursorIntensities = std::vector<double>(mExtractor.mSipPSM->scanNumbers.size(), 0);
        mExtractor.extractFeaturesOfEachPSM();
    }
}

void PSMfeatureExtractor::extractPSMfeatureParallel(const std::string &Spe2PepFilePath, const int topN,
                                                    const std::string &ftFilepath, const int threadNumber)
{
    mSpe2PepFileReader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(Spe2PepFilePath, topN);
    extractPSMfeatureParallel(ftFilepath, threadNumber);
}

void PSMfeatureExtractor::extractPSMfeatureParallel(const std::string &targetPath, const std::string &decoyPath, const int topN,
                                                    const std::string &ftFilepath, const int threadNumber)
{
    mSpe2PepFileReader.readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(
        targetPath, decoyPath, topN);
    extractPSMfeatureParallel(ftFilepath, threadNumber);
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
            for (size_t k = j; k < std::min(j + chunkSize, sipPSMs[i].fileNames.size()); ++k)
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

void PSMfeatureExtractor::writePecorlatorPin(const std::string &fileName, bool doProteinInference)
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
    ss << "SpecId"
       << "\t"
       << "Label"
       << "\t"
       << "ScanNr"
       << "\t"
       << "ExpMass"
       << "\t"
       << "retentiontime"
       << "\t"
       << "ranks"
       << "\t"
       << "parentCharges"
       << "\t"
       << "massErrors"
       << "\t"
       << "isotopicMassWindowShifts"
       << "\t"
       << "mzShiftFromisolationWindowCenters"
       << "\t"
       << "peptideLengths"
       << "\t"
       << "missCleavageSiteNumbers"
       << "\t"
       << "PTMnumbers"
       << "\t"
       << "istopicPeakNumbers"
       << "\t"
       << "MS1IsotopicAbundances"
       << "\t"
       << "MS2IsotopicAbundances"
       << "\t"
       << "isotopicAbundanceDiffs"
       << "\t"
       << "WDPscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "MVHscores"
       << "\t"
       << "diffScores"
       << "\t"
       << "log10_precursorIntensities"
       << "\t"
       << "Peptide"
       << "\t"
       << "Proteins"
       << "\n";
    std::string proteinName;
    std::string peptideSeq;
    std::vector<sipPSM> &sipPSMs = mSpe2PepFileReader.sipPSMs;
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        for (size_t j = 0; j < sipPSMs[i].fileNames.size(); j += chunkSize)
        {
            for (size_t k = j; k < std::min(j + chunkSize, sipPSMs[i].fileNames.size()); ++k)
            {
                ss << sipPSMs[i].fileNames[k] << "." << sipPSMs[i].scanNumbers[k] << "." << sipPSMs[i].ranks[k] << "\t";
                ss << (sipPSMs[i].isDecoys[k] ? -1 : 1) << "\t";
                ss << sipPSMs[i].scanNumbers[k] << "\t"
                   << sipPSMs[i].calculatedParentMasses[k] << "\t"
                   << sipPSMs[i].retentionTimes[k] << "\t"
                   << sipPSMs[i].ranks[k] << "\t"
                   << sipPSMs[i].parentCharges[k] << "\t"
                   << sipPSMs[i].massErrors[k] << "\t"
                   << sipPSMs[i].isotopicMassWindowShifts[k] << "\t"
                   << sipPSMs[i].mzShiftFromisolationWindowCenters[k] << "\t"
                   << sipPSMs[i].peptideLengths[k] << "\t"
                   << sipPSMs[i].missCleavageSiteNumbers[k] << "\t"
                   << sipPSMs[i].PTMnumbers[k] << "\t"
                   << sipPSMs[i].isotopicPeakNumbers[k] << "\t"
                   << sipPSMs[i].MS1IsotopicAbundances[k] << "\t"
                   << sipPSMs[i].MS2IsotopicAbundances[k] << "\t"
                   << sipPSMs[i].isotopicAbundanceDiffs[k] << "\t"
                   << sipPSMs[i].WDPscores[k] << "\t"
                   << sipPSMs[i].XcorrScores[k] << "\t"
                   << sipPSMs[i].MVHscores[k] << "\t"
                   << sipPSMs[i].diffScores[k] << "\t"
                   << (sipPSMs[i].precursorIntensities[k] > 0
                           ? std::log10(sipPSMs[i].precursorIntensities[k])
                           : 0)
                   << "\t";

                if (doProteinInference)
                {
                    // for percolator protein inference format
                    peptideSeq = sipPSMs[i].originalPeptides[k];
                    if (peptideSeq.front() == '[')
                        peptideSeq.insert(peptideSeq.begin(), 'n');
                    if (peptideSeq.back() == ']')
                        peptideSeq.push_back('n');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), '[', '.');
                    std::replace(peptideSeq.begin(), peptideSeq.end(), ']', '.');
                }
                else
                    peptideSeq = sipPSMs[i].identifiedPeptides[k];
                ss << peptideSeq << "\t";

                proteinName = sipPSMs[i].proteinNames[k];
                if (doProteinInference)
                {
                    // for percolator protein inference format
                    proteinName = proteinName.substr(1, sipPSMs[i].proteinNames[k].length() - 2);
                    std::replace(proteinName.begin(), proteinName.end(), ',', '\t');
                }
                ss << proteinName;
                ss << "\n";
            }
            file << ss.str();
            ss.str(std::string()); // Clear the stringstream
        }
    }
    file.close();
}