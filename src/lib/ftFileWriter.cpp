#include "ftFileWriter.h"

ftFileWriter::ftFileWriter() {}

ftFileWriter::ftFileWriter(std::string file, std::string mInstrument,
                           std::string mScanType, std::string mScanFilter,
                           bool mHasCharge, int mMSlevel) : ftFileName(file),
                                                            instrument(mInstrument),
                                                            scanType(mScanType),
                                                            scanFilter(mScanFilter),
                                                            hasCharge(mHasCharge),
                                                            MSlevel(mMSlevel)
{
    // speedup writing
    // setlocale(LC_ALL, "C");
    // ios_base::sync_with_stdio(false);
    ftFileStream = std::stringstream("", std::ios_base::app | std::ios_base::out);
    // in case overwrite existing file
    if (!fs::exists(ftFileName))
    {
        ftFileStreamtoFile.open(ftFileName.c_str(), std::ios::out);
        if (ftFileStreamtoFile.is_open())
        {
            isWriteable = true;
        }
        else
            std::cout << "Cannot write " << ftFileName << std::endl;
    }
    else
        std::cout << "Old " << ftFileName << " exists" << std::endl;
}

ftFileWriter::~ftFileWriter()
{
    if (ftFileStreamtoFile.is_open())
        ftFileStreamtoFile.close();
}

std::vector<int> ftFileWriter::rankify(const std::vector<double> &A)
{
    int n = A.size();
    std::vector<int> R(n);
    // T[][0] is the data and T[][1] is
    // the index of data in A
    std::vector<std::pair<double, int>> T(n);
    for (int i = 0; i < n; i++)
    {
        T[i].first = A[i];
        T[i].second = i;
    }
    // Sort T according to first element
    std::sort(T.begin(), T.end(), [](const std::pair<double, int> &a, const std::pair<double, int> &b) -> bool
              { return (a.first > b.first); });
    for (int i = 0; i < n; i++)
    {
        R[T[i].second] = i;
    }
    return (R);
}

void ftFileWriter::getMeanRanksOfPeaks(std::vector<float> &ranks,
                                       Scan const &mScan, const double &window,
                                       const double &step)
{
    int n = mScan.mz.size();
    // in case devide by 0
    std::vector<int> windowsCount(n, 1);
    double start = mScan.mz.front();
    double end = start + window;
    int startIX = 0, IX = 0;
    bool moveWindow = true;
    // get sum and number of ranks of each moving window
    while (end < mScan.mz.back() + window)
    {
        if (moveWindow)
        {
            while (mScan.mz[IX] < end && IX < n)
            {
                IX++;
                // moveWindow = true;
            }
            std::vector<int> tempRanks(rankify(
                std::vector<double>(mScan.intensity.begin() + startIX, mScan.intensity.begin() + IX)));
            for (size_t i = 0; i < tempRanks.size(); i++)
            {
                ranks[startIX + i] += (float)tempRanks[i];
                windowsCount[startIX + i]++;
            }
            start += step;
            end = start + window;
            moveWindow = false;
        }
        startIX++;
        IX = startIX;
        if (mScan.mz[IX] >= start)
            moveWindow = true;
    }
    // average the ranks
    for (size_t i = 0; i < ranks.size(); i++)
    {
        ranks[i] = ranks[i] / (float)windowsCount[i];
    }
}

void ftFileWriter::denoiseMS2ScanHasCharge(Scan &mScan, const double &window, const double &step,
                                           const float &threshold)
{
    int n = mScan.mz.size();
    std::vector<float> ranks(n);
    getMeanRanksOfPeaks(ranks, mScan, window, step);
    std::vector<double> mz, intensity;
    mz.reserve(n);
    intensity.reserve(n);
    std::vector<int> charge, resolution;
    charge.reserve(n);
    resolution.reserve(n);
    std::vector<float> baseLine, signalToNoise;
    baseLine.reserve(n);
    signalToNoise.reserve(n);
    // filter out the top threshold peaks
    for (size_t i = 0; i < ranks.size(); i++)
    {
        if (ranks[i] < threshold)
        {
            mz.push_back(mScan.mz[i]);
            intensity.push_back(mScan.intensity[i]);
            charge.push_back(mScan.charge[i]);
            resolution.push_back(mScan.resolution[i]);
            baseLine.push_back(mScan.baseLine[i]);
            signalToNoise.push_back(mScan.signalToNoise[i]);
            // for debug
            // signalToNoise.push_back(ranks[i]);
        }
    }
    mScan.mz = mz;
    mScan.intensity = intensity;
    mScan.charge = charge;
    mScan.resolution = resolution;
    mScan.baseLine = baseLine;
    mScan.signalToNoise = signalToNoise;
}

void ftFileWriter::writeHeaderHasCharge()
{
    ftFileStream << "H\tExtractor\tAerith V0.1" << std::endl;
    ftFileStream << "H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge" << std::endl;
    ftFileStream << "H\tInstrument Model\t" << instrument << std::endl;
}

void ftFileWriter::writeMS2ScanHasCharge(Scan &mScan)
{
    ftFileStream << std::fixed << std::setprecision(6);
    ftFileStream << "S\t" << std::to_string(mScan.scanNumber) << "\t" << mScan.isolationWindowCenterMZ
                 << std::setprecision(2) << "\t" << mScan.TIC << std::endl;
    ftFileStream << "Z\t" << std::to_string(mScan.precursorCharge)
                 << std::setprecision(6) << "\t" << mScan.isolationWindowCenterMZ * mScan.precursorCharge;
    for (size_t i = 0; i < mScan.precursorCharges.size(); i++)
    {
        ftFileStream << "\t" << std::to_string(mScan.precursorCharges[i]);
        ftFileStream << "\t" << std::setprecision(6) << mScan.precursorMZs[i];
    }
    ftFileStream << std::endl;
    ftFileStream << "I\tRetentionTime\t" << mScan.retentionTime << std::endl;
    ftFileStream << "I\tScanType\t" << scanType << std::endl;
    ftFileStream << "I\tScanFilter\t" << scanFilter << std::endl;
    ftFileStream << "D\tParentScanNumber\t" << std::to_string(mScan.precursorScanNumber) << std::endl;
    for (size_t i = 0; i < mScan.mz.size(); i++)
    {
        ftFileStream << std::setprecision(6) << mScan.mz[i]
                     << "\t" << std::setprecision(2) << mScan.intensity[i]
                     << "\t" << std::to_string(mScan.resolution[i])
                     << "\t" << mScan.baseLine[i]
                     << "\t" << mScan.signalToNoise[i] << "\t" << std::to_string(mScan.charge[i])
                     << std::endl;
    }
}

void ftFileWriter::writeMS1ScanHasCharge(Scan &mScan)
{
    ftFileStream << std::fixed << std::setprecision(6);
    ftFileStream << "S\t" << std::to_string(mScan.scanNumber)
                 << std::setprecision(2) << "\t" << mScan.TIC << std::endl;
    ftFileStream << "I\tRetentionTime\t" << mScan.retentionTime << std::endl;
    ftFileStream << "I\tScanType\t" << scanType << std::endl;
    ftFileStream << "I\tScanFilter\t" << scanFilter << std::endl;
    for (size_t i = 0; i < mScan.mz.size(); i++)
    {
        ftFileStream << std::setprecision(6) << mScan.mz[i]
                     << "\t" << std::setprecision(2) << mScan.intensity[i]
                     << "\t" << std::to_string(mScan.resolution[i])
                     << "\t" << mScan.baseLine[i]
                     << "\t" << mScan.signalToNoise[i] << "\t" << std::to_string(mScan.charge[i])
                     << std::endl;
    }
}

void ftFileWriter::writeAllScanMS1(std::vector<Scan> &mScans)
{
    if (!isWriteable)
        return;
    if (hasCharge && MSlevel == 1)
    {
        writeHeaderHasCharge();
        for (size_t i = 0; i < mScans.size(); i++)
        {
            writeMS1ScanHasCharge(mScans[i]);
            // output buffer to FT file
            ftFileStreamtoFile << ftFileStream.str();
            ftFileStream.str("");
        }
    }
    else
        std::cout << "This format is not implemented!" << std::endl;
}

void ftFileWriter::writeAllScanMS2(std::vector<Scan> &mScans)
{
    if (!isWriteable)
        return;
    if (hasCharge && MSlevel == 2)
    {
        writeHeaderHasCharge();
        for (size_t i = 0; i < mScans.size(); i++)
        {
            writeMS2ScanHasCharge(mScans[i]);
            // output buffer to FT file
            ftFileStreamtoFile << ftFileStream.str();
            ftFileStream.str("");
        }
    }
    else
        std::cout << "This format is not implemented!" << std::endl;
}