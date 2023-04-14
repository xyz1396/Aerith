#include "ftFileWriter.h"

ftFileWriter::ftFileWriter() {}

ftFileWriter::ftFileWriter(string file, string mInstrument,
                           string mScanType, string mScanFilter,
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
    ftFileStream = stringstream("", ios_base::app | ios_base::out);
    // in case overwrite existing file
    if (!fs::exists(ftFileName))
    {
        ftFileStreamtoFile.open(ftFileName.c_str(), ios::out);
        if (ftFileStreamtoFile.is_open())
        {
            isWriteable = true;
        }
        else
            cout << "Cannot write " << ftFileName << endl;
    }
    else
        cout << "Old " << ftFileName << " exists" << endl;
}

ftFileWriter::~ftFileWriter()
{
    if (ftFileStreamtoFile.is_open())
        ftFileStreamtoFile.close();
}

vector<int> ftFileWriter::rankify(const vector<double> &A)
{
    int n = A.size();
    vector<int> R(n);
    // T[][0] is the data and T[][1] is
    // the index of data in A
    vector<pair<double, int>> T(n);
    for (int i = 0; i < n; i++)
    {
        T[i].first = A[i];
        T[i].second = i;
    }
    // Sort T according to first element
    sort(T.begin(), T.end(), [](const pair<double, int> &a, const pair<double, int> &b) -> bool
         { return (a.first > b.first); });
    for (int i = 0; i < n; i++)
    {
        R[T[i].second] = i;
    }
    return (R);
}

void ftFileWriter::getMeanRanksOfPeaks(vector<float> &ranks,
                                       Scan const &mScan, const double &window,
                                       const double &step)
{
    int n = mScan.mz.size();
    // in case devide by 0
    vector<int> windowsCount(n, 1);
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
            vector<int> tempRanks(rankify(
                vector<double>(mScan.intensity.begin() + startIX, mScan.intensity.begin() + IX)));
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
    vector<float> ranks(n);
    getMeanRanksOfPeaks(ranks, mScan, window, step);
    vector<double> mz, intensity;
    mz.reserve(n);
    intensity.reserve(n);
    vector<int> charge, resolution;
    charge.reserve(n);
    resolution.reserve(n);
    vector<float> baseLine, signalToNoise;
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
    ftFileStream << "H\tExtractor\tAerith V0.1" << endl;
    ftFileStream << "H\tm/z\tIntensity\tResolution\tBaseline\tNoise\tCharge" << endl;
    ftFileStream << "H\tInstrument Model\t" << instrument << endl;
}

void ftFileWriter::writeMS2ScanHasCharge(Scan &mScan)
{
    ftFileStream << fixed << setprecision(6);
    ftFileStream << "S\t" << to_string(mScan.scanNumber) << "\t" << mScan.precursorMz
                 << setprecision(2) << "\t" << mScan.TIC << endl;
    ftFileStream << "Z\t" << to_string(mScan.precursorCharge)
                 << setprecision(6) << "\t" << mScan.precursorMz * mScan.precursorCharge << endl;
    ftFileStream << "I\tRetentionTime\t" << mScan.retentionTime << endl;
    ftFileStream << "I\tScanType\t" << scanType << endl;
    ftFileStream << "I\tScanFilter\t" << scanFilter << endl;
    ftFileStream << "D\tParentScanNumber\t" << to_string(mScan.precursorScanNumber) << endl;
    for (size_t i = 0; i < mScan.mz.size(); i++)
    {
        ftFileStream << setprecision(6) << mScan.mz[i]
                     << "\t" << setprecision(2) << mScan.intensity[i]
                     << "\t" << to_string(mScan.resolution[i])
                     << "\t" << mScan.baseLine[i]
                     << "\t" << mScan.signalToNoise[i] << "\t" << to_string(mScan.charge[i])
                     << endl;
    }
}

void ftFileWriter::writeAllScanMS2(vector<Scan> &mScans)
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
        cout << "This format is not implemented!" << endl;
}