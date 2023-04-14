#pragma once
#include "ftFileReader.h"
#include <algorithm>
#include <sstream>

class ftFileWriter
{
private:
public:
    ftFileWriter();
    ftFileWriter(string file, string mInstrument, string mScanType,
                 string mScanFilter, bool mHasCharge, int mMSlevel);
    ~ftFileWriter();
    string ftFileName, instrument, scanType, scanFilter;
    // buffer for each scan;
    stringstream ftFileStream;
    ofstream ftFileStreamtoFile;
    // in case write error
    bool isWriteable = false, hasCharge = false;
    int MSlevel = 1;
    // return rank of all elements in vector from high to low
    vector<int> rankify(const vector<double> &A);
    // length of ranks should be as same as length of mScan.mz
    void getMeanRanksOfPeaks(vector<float> &ranks, Scan const &mScan,
                                      const double &window, const double &step);
    void denoiseMS2ScanHasCharge(Scan &mScan, const double &window, const double &step,
                                 const float &threshold);
    void writeHeaderHasCharge();
    void writeMS2ScanHasCharge(Scan &mScan);
    void writeAllScanMS2(vector<Scan> &mScans);
};
