#pragma once
#include "ftFileReader.h"
#include <algorithm>
#include <sstream>

class ftFileWriter
{
private:
public:
    ftFileWriter();
    ftFileWriter(std::string file, std::string mInstrument, std::string mScanType,
                 std::string mScanFilter, bool mHasCharge, int mMSlevel);
    ~ftFileWriter();
    std::string ftFileName, instrument, scanType, scanFilter;
    // buffer for each scan;
    std::stringstream ftFileStream;
    std::ofstream ftFileStreamtoFile;
    // in case write error
    bool isWriteable = false, hasCharge = false;
    int MSlevel = 1;
    // return rank of all elements in vector from high to low
    std::vector<int> rankify(const std::vector<double> &A);
    // length of ranks should be as same as length of mScan.mz
    void getMeanRanksOfPeaks(std::vector<float> &ranks, Scan const &mScan,
                             const double &window, const double &step);
    void denoiseMS2ScanHasCharge(Scan &mScan, const double &window, const double &step,
                                 const float &threshold);
    void writeHeaderHasCharge();
    void writeMS2ScanHasCharge(Scan &mScan);
    void writeMS1ScanHasCharge(Scan &mScan);
    void writeAllScanMS2(std::vector<Scan> &mScans);
    void writeAllScanMS1(std::vector<Scan> &mScans);
};
