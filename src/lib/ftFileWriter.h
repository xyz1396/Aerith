#pragma once
#include "ftFileReader.h"
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
    void writeHeaderHasCharge();
    void writeMS2ScanHasCharge(Scan &mScan);
    void writeAllScanMS2(vector<Scan> &mScans);
};
