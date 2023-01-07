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