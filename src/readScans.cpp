#include "lib/ftFileReader.h"
#include <Rcpp.h>

using namespace Rcpp;

//' read FT file header
//' @param ftFile a ft1 file's full path
//' @return a list of ft file header
//' @examples
//' header <- readFTheader("demo.ft1")
//' @export
// [[Rcpp::export]]
List readFTheader(CharacterVector ftFile)
{
    ftFileReader reader(as<string>(ftFile));
    reader.detectPrecursorAndCharge();
    List header = List::create(Named("instrument") = reader.instrument,
                               _["scanType"] = reader.scanType,
                               _["scanFilter"] = reader.scanFilter,
                               _["hasPrecursor"] = reader.hasPrecursor,
                               _["hasCharge"] = reader.hasCharge);
    return header;
}

//' read MS1 scans with scanNumber as index
//' @param ftFile a ft1 file's full path
//' @return a list of MS1 scans with names of scan number
//' @examples
//' ft1 <- readAllScanMS1("demo.ft1")
//' @export
// [[Rcpp::export]]
List readScansMS1(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<string>(ftFile));
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    size_t scanCountInt = as<int>(scanCount);
    reader.readScans(scanCountInt);
    // in case there is not enough scan in ft file
    scanCountInt = scanCountInt < reader.Scans.size() ? scanCountInt : reader.Scans.size();
    List scanList(scanCountInt);
    Scan *mScan;
    for (size_t i = 0; i < scanCountInt; i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = move(mScan->mz),
                                             _["intensity"] = move(mScan->intensity),
                                             _["resolution"] = move(mScan->resolution),
                                             _["baseLine"] = move(mScan->baseLine),
                                             _["signalToNoise"] = move(mScan->signalToNoise),
                                             _["charge"] = move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["TIC"] = mScan->TIC,
                                      _["peaks"] = move(peakDf));
        scanList[i] = move(mScanList);
    }
    return scanList;
}

//' readAllScanMS1
//' @param ftFile a ft1 file's full path
//' @export
// [[Rcpp::export]]
List readAllScanMS1(CharacterVector ftFile)
{
    ftFileReader reader(as<string>(ftFile));
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    reader.readAllScan();
    List scanList(reader.Scans.size());
    CharacterVector scanNumbers(reader.Scans.size());
    Scan *mScan;
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        // use move to speed up copy vector
        // use wrap to convert vector to numericVector
        DataFrame peakDf = DataFrame::create(Named("mz") = wrap(move(mScan->mz)),
                                             _["intensity"] = wrap(move(mScan->intensity)),
                                             _["resolution"] = wrap(move(mScan->resolution)),
                                             _["baseLine"] = wrap(move(mScan->baseLine)),
                                             _["signalToNoise"] = wrap(move(mScan->signalToNoise)),
                                             _["charge"] = wrap(move(mScan->charge)));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["TIC"] = mScan->TIC,
                                      _["peaks"] = move(peakDf));
        scanList[i] = move(mScanList);
        scanNumbers[i] = to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}

//' readScansMS2
//' @param ftFile a ft2 file's full path
//' @param scanNumber the scanNumber th scan
//' @export
// [[Rcpp::export]]
List readScansMS2(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<string>(ftFile));
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    size_t scanCountInt = as<int>(scanCount);
    reader.readScans(scanCountInt);
    // in case there is not enough scan in ft file
    scanCountInt = scanCountInt < reader.Scans.size() ? scanCountInt : reader.Scans.size();
    List scanList(scanCountInt);
    Scan *mScan;
    for (size_t i = 0; i < scanCountInt; i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = move(mScan->mz),
                                             _["intensity"] = move(mScan->intensity),
                                             _["resolution"] = move(mScan->resolution),
                                             _["baseLine"] = move(mScan->baseLine),
                                             _["signalToNoise"] = move(mScan->signalToNoise),
                                             _["charge"] = move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorMz"] = mScan->precursorMz,
                                      _["TIC"] = mScan->TIC,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["peaks"] = move(peakDf));
        scanList[i] = move(mScanList);
    }
    return scanList;
}

//' read MS2 scans with scanNumber as index
//' @param ftFile a ft2 file's full path
//' @return a list of MS2 scans with names of scan number
//' @examples
//' ft2 <- readAllScanMS2("demo.ft2")
//' @export
// [[Rcpp::export]]
List readAllScanMS2(CharacterVector ftFile)
{
    ftFileReader reader(as<string>(ftFile));
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    reader.readAllScan();
    List scanList(reader.Scans.size());
    Scan *mScan;
    CharacterVector scanNumbers(reader.Scans.size());
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        // use move to speed up copy vector
        DataFrame peakDf = DataFrame::create(Named("mz") = wrap(move(mScan->mz)),
                                             _["intensity"] = wrap(move(mScan->intensity)),
                                             _["resolution"] = wrap(move(mScan->resolution)),
                                             _["baseLine"] = wrap(move(mScan->baseLine)),
                                             _["signalToNoise"] = wrap(move(mScan->signalToNoise)),
                                             _["charge"] = wrap(move(mScan->charge)));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorMz"] = mScan->precursorMz,
                                      _["TIC"] = mScan->TIC,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["peaks"] = move(peakDf));
        scanList[i] = move(mScanList);
        scanNumbers[i] = to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}
