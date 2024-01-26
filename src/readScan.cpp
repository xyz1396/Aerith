#include "lib/ftFileReader.h"
#include <Rcpp.h>
using namespace Rcpp;

//' readOneScanMS2
//' @param ftFile a ft2 file's full path
//' @param scanNumber the scan at scanNumber
//' @return a list of MS2 scan
//' @examples
//' ft2 <- readOneScanMS2("demo.ft2", 2)
//' @export
// [[Rcpp::export]]
List readOneScanMS2(const String &ftFile, const size_t scanNumber)
{
    ftFileReader reader(ftFile);
    Scan mScan = reader.readOneScan(scanNumber);
    DataFrame peakDf = DataFrame::create(Named("mz") = mScan.mz,
                                         _["intensity"] = mScan.intensity,
                                         _["resolution"] = mScan.resolution,
                                         _["baseLine"] = mScan.baseLine,
                                         _["signalToNoise"] = mScan.signalToNoise,
                                         _["charge"] = mScan.charge);
    List mScanList = List::create(Named("scanNumber") = mScan.scanNumber,
                                  _["retentionTime"] = mScan.retentionTime,
                                  _["precursorScanNumber"] = mScan.precursorScanNumber,
                                  _["precursorMz"] = mScan.precursorMz,
                                  _["precursorCharge"] = mScan.precursorCharge,
                                  _["peaks"] = peakDf);
    return mScanList;
}

//' readOneScanMS1
//' @param ftFile a ft1 file's full path
//' @param scanNumber the scan at scanNumber
//' @return a list of MS1 scan
//' @examples
//' ft1 <- readOneScanMS1("demo.ft1", 2)
//' @export
// [[Rcpp::export]]
List readOneScanMS1(const String &ftFile, const size_t scanNumber)
{
    ftFileReader reader(ftFile);
    Scan mScan = reader.readOneScan(scanNumber);
    DataFrame peakDf = DataFrame::create(Named("mz") = mScan.mz,
                                         _["intensity"] = mScan.intensity,
                                         _["resolution"] = mScan.resolution,
                                         _["baseLine"] = mScan.baseLine,
                                         _["signalToNoise"] = mScan.signalToNoise,
                                         _["charge"] = mScan.charge);
    List mScanList = List::create(Named("scanNumber") = mScan.scanNumber,
                                  _["retentionTime"] = mScan.retentionTime,
                                  _["peaks"] = peakDf);
    return mScanList;
}