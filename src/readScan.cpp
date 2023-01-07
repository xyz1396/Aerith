#include "lib/ftFileReader.h"
#include <Rcpp.h>
using namespace Rcpp;

//' readOneScanMS2
//' @param ftFile a ft2 file's full path
//' @param scanCount the scanCount'th scan
//' @export
// [[Rcpp::export]]
List readOneScanMS2(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<std::string>(ftFile));
    Scan mScan = reader.readOneScan(as<int>(scanCount));
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
//' @param scanCount the scanCount'th scan
//' @export
// [[Rcpp::export]]
List readOneScanMS1(CharacterVector ftFile, NumericVector scanCount)
{
    ftFileReader reader(as<std::string>(ftFile));
    Scan mScan = reader.readOneScan(as<int>(scanCount));
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