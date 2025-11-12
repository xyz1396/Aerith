#include <Rcpp.h>

#include "lib/ftFileReader.h"
using namespace Rcpp;

//' readOneScanMS2
//' @param ftFile a ft2 file's full path
//' @param scanNumber the scan at scanNumber
//' @return a list of MS2 scan
//' @examples
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' ft2 <- readOneScanMS2(demo_file, 1633)
//' @export
// [[Rcpp::export]]
List readOneScanMS2(const String &ftFile, const size_t scanNumber) {
    ftFileReader reader(ftFile);
    Scan mScan = reader.readOneScan(scanNumber);
    DataFrame peakDf = DataFrame::create(
        Named("mz") = mScan.mz, _["intensity"] = mScan.intensity,
        _["resolution"] = mScan.resolution, _["baseLine"] = mScan.baseLine,
        _["signalToNoise"] = mScan.signalToNoise, _["charge"] = mScan.charge);
    List mScanList = List::create(
        Named("scanNumber") = mScan.scanNumber,
        _["retentionTime"] = mScan.retentionTime,
        _["precursorScanNumber"] = mScan.precursorScanNumber,
        _["precursorCharge"] = mScan.precursorCharge,
        _["isolationWindowCenterMZ"] = mScan.isolationWindowCenterMZ,
        _["TIC"] = mScan.TIC, _["precursorCharges"] = mScan.precursorCharges,
        _["precursorMZs"] = mScan.precursorMZs, _["peaks"] = std::move(peakDf));
    return mScanList;
}

//' readOneScanMS1
//' @param ftFile a ft1 file's full path
//' @param scanNumber the scan at scanNumber
//' @return a list of MS1 scan
//' @examples
//' rds <- system.file("extdata", "demo.FT1.rds", package = "Aerith")
//' demo_file <- tempfile(fileext = ".FT1")
//' writeLines(readRDS(rds), demo_file)
//' ft1 <- readOneScanMS1(demo_file, 1588)
//' @export
// [[Rcpp::export]]
List readOneScanMS1(const String &ftFile, const size_t scanNumber) {
    ftFileReader reader(ftFile);
    Scan mScan = reader.readOneScan(scanNumber);
    DataFrame peakDf = DataFrame::create(
        Named("mz") = mScan.mz, _["intensity"] = mScan.intensity,
        _["resolution"] = mScan.resolution, _["baseLine"] = mScan.baseLine,
        _["signalToNoise"] = mScan.signalToNoise, _["charge"] = mScan.charge);
    List mScanList = List::create(Named("scanNumber") = mScan.scanNumber,
                                  _["retentionTime"] = mScan.retentionTime,
                                  _["TIC"] = mScan.TIC, _["peaks"] = peakDf);
    return mScanList;
}