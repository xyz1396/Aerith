#include "lib/ftFileWriter.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rankyfify(NumericVector a)
{
    ftFileWriter writer;
    vector<double> b = as<vector<double>>(a);
    return wrap(writer.rankify(b));
}

//' denoise one MS2 scan has charge
//' @param scanList a list of one MS2 scan has charge
//' @param window a float of mz window size for denoise
//' @param step a float of mz step for denoise
//' @param threshold a float of top N threshold for denoise
//' @return a denoised MS2 scan has charge
//' @examples
//' ft2 <- readAllScanMS2("demo.ft2")
//' ms2 <- denoiseOneMS2ScanHasCharge(ft2[1], 100, 10, 25)
//' @export
// [[Rcpp::export]]
List denoiseOneMS2ScanHasCharge(List scanList, float window, float step, float threshold)
{
    ftFileWriter writer;
    Scan scan;
    DataFrame peaks;
    List newScanList = clone(scanList);
    peaks = scanList["peaks"];
    scan.mz = as<vector<double>>(peaks["mz"]);
    scan.intensity = as<vector<double>>(peaks["intensity"]);
    scan.resolution = as<vector<int>>(peaks["resolution"]);
    scan.baseLine = as<vector<float>>(peaks["baseLine"]);
    scan.signalToNoise = as<vector<float>>(peaks["signalToNoise"]);
    scan.charge = as<vector<int>>(peaks["charge"]);
    writer.denoiseMS2ScanHasCharge(scan, window, step, threshold);
    DataFrame newPeaks = DataFrame::create(Named("mz") = scan.mz,
                                         _["intensity"] = scan.intensity,
                                         _["resolution"] = scan.resolution,
                                         _["baseLine"] = scan.baseLine,
                                         _["signalToNoise"] = scan.signalToNoise,
                                         _["charge"] = scan.charge);
    newScanList["peaks"] = newPeaks;
    return newScanList;
}

//' write all MS2 scans has charge
//' @param header a list of FT file header
//' @param scans a list of scans for output
//' @param ftFile a ft2 file's output path
//' @return void
//' @examples
//' header <- readFTheader("demo.ft2")
//' ft2 <- readAllScanMS2("demo.ft2")
//' writeAllScanMS2(header,ft2[1:10],"demo10.ft2")
//' @export
// [[Rcpp::export]]
bool writeAllScanMS2(List header, List scansList, CharacterVector ftFile)
{
    ftFileWriter writer(as<string>(ftFile),
                        as<string>(header["instrument"]),
                        as<string>(header["scanType"]),
                        as<string>(header["scanFilter"]),
                        as<bool>(header["hasCharge"]), 2);
    if (!writer.isWriteable || !writer.hasCharge)
        return false;
    vector<Scan> scans;
    scans.reserve(scansList.size());
    Scan scan;
    List scanList;
    DataFrame peaks;
    for (int i = 0; i < scansList.size(); i++)
    {
        scanList = scansList[i];
        scan.scanNumber = as<int>(scanList["scanNumber"]);
        scan.retentionTime = as<float>(scanList["retentionTime"]);
        scan.precursorScanNumber = as<int>(scanList["precursorScanNumber"]);
        scan.precursorMz = as<double>(scanList["precursorMz"]);
        scan.TIC = as<float>(scanList["TIC"]);
        scan.precursorCharge = as<int>(scanList["precursorCharge"]);

        peaks = scanList["peaks"];
        scan.mz = as<vector<double>>(peaks["mz"]);
        scan.intensity = as<vector<double>>(peaks["intensity"]);
        scan.resolution = as<vector<int>>(peaks["resolution"]);
        scan.baseLine = as<vector<float>>(peaks["baseLine"]);
        scan.signalToNoise = as<vector<float>>(peaks["signalToNoise"]);
        scan.charge = as<vector<int>>(peaks["charge"]);
        scans.push_back(scan);
    }
    // cout << "convert down" << endl;
    writer.writeAllScanMS2(scans);
    return true;
}