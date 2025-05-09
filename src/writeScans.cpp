#include "lib/ftFileWriter.h"
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rankyfify(NumericVector a)
{
    ftFileWriter writer;
    std::vector<double> b = as<std::vector<double>>(a);
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
    scan.mz = as<std::vector<double>>(peaks["mz"]);
    scan.intensity = as<std::vector<double>>(peaks["intensity"]);
    scan.resolution = as<std::vector<int>>(peaks["resolution"]);
    scan.baseLine = as<std::vector<float>>(peaks["baseLine"]);
    scan.signalToNoise = as<std::vector<float>>(peaks["signalToNoise"]);
    scan.charge = as<std::vector<int>>(peaks["charge"]);
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

//' write all MS1 scans has charge
//' @param header a list of FT file header
//' @param scans a list of scans for output
//' @param ftFile a ft1 file's output path
//' @return void
//' @examples
//' header <- readFTheader("demo.ft1")
//' ft1 <- readAllScanMS1("demo.ft1")
//' writeAllScanMS1(header,ft1[1:10],"demo10.ft1")
//' @export
// [[Rcpp::export]]
bool writeAllScanMS1(List header, List scansList, const String &ftFile)
{
    ftFileWriter writer(ftFile,
                        as<std::string>(header["instrument"]),
                        as<std::string>(header["scanType"]),
                        as<std::string>(header["scanFilter"]),
                        as<bool>(header["hasCharge"]), 1);
    if (!writer.isWriteable || !writer.hasCharge)
        return false;
    std::vector<Scan> scans;
    scans.reserve(scansList.size());
    Scan scan;
    List scanList;
    DataFrame peaks;
    for (int i = 0; i < scansList.size(); i++)
    {
        scanList = scansList[i];
        scan.scanNumber = as<int>(scanList["scanNumber"]);
        scan.retentionTime = as<float>(scanList["retentionTime"]);
        scan.TIC = as<float>(scanList["TIC"]);

        peaks = scanList["peaks"];
        scan.mz = as<std::vector<double>>(peaks["mz"]);
        scan.intensity = as<std::vector<double>>(peaks["intensity"]);
        scan.resolution = as<std::vector<int>>(peaks["resolution"]);
        scan.baseLine = as<std::vector<float>>(peaks["baseLine"]);
        scan.signalToNoise = as<std::vector<float>>(peaks["signalToNoise"]);
        scan.charge = as<std::vector<int>>(peaks["charge"]);
        scans.push_back(scan);
    }
    writer.writeAllScanMS1(scans);
    return true;
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
bool writeAllScanMS2(List header, List scansList, const String &ftFile)
{
    ftFileWriter writer(ftFile,
                        as<std::string>(header["instrument"]),
                        as<std::string>(header["scanType"]),
                        as<std::string>(header["scanFilter"]),
                        as<bool>(header["hasCharge"]), 2);
    if (!writer.isWriteable || !writer.hasCharge)
        return false;
    std::vector<Scan> scans;
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
        scan.isolationWindowCenterMZ = as<double>(scanList["isolationWindowCenterMZ"]);
        scan.TIC = as<float>(scanList["TIC"]);
        scan.precursorCharge = as<int>(scanList["precursorCharge"]);
        scan.precursorCharges = as<std::vector<int>>(scanList["precursorCharges"]);
        scan.precursorMZs = as<std::vector<double>>(scanList["precursorMZs"]);

        peaks = scanList["peaks"];
        scan.mz = as<std::vector<double>>(peaks["mz"]);
        scan.intensity = as<std::vector<double>>(peaks["intensity"]);
        scan.resolution = as<std::vector<int>>(peaks["resolution"]);
        scan.baseLine = as<std::vector<float>>(peaks["baseLine"]);
        scan.signalToNoise = as<std::vector<float>>(peaks["signalToNoise"]);
        scan.charge = as<std::vector<int>>(peaks["charge"]);
        scans.push_back(scan);
    }
    // cout << "convert down" << endl;
    writer.writeAllScanMS2(scans);
    return true;
}