#include <Rcpp.h>

#include "lib/ftFileWriter.h"

using namespace Rcpp;

//' rankify numeric vector via ftFileWriter
//' @param a A numeric vector whose values will be rank-transformed.
//' @return A numeric vector containing the ranks of the input values.
//' @examples
//' demo_vec <- c(12.5, 3.2, 7.7, 3.2)
//' rankyfify(demo_vec)
//' @export
// [[Rcpp::export]]
NumericVector rankyfify(NumericVector a) {
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
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' ft2 <- readAllScanMS2(demo_file)
//' plot(getRealScanFromList(ft2[["1346"]]))
//' ms2 <- denoiseOneMS2ScanHasCharge(ft2[["1346"]], 100, 10, 5)
//' plot(getRealScanFromList(ms2))
//' @export
// [[Rcpp::export]]
List denoiseOneMS2ScanHasCharge(List scanList, float window, float step,
                                float threshold) {
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
    DataFrame newPeaks = DataFrame::create(
        Named("mz") = scan.mz, _["intensity"] = scan.intensity,
        _["resolution"] = scan.resolution, _["baseLine"] = scan.baseLine,
        _["signalToNoise"] = scan.signalToNoise, _["charge"] = scan.charge);
    newScanList["peaks"] = newPeaks;
    return newScanList;
}

//' write all MS1 scans has charge
//' @param header a list of FT file header
//' @param scansList a list of scans for output
//' @param ftFile a ft1 file's output path
//' @return TRUE if the file was written successfully, FALSE otherwise
//' @examples
//' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
//' header <- readFTheader(demo_file)
//' ft1 <- readAllScanMS1(demo_file)
//' tmp <- tempdir()
//' writeAllScanMS1(header, ft1[1:10], file.path(tmp, "demo10.FT1"))
//' list.files(tmp, pattern = "demo10.FT1", full.names = TRUE)
//' @export
// [[Rcpp::export]]
bool writeAllScanMS1(List header, List scansList, const String &ftFile) {
    ftFileWriter writer(ftFile, as<std::string>(header["instrument"]),
                        as<std::string>(header["scanType"]),
                        as<std::string>(header["scanFilter"]),
                        as<bool>(header["hasCharge"]), 1);
    if (!writer.isWriteable || !writer.hasCharge) return false;
    std::vector<Scan> scans;
    scans.reserve(scansList.size());
    Scan scan;
    List scanList;
    DataFrame peaks;
    for (int i = 0; i < scansList.size(); i++) {
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
//' @param scansList a list of scans for output
//' @param ftFile a ft2 file's output path
//' @return TRUE if the file was written successfully, FALSE otherwise
//' @examples
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' header <- readFTheader(demo_file)
//' ft2 <- readAllScanMS2(demo_file)
//' tmp <- tempdir()
//' writeAllScanMS2(header,ft2[1:10],file.path(tmp, "demo10.FT2"))
//' list.files(tmp, pattern = "demo10.FT2", full.names = TRUE)
//' @export
// [[Rcpp::export]]
bool writeAllScanMS2(List header, List scansList, const String &ftFile) {
    ftFileWriter writer(ftFile, as<std::string>(header["instrument"]),
                        as<std::string>(header["scanType"]),
                        as<std::string>(header["scanFilter"]),
                        as<bool>(header["hasCharge"]), 2);
    if (!writer.isWriteable || !writer.hasCharge) return false;
    std::vector<Scan> scans;
    scans.reserve(scansList.size());
    Scan scan;
    List scanList;
    DataFrame peaks;
    for (int i = 0; i < scansList.size(); i++) {
        scanList = scansList[i];
        scan.scanNumber = as<int>(scanList["scanNumber"]);
        scan.retentionTime = as<float>(scanList["retentionTime"]);
        scan.precursorScanNumber = as<int>(scanList["precursorScanNumber"]);
        scan.isolationWindowCenterMZ =
            as<double>(scanList["isolationWindowCenterMZ"]);
        scan.TIC = as<float>(scanList["TIC"]);
        scan.precursorCharge = as<int>(scanList["precursorCharge"]);
        scan.precursorCharges =
            as<std::vector<int>>(scanList["precursorCharges"]);
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