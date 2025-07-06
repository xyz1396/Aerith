#include "lib/ftFileReader.h"
#include <Rcpp.h>

using namespace Rcpp;

//' read FT file header
//' @param ftFile a ft1 file's full path
//' @return a list of ft file header
//' @examples
//' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
//' header <- readFTheader(demo_file)
//' @export
// [[Rcpp::export]]
List readFTheader(String ftFile)
{
    ftFileReader reader(ftFile);
    reader.detectPrecursorAndCharge();
    List header = List::create(Named("instrument") = reader.instrument,
                               _["scanType"] = reader.scanType,
                               _["scanFilter"] = reader.scanFilter,
                               _["hasPrecursor"] = reader.hasPrecursor,
                               _["hasCharge"] = reader.hasCharge);
    return header;
}

//' read MS1 scans with scanNumber as index in a range
//' @param ftFile a ft1 file's full path
//' @param startScanNumber read scans starting from this scanNumber
//' @param endScanNumber read scans ending at this scanNumber
//' @return a list of MS1 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
//' ft1 <- readScansMS1(demo_file, 1398, 1503)
//' @export
// [[Rcpp::export]]
List readScansMS1(const String ftFile, const size_t startScanNumber, const size_t endScanNumber)
{
    ftFileReader reader(ftFile);
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    reader.readScans(startScanNumber, endScanNumber);
    List scanList(reader.Scans.size());
    Scan *mScan;
    CharacterVector scanNumbers(reader.Scans.size());
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = std::move(mScan->mz),
                                             _["intensity"] = std::move(mScan->intensity),
                                             _["resolution"] = std::move(mScan->resolution),
                                             _["baseLine"] = std::move(mScan->baseLine),
                                             _["signalToNoise"] = std::move(mScan->signalToNoise),
                                             _["charge"] = std::move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["TIC"] = mScan->TIC,
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbers[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}

//' read MS1 scans with scanNumber as index in a vector
//' @param ftFile a ft1 file's full path
//' @param startScanNumber read scans starting from this scanNumber
//' @param endScanNumber read scans ending at this scanNumber
//' @return a list of MS1 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
//' ft1 <- readScansMS1Vector(demo_file, c(1398, 1503, 1508))
//' @export
// [[Rcpp::export]]
List readScansMS1Vector(const String ftFile, const NumericVector scanNumbersVector)
{
    ftFileReader reader(ftFile);
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    std::vector<std::size_t> scanNumbers = as<std::vector<std::size_t>>(scanNumbersVector);
    reader.readScans(scanNumbers);
    List scanList(reader.Scans.size());
    Scan *mScan;
    CharacterVector scanNumbersStr(reader.Scans.size());
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = std::move(mScan->mz),
                                             _["intensity"] = std::move(mScan->intensity),
                                             _["resolution"] = std::move(mScan->resolution),
                                             _["baseLine"] = std::move(mScan->baseLine),
                                             _["signalToNoise"] = std::move(mScan->signalToNoise),
                                             _["charge"] = std::move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["TIC"] = mScan->TIC,
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbersStr[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbersStr;
    return scanList;
}

//' read MS1 scans with scanNumber as index
//' @param ftFile a ft1 file's full path
//' @return a list of MS1 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT1", package = "Aerith")
//' ft1 <- readAllScanMS1(demo_file)
//' @export
// [[Rcpp::export]]
List readAllScanMS1(const String ftFile)
{
    ftFileReader reader(ftFile);
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
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbers[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}

//' read MS2 scans with scanNumber as index in a range
//' @param ftFile a ft2 file's full path
//' @param startScanNumber read scans starting from this scanNumber
//' @param endScanNumber read scans ending at this scanNumber
//' @return a list of MS2 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' ft2 <- readScansMS2(demo_file, 1350, 1355)
//' @export
// [[Rcpp::export]]
List readScansMS2(const String ftFile, const size_t startScanNumber, const size_t endScanNumber)
{
    ftFileReader reader(ftFile);
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    reader.readScans(startScanNumber, endScanNumber);
    List scanList(reader.Scans.size());
    Scan *mScan;
    CharacterVector scanNumbers(reader.Scans.size());
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = std::move(mScan->mz),
                                             _["intensity"] = std::move(mScan->intensity),
                                             _["resolution"] = std::move(mScan->resolution),
                                             _["baseLine"] = std::move(mScan->baseLine),
                                             _["signalToNoise"] = std::move(mScan->signalToNoise),
                                             _["charge"] = std::move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["isolationWindowCenterMZ"] = mScan->isolationWindowCenterMZ,
                                      _["TIC"] = mScan->TIC,
                                      _["precursorCharges"] = mScan->precursorCharges,
                                      _["precursorMZs"] = mScan->precursorMZs,
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbers[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}

//' read MS2 scans with scanNumber as index in a vector
//' @param ftFile a ft2 file's full path
//' @param scanNumbersVector read scans starting of these scanNumbers
//' @return a list of MS2 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' ft2 <- readScansMS2Vector(demo_file, c(1350, 1355, 1359))
//' @export
// [[Rcpp::export]]
List readScansMS2Vector(const String ftFile, const NumericVector scanNumbersVector)
{
    ftFileReader reader(ftFile);
    // avoid empty scanList crash
    if (reader.isEmpty)
        return 0;
    std::vector<std::size_t> scanNumbers = as<std::vector<std::size_t>>(scanNumbersVector);
    reader.readScans(scanNumbers);
    List scanList(reader.Scans.size());
    Scan *mScan;
    CharacterVector scanNumbersStr(reader.Scans.size());
    for (size_t i = 0; i < reader.Scans.size(); i++)
    {
        mScan = &reader.Scans[i];
        DataFrame peakDf = DataFrame::create(Named("mz") = std::move(mScan->mz),
                                             _["intensity"] = std::move(mScan->intensity),
                                             _["resolution"] = std::move(mScan->resolution),
                                             _["baseLine"] = std::move(mScan->baseLine),
                                             _["signalToNoise"] = std::move(mScan->signalToNoise),
                                             _["charge"] = std::move(mScan->charge));
        List mScanList = List::create(Named("scanNumber") = mScan->scanNumber,
                                      _["retentionTime"] = mScan->retentionTime,
                                      _["precursorScanNumber"] = mScan->precursorScanNumber,
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["isolationWindowCenterMZ"] = mScan->isolationWindowCenterMZ,
                                      _["TIC"] = mScan->TIC,
                                      _["precursorCharges"] = mScan->precursorCharges,
                                      _["precursorMZs"] = mScan->precursorMZs,
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbersStr[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbersStr;
    return scanList;
}

//' read MS2 scans with scanNumber as index
//' @param ftFile a ft2 file's full path
//' @return a list of MS2 scans with names of scan number
//' @examples
//' demo_file <- system.file("extdata", "demo.FT2", package = "Aerith")
//' ft2 <- readAllScanMS2(demo_file)
//' @export
// [[Rcpp::export]]
List readAllScanMS2(const String ftFile)
{
    ftFileReader reader(ftFile);
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
                                      _["precursorCharge"] = mScan->precursorCharge,
                                      _["isolationWindowCenterMZ"] = mScan->isolationWindowCenterMZ,
                                      _["TIC"] = mScan->TIC,
                                      _["precursorCharges"] = mScan->precursorCharges,
                                      _["precursorMZs"] = mScan->precursorMZs,
                                      _["peaks"] = std::move(peakDf));
        scanList[i] = std::move(mScanList);
        scanNumbers[i] = std::to_string(mScan->scanNumber);
    }
    // change list's names to scanNumber
    scanList.names() = scanNumbers;
    return scanList;
}
