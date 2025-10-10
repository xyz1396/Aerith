#include "lib/sipFileReader.h"
#include <Rcpp.h>
#include <regex>
using namespace Rcpp;

//' readSip
//' @param sipFile a .sip file's full path
//' @description Read a single .sip file and convert its peptide-spectrum matches (PSMs) into an R list.
//' @details The reader parses the provided .sip file and returns metadata together with a data.frame of PSM-level attributes.
//' @param sipFile A character vector of length one containing the full path to a .sip file.
//' @return An R list with file-level metadata (`fileName`, `scanType`, `searchName`, `scoringFunction`) and a `PSM` data frame containing scan numbers, charges, masses, scores, ranks, peptides, and protein names.
//' @examples
//' demo_file <- system.file("extdata", "demo.sip", package = "Aerith")
//' re <- readSip(demo_file)
//' head(re$PSM)
//' @export
// [[Rcpp::export]]
List readSip(CharacterVector sipFile)
{
    sipFileReader reader;
    reader.readOneFile(as<std::string>(sipFile));
    DataFrame psmDf = DataFrame::create(Named("scanNumbers") = reader.currentSipPSM.scanNumbers,
                                        _["parentCharges"] = reader.currentSipPSM.parentCharges,
                                        _["measuredParentMasses"] = reader.currentSipPSM.measuredParentMasses,
                                        _["calculatedParentMasses"] = reader.currentSipPSM.calculatedParentMasses,
                                        _["scores"] = reader.currentSipPSM.scores,
                                        _["ranks"] = reader.currentSipPSM.ranks,
                                        _["identifiedPeptides"] = reader.currentSipPSM.identifiedPeptides,
                                        _["originalPeptides"] = reader.currentSipPSM.originalPeptides,
                                        _["proteinNames"] = reader.currentSipPSM.proteinNames);
    List mSipList = List::create(Named("fileName") = reader.currentSipPSM.fileName,
                                 _["scanType"] = reader.currentSipPSM.scanType,
                                 _["searchName"] = reader.currentSipPSM.searchName,
                                 _["scoringFunction"] = reader.currentSipPSM.scoringFunction,
                                 _["PSM"] = psmDf);
    return mSipList;
}

//' readSips
//' @description Read every `.sip` file found in the provided directory and transform each file's peptide-spectrum matches into an R list entry.
//' @details This function constructs a `sipFileReader` for the supplied folder, loads all `.sip` files, and for each file returns a list containing metadata and a `data.frame` of peptide-spectrum matches (PSMs). The resulting object is an R list whose elements correspond to individual input files.
//' @param workingPath Character vector of length one giving the directory that contains `.sip` files.
//' @return An R list; each element represents one `.sip` file and provides file-level descriptors (`fileName`, `scanType`, `searchName`, `scoringFunction`) together with a `PSM` data frame holding scan numbers, precursor charges, observed/calculated masses, scores, ranks, peptide sequences, and protein assignments.
//' @examples
//' demo_dir <- system.file("extdata", package = "Aerith")
//' re <- readSips(demo_dir)
//' head(re[[1]]$PSM)
//' @export
// [[Rcpp::export]]
List readSips(CharacterVector workingPath)
{
    sipFileReader reader(as<std::string>(workingPath));
    reader.readAllFiles();
    List psmList(reader.sipPSMs.size());
    sipPSM *mPSM;
    for (size_t i = 0; i < reader.sipPSMs.size(); i++)
    {
        mPSM = &reader.sipPSMs[i];
        // use std::move to speed up copy vector
        DataFrame psmDf = DataFrame::create(Named("scanNumbers") = std::move(mPSM->scanNumbers),
                                            _["parentCharges"] = std::move(mPSM->parentCharges),
                                            _["measuredParentMasses"] = std::move(mPSM->measuredParentMasses),
                                            _["calculatedParentMasses"] = std::move(mPSM->calculatedParentMasses),
                                            _["scores"] = std::move(mPSM->scores),
                                            _["ranks"] = std::move(mPSM->ranks),
                                            _["identifiedPeptides"] = std::move(mPSM->identifiedPeptides),
                                            _["originalPeptides"] = std::move(mPSM->originalPeptides),
                                            _["proteinNames"] = std::move(mPSM->proteinNames));
        List mSipList = List::create(Named("fileName") = mPSM->fileName,
                                     _["scanType"] = mPSM->scanType,
                                     _["searchName"] = mPSM->searchName,
                                     _["scoringFunction"] = mPSM->scoringFunction,
                                     _["PSM"] = psmDf);
        psmList[i] = std::move(mSipList);
    }
    return psmList;
}

//' readFilesScansTopPSMs
//' @description Aggregate the top-ranked peptide-spectrum matches (PSMs) from each scan across all `.sip` files found in a directory.
//' @details This helper constructs an internal `sipFileReader`, loads every `.sip` file located in the supplied path, and extracts the best `topN` PSMs per scan. The result is returned as a single `data.frame` with file identifiers, scan numbers, charge states, mass measurements, scores, peptide assignments, and protein annotations.
//' @param workingPath Character vector of length one giving the directory that contains `.sip` files.
//' @param topN Integer specifying how many top-ranked PSMs per scan to retain for each `.sip` file.
//' @return An R `data.frame` with one row per retained PSM and the columns `fileNames`, `scanNumbers`, `parentCharges`, `measuredParentMasses`, `calculatedParentMasses`, `searchNames`, `scores`, `identifiedPeptides`, `originalPeptides`, and `proteinNames`.
//' @examples
//' demo_dir <- system.file("extdata", package = "Aerith")
//' re <- readFilesScansTopPSMs(demo_dir, 10)
//' head(re)
//' @export
// [[Rcpp::export]]
DataFrame readFilesScansTopPSMs(CharacterVector workingPath, size_t topN)
{
    sipFileReader reader(as<std::string>(workingPath));
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    DataFrame psmDf = DataFrame::create(Named("fileNames") = std::move(topPSMs.fileNames),
                                        _["scanNumbers"] = std::move(topPSMs.scanNumbers),
                                        _["parentCharges"] = std::move(topPSMs.parentCharges),
                                        _["measuredParentMasses"] = std::move(topPSMs.measuredParentMasses),
                                        _["calculatedParentMasses"] = std::move(topPSMs.calculatedParentMasses),
                                        _["searchNames"] = std::move(topPSMs.searchNames),
                                        _["scores"] = std::move(topPSMs.scores),
                                        _["identifiedPeptides"] = std::move(topPSMs.identifiedPeptides),
                                        _["originalPeptides"] = std::move(topPSMs.originalPeptides),
                                        _["proteinNames"] = std::move(topPSMs.proteinNames));
    return psmDf;
}

//' readFilesScansTopPSMsFromOneFT2 read each scan's top PSMs from multiple .sip files of one .FT2 file
//' @param workingPath a full path with .sip files in it
//' @param pattern a regex pattern of the .FT2 file
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @return a dataframe of top N PSMs
//' @examples
//' demo_dir <- system.file("extdata", package = "Aerith")
//' re <- readFilesScansTopPSMsFromOneFT2(demo_dir, ".*demo.*", 3)
//' head(re)
//' @export
// [[Rcpp::export]]
DataFrame readFilesScansTopPSMsFromOneFT2(String workingPath, String pattern, size_t topN)
{
    sipFileReader reader(workingPath);
    std::vector<std::string> matchedNames;
    std::regex ePattern((std::string)pattern);
    for (size_t i = 0; i < reader.sipFileNames.size(); i++)
    {
        if (regex_match(reader.sipFileNames[i], ePattern))
            matchedNames.push_back(reader.sipFileNames[i]);
    }
    if (matchedNames.size() == 0)
    {
        Rcout << "No .sip file was matched!" << std::endl;
        return DataFrame();
    }
    reader.sipFileNames = matchedNames;
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    DataFrame psmDf = DataFrame::create(Named("FileName") = std::move(topPSMs.fileNames),
                                        _["ScanNumber"] = std::move(topPSMs.scanNumbers),
                                        _["ParentCharge"] = std::move(topPSMs.parentCharges),
                                        _["MeasuredParentMass"] = std::move(topPSMs.measuredParentMasses),
                                        _["CalculatedParentMass"] = std::move(topPSMs.calculatedParentMasses),
                                        _["SearchName"] = std::move(topPSMs.searchNames),
                                        _["Rank"] = std::move(topPSMs.ranks),
                                        _["Score"] = std::move(topPSMs.scores),
                                        _["IdentifiedPeptide"] = std::move(topPSMs.identifiedPeptides),
                                        _["OriginalPeptide"] = std::move(topPSMs.originalPeptides),
                                        _["ProteinNames"] = std::move(topPSMs.proteinNames));
    return psmDf;
}