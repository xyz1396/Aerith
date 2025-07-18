#include "lib/Spe2PepFileReader.h"
#include <Rcpp.h>
#include <regex>
using namespace Rcpp;

//' readSpe2Pep
//' @param Spe2PepFile a .Spe2PepFile file's full path
//' @return the PSMs in a dataframe in a list
//' @examples
//' target_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' psm <- readSpe2Pep(target_file)
//' psm <- psm$PSM
//' @export
// [[Rcpp::export]]
List readSpe2Pep(String Spe2PepFile)
{
    Spe2PepFileReader reader;
    reader.readOneEntireFile(Spe2PepFile);
    DataFrame psmDf = DataFrame::create(Named("scanNumbers") = reader.currentSipPSM.scanNumbers,
                                        _["precursorScanNumbers"] = reader.currentSipPSM.precursorScanNumbers,
                                        _["parentCharges"] = reader.currentSipPSM.parentCharges,
                                        _["isolationWindowCenterMZs"] = reader.currentSipPSM.isolationWindowCenterMZs,
                                        _["measuredParentMasses"] = reader.currentSipPSM.measuredParentMasses,
                                        _["calculatedParentMasses"] = reader.currentSipPSM.calculatedParentMasses,
                                        _["retentionTimes"] = reader.currentSipPSM.retentionTimes,
                                        _["WDPscores"] = reader.currentSipPSM.WDPscores,
                                        _["XcorrScores"] = reader.currentSipPSM.XcorrScores,  
                                        _["MVHscores"] = reader.currentSipPSM.MVHscores,                                                                          
                                        _["ranks"] = reader.currentSipPSM.ranks,
                                        _["identifiedPeptides"] = reader.currentSipPSM.identifiedPeptides,
                                        _["originalPeptides"] = reader.currentSipPSM.originalPeptides,
                                        _["nakePeptides"] = reader.currentSipPSM.nakePeptides,
                                        _["proteinNames"] = reader.currentSipPSM.proteinNames);
    List mSipList = List::create(Named("fileName") = reader.currentSipPSM.fileName,
                                 _["scanType"] = reader.currentSipPSM.scanType,
                                 _["searchName"] = reader.currentSipPSM.searchName,
                                 _["PSM"] = psmDf);
    return mSipList;
}

//' readSpe2Peps
//' @param workingPath a full path with .Spe2Pep.txt files in it
//' @return the PSMs dataframes in lists
//' @examples
//' tmp <- tempdir()
//' sip_dir <- file.path(tmp, "sip")
//' dir.create(sip_dir)
//' demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
//' demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
//' list.files(sip_dir, full.names = TRUE)
//' psm <- readSpe2Peps(sip_dir)
//' psm <- psm[[1]]$PSM
//' @export
// [[Rcpp::export]]
List readSpe2Peps(String workingPath)
{
    Spe2PepFileReader reader(workingPath);
    reader.readAllFiles();
    List psmList(reader.sipPSMs.size());
    sipPSM *mPSM;
    for (size_t i = 0; i < reader.sipPSMs.size(); i++)
    {
        mPSM = &reader.sipPSMs[i];
        // use std::move to speed up copy vector
        DataFrame psmDf = DataFrame::create(Named("scanNumbers") = std::move(mPSM->scanNumbers),
                                            _["precursorScanNumbers"] = std::move(mPSM->precursorScanNumbers),
                                            _["parentCharges"] = std::move(mPSM->parentCharges),
                                            _["isolationWindowCenterMZs"] = std::move(mPSM->isolationWindowCenterMZs),
                                            _["measuredParentMasses"] = std::move(mPSM->measuredParentMasses),
                                            _["calculatedParentMasses"] = std::move(mPSM->calculatedParentMasses),
                                            _["retentionTimes"] = std::move(mPSM->retentionTimes),
                                            _["MVHscores"] = std::move(mPSM->MVHscores),
                                            _["XcorrScores"] = std::move(mPSM->XcorrScores),
                                            _["WDPscores"] = std::move(mPSM->WDPscores),
                                            _["ranks"] = std::move(mPSM->ranks),
                                            _["identifiedPeptides"] = std::move(mPSM->identifiedPeptides),
                                            _["originalPeptides"] = std::move(mPSM->originalPeptides),
                                            _["nakePeptides"] =  std::move(mPSM->nakePeptides),
                                            _["proteinNames"] = std::move(mPSM->proteinNames));
        List mSipList = List::create(Named("fileName") = mPSM->fileName,
                                     _["scanType"] = mPSM->scanType,
                                     _["searchName"] = mPSM->searchName,
                                     _["PSM"] = psmDf);
        psmList[i] = std::move(mSipList);
    }
    return psmList;
}

//' readSpe2PepFilesScansTopPSMs read each scan's top PSMs from multiple .Spe2Pep.txt files
//' @param workingPath a full path with .Spe2Pep.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @return the PSMs in a dataframe in a list
//' @examples
//' tmp <- tempdir()
//' sip_dir <- file.path(tmp, "sip")
//' dir.create(sip_dir)
//' demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
//' demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
//' list.files(sip_dir, full.names = TRUE)
//' psm <- readSpe2PepFilesScansTopPSMs(sip_dir, 3)
//' @export
// [[Rcpp::export]]
DataFrame readSpe2PepFilesScansTopPSMs(String workingPath, size_t topN = 5)
{
    Spe2PepFileReader reader(workingPath);
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    DataFrame psmDf = DataFrame::create(Named("fileNames") = std::move(topPSMs.fileNames),
                                        _["scanNumbers"] = std::move(topPSMs.scanNumbers),
                                        _["precursorScanNumbers"] = std::move(topPSMs.precursorScanNumbers),
                                        _["parentCharges"] = std::move(topPSMs.parentCharges),
                                        _["isolationWindowCenterMZs"] = std::move(topPSMs.isolationWindowCenterMZs),
                                        _["measuredParentMasses"] = std::move(topPSMs.measuredParentMasses),
                                        _["calculatedParentMasses"] = std::move(topPSMs.calculatedParentMasses),
                                        _["searchNames"] = std::move(topPSMs.searchNames),
                                        _["retentionTimes"] = std::move(topPSMs.retentionTimes),
                                        _["WDPscores"] = std::move(topPSMs.WDPscores),
                                        _["XcorrScores"] = std::move(topPSMs.XcorrScores), 
                                        _["MVHscores"] = std::move(topPSMs.MVHscores),                                                                               
                                        _["ranks"] = std::move(topPSMs.ranks),
                                        _["identifiedPeptides"] = std::move(topPSMs.identifiedPeptides),
                                        _["originalPeptides"] = std::move(topPSMs.originalPeptides),
                                        _["nakePeptides"] = std::move(topPSMs.nakePeptides),
                                        _["proteinNames"] = std::move(topPSMs.proteinNames));
    return psmDf;
}

//' readSpe2PepFilesScansTopPSMsFromOneFT2 read each scan's top PSMs from multiple .Spe2PepFile.txt files of one .FT2 file
//' @param workingPath a full path with .Spe2PepFile.txt files in it
//' @param pattern a regex pattern of the .FT2 file
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @return a dataframe of top N PSMs
//' @examples
//' tmp <- tempdir()
//' sip_dir <- file.path(tmp, "sip")
//' dir.create(sip_dir)
//' demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
//' demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
//' list.files(sip_dir, full.names = TRUE)
//' top3 <-  readSpe2PepFilesScansTopPSMsFromOneFT2(sip_dir, ".*X13.*", 3)
//' @export
// [[Rcpp::export]]
DataFrame readSpe2PepFilesScansTopPSMsFromOneFT2(String workingPath, String pattern, size_t topN = 5)
{
    Spe2PepFileReader reader(workingPath);
    std::vector<std::string> matchedNames;
    std::regex ePattern((std::string)pattern);
    for (size_t i = 0; i < reader.sipFileNames.size(); i++)
    {
        if (regex_match(reader.sipFileNames[i], ePattern))
            matchedNames.push_back(reader.sipFileNames[i]);
    }
    if (matchedNames.size() == 0)
    {
        Rcout << "No .Spe2PepFile file was matched!" << std::endl;
        return DataFrame();
    }
    reader.sipFileNames = matchedNames;
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    DataFrame psmDf = DataFrame::create(Named("fileNames") = std::move(topPSMs.fileNames),
                                        _["scanNumbers"] = std::move(topPSMs.scanNumbers),
                                        _["precursorScanNumbers"] = std::move(topPSMs.precursorScanNumbers),
                                        _["parentCharges"] = std::move(topPSMs.parentCharges),
                                        _["isolationWindowCenterMZs"] = std::move(topPSMs.isolationWindowCenterMZs),
                                        _["measuredParentMasses"] = std::move(topPSMs.measuredParentMasses),
                                        _["calculatedParentMasses"] = std::move(topPSMs.calculatedParentMasses),
                                        _["searchNames"] = std::move(topPSMs.searchNames),
                                        _["retentionTimes"] = std::move(topPSMs.retentionTimes),
                                        _["MVHscores"] = std::move(topPSMs.MVHscores),
                                        _["XcorrScores"] = std::move(topPSMs.XcorrScores),
                                        _["WDPscores"] = std::move(topPSMs.WDPscores),
                                        _["ranks"] = std::move(topPSMs.ranks),
                                        _["identifiedPeptides"] = std::move(topPSMs.identifiedPeptides),
                                        _["originalPeptides"] = std::move(topPSMs.originalPeptides),
                                        _["nakePeptides"] = std::move(topPSMs.nakePeptides),
                                        _["proteinNames"] = std::move(topPSMs.proteinNames));
    return psmDf;
}

//' readSpe2PepFilesScansTopPSMsFromEachFT2Parallel read each scan's top PSMs from multiple .Spe2PepFile.txt files of each .FT2 file
//' @param workingPath a full path with .Spe2PepFile.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @return a dataframe of top N PSMs
//' @examples
//' tmp <- tempdir()
//' sip_dir <- file.path(tmp, "sip")
//' dir.create(sip_dir)
//' demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
//' demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
//' list.files(sip_dir, full.names = TRUE)
//' top3 <-  readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(sip_dir, 3)
//' @export
// [[Rcpp::export]]
List readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(String workingPath, size_t topN = 5)
{
    Spe2PepFileReader reader;
    reader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(workingPath, topN);
    List psmList(reader.sipPSMs.size());
    for (size_t i = 0; i < reader.sipPSMs.size(); i++)
    {
        DataFrame psmDf = DataFrame::create(Named("fileNames") = reader.sipPSMs[i].fileNames,
                                            _["scanNumbers"] = reader.sipPSMs[i].scanNumbers,
                                            _["precursorScanNumbers"] = reader.sipPSMs[i].precursorScanNumbers,
                                            _["parentCharges"] = reader.sipPSMs[i].parentCharges,
                                            _["isolationWindowCenterMZs"] = reader.sipPSMs[i].isolationWindowCenterMZs,
                                            _["measuredParentMasses"] = reader.sipPSMs[i].measuredParentMasses,
                                            _["calculatedParentMasses"] = reader.sipPSMs[i].calculatedParentMasses,
                                            _["searchNames"] = reader.sipPSMs[i].searchNames,
                                            _["retentionTimes"] = reader.sipPSMs[i].retentionTimes,
                                            _["WDPscores"] = reader.sipPSMs[i].WDPscores,
                                            _["XcorrScores"] = reader.sipPSMs[i].XcorrScores,
                                            _["MVHscores"] = reader.sipPSMs[i].MVHscores,
                                            _["ranks"] = reader.sipPSMs[i].ranks,
                                            _["identifiedPeptides"] = reader.sipPSMs[i].identifiedPeptides,
                                            _["originalPeptides"] = reader.sipPSMs[i].originalPeptides,
                                            _["nakePeptides"] = reader.sipPSMs[i].nakePeptides,
                                            _["proteinNames"] = reader.sipPSMs[i].proteinNames);

        psmList[i] = psmDf;
    }
    psmList.names() = reader.FT2s;
    return psmList;
}

//' readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParalle read each scan's top PSMs from multiple .Spe2PepFile.txt files of each .FT2 file
//' @param targetPath a full path with target .Spe2PepFile.txt files in it
//' @param decoyPath a full path with decoy .Spe2PepFile.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @return a dataframe of top N PSMs
//' @examples
//' tmp <- tempdir()
//' target_dir <- file.path(tmp, "target")
//' dir.create(target_dir, showWarnings = FALSE)
//' target_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(target_file, file.path(target_dir, "Pan_052322_X13.SIP_C13_050_000Pct.Spe2Pep.txt"))
//' decoy_dir <- file.path(tmp, "decoy")
//' dir.create(decoy_dir, showWarnings = FALSE)
//' decoy_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(decoy_file, file.path(decoy_dir, "Pan_052322_X13.SIP_C13_050_000Pct.Spe2Pep.txt"))
//' list.files(target_dir, full.names = TRUE)
//' list.files(decoy_dir, full.names = TRUE)
//' top3 <- readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(target_dir, decoy_dir, 3)
//' @export
// [[Rcpp::export]]
List readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(String targetPath, String decoyPath, size_t topN = 5)
{
    Spe2PepFileReader reader;
    reader.readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(targetPath, decoyPath, topN);
    List psmList(reader.sipPSMs.size());
    for (size_t i = 0; i < reader.sipPSMs.size(); i++)
    {
        DataFrame psmDf = DataFrame::create(Named("fileNames") = reader.sipPSMs[i].fileNames,
                                            _["isDecoys"] = reader.sipPSMs[i].isDecoys,
                                            _["scanNumbers"] = reader.sipPSMs[i].scanNumbers,
                                            _["precursorScanNumbers"] = reader.sipPSMs[i].precursorScanNumbers,
                                            _["parentCharges"] = reader.sipPSMs[i].parentCharges,
                                            _["isolationWindowCenterMZs"] = reader.sipPSMs[i].isolationWindowCenterMZs,
                                            _["measuredParentMasses"] = reader.sipPSMs[i].measuredParentMasses,
                                            _["calculatedParentMasses"] = reader.sipPSMs[i].calculatedParentMasses,
                                            _["searchNames"] = reader.sipPSMs[i].searchNames,
                                            _["retentionTimes"] = reader.sipPSMs[i].retentionTimes,
                                            _["WDPscores"] = reader.sipPSMs[i].WDPscores,                                            
                                            _["XcorrScores"] = reader.sipPSMs[i].XcorrScores,
                                            _["MVHscores"] = reader.sipPSMs[i].MVHscores,
                                            _["ranks"] = reader.sipPSMs[i].ranks,
                                            _["identifiedPeptides"] = reader.sipPSMs[i].identifiedPeptides,
                                            _["originalPeptides"] = reader.sipPSMs[i].originalPeptides,
                                            _["nakePeptides"] = reader.sipPSMs[i].nakePeptides,
                                            _["proteinNames"] = reader.sipPSMs[i].proteinNames);

        psmList[i] = psmDf;
    }
    psmList.names() = reader.FT2s;
    return psmList;
}

//' writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel read each scan's top PSMs from multiple .Spe2PepFile.txt
//' files of each .FT2 file and write them to one tsv file
//' @param workingPath a full path with .Spe2PepFile.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @param fileName the output path
//' @return nothing but write a tsv of top N PSMs
//' @examples
//' tmp <- tempdir()
//' sip_dir <- file.path(tmp, "sip")
//' dir.create(sip_dir)
//' demo_file <- system.file("extdata", "demo_target.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000target.Spe2Pep.txt"))
//' demo_file <- system.file("extdata", "demo_decoy.Spe2Pep.txt", package = "Aerith")
//' file.copy(demo_file, file.path(sip_dir, "Pan_052322_X13.SIP_C13_050_000decoy.Spe2Pep.txt"))
//' writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(sip_dir, 3, file.path(sip_dir, "top3.tsv"))
//' list.files(sip_dir, full.names = TRUE)
//' @export
// [[Rcpp::export]]
void writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(String workingPath, size_t topN = 5, String fileName = "a.tsv")
{
    Spe2PepFileReader reader;
    reader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(workingPath, topN);
    reader.writeTSV(fileName);
}