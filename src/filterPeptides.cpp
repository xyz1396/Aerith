#include "lib/peptidesFiltrator.h"
#include <Rcpp.h>
using namespace Rcpp;

//' getUnfilteredPeptides
//' @param workingPath a full path with .sip files in it
//' @return a dataframe of unique peptides and whether it is decoy sequence
//' @export
// [[Rcpp::export]]
DataFrame getUnfilteredPeptides(CharacterVector workingPath)
{
    sipFileReader reader(as<std::string>(workingPath));
    reader.readAllFiles();
    peptidesFiltrator filtrator(reader.sipPSMs, 0.01);
    filtrator.peptideMap.merge(filtrator.peptideMapCharge2);
    filtrator.peptideMap.merge(filtrator.peptideMapCharge3);
    filtrator.peptideMap.merge(filtrator.peptideMapChargeLargerThan3);
    size_t pepCount = filtrator.peptideMap.size();
    std::vector<std::string> identifiedPeptides(pepCount);
    std::vector<float> bestScores(pepCount);
    std::vector<bool> isDecoy(pepCount);
    size_t i = 0;
    for (auto pepIX : filtrator.peptideMap)
    {
        identifiedPeptides[i] = pepIX.first;
        bestScores[i] = pepIX.second.bestScore;
        isDecoy[i] = pepIX.second.isDecoy;
        i++;
    }
    return DataFrame::create(Named("identifiedPeptides") = move(identifiedPeptides),
                             _["bestScores"] = move(bestScores),
                             _["isDecoy"] = move(isDecoy));
}

//' getFilterThreshold
//' @param workingPath a full path with .sip files in it
//' @param OverallThreshold FDR thredhold of peptides
//' @return a dataframe about filter threshold and FDR results
//' @export
// [[Rcpp::export]]
DataFrame getFilterThreshold(CharacterVector workingPath, NumericVector OverallThreshold)
{
    sipFileReader reader(as<std::string>(workingPath));
    reader.readAllFiles();
    peptidesFiltrator filtrator(reader.sipPSMs, as<float>(OverallThreshold));
    filtrator.filterPeptideMap();
    std::vector<int> decoyCount{filtrator.decoyCountCharge2, filtrator.decoyCountCharge3,
                           filtrator.decoyCountChargeLargerThan3};
    std::vector<int> pepCount{filtrator.pepCountCharge2, filtrator.pepCountCharge3,
                         filtrator.pepCountChargeLargerThan3};
    std::vector<float> scoreThreshold{filtrator.scoreThresholdCharge2, filtrator.scoreThresholdCharge3,
                                 filtrator.scoreThresholdChargeLargerThan3};
    return DataFrame::create(
        Named("decoyCount") = decoyCount,
        _["pepCount"] = pepCount,
        _["scoreThreshold"] = scoreThreshold);
}

//' getFilterThresholdTopPSMs get filter threshold of top PSMs of each scan from multiple .sip file
//' @param workingPath a full path with .sip files in it
//' @param OverallThreshold FDR thredhold of peptides
//' @param topN store top N PSMs of each scan of one .FT file
//' @return a dataframe about filter threshold and FDR results
//' @export
// [[Rcpp::export]]
List getFilterThresholdTopPSMs(CharacterVector workingPath, NumericVector OverallThreshold, size_t topN)
{
    sipFileReader reader(as<std::string>(workingPath));
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    sipPSM topPSMs = reader.convertFilesScansTopPSMs();
    std::vector<sipPSM> topPSMss{topPSMs};
    DataFrame psmDf = DataFrame::create(Named("fileNames") = move(topPSMs.fileNames),
                                        _["scanNumbers"] = move(topPSMs.scanNumbers),
                                        _["parentCharges"] = move(topPSMs.parentCharges),
                                        _["measuredParentMasses"] = move(topPSMs.measuredParentMasses),
                                        _["calculatedParentMasses"] = move(topPSMs.calculatedParentMasses),
                                        _["searchNames"] = move(topPSMs.searchNames),
                                        _["scores"] = move(topPSMs.scores),
                                        _["identifiedPeptides"] = move(topPSMs.identifiedPeptides),
                                        _["originalPeptides"] = move(topPSMs.originalPeptides),
                                        _["proteinNames"] = move(topPSMs.proteinNames));
    peptidesFiltrator filtrator(topPSMss, as<float>(OverallThreshold));
    filtrator.filterPeptideMap();
    std::vector<int> decoyCount{filtrator.decoyCountCharge2, filtrator.decoyCountCharge3,
                           filtrator.decoyCountChargeLargerThan3};
    std::vector<int> pepCount{filtrator.pepCountCharge2, filtrator.pepCountCharge3,
                         filtrator.pepCountChargeLargerThan3};
    std::vector<float> scoreThreshold{filtrator.scoreThresholdCharge2, filtrator.scoreThresholdCharge3,
                                 filtrator.scoreThresholdChargeLargerThan3};
    return List::create(Named("threshold") = DataFrame::create(
                            Named("decoyCount") = decoyCount,
                            _["pepCount"] = pepCount,
                            _["scoreThreshold"] = scoreThreshold),
                        _["topPSMs"] = std::move(psmDf));
}