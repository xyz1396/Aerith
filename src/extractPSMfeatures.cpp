#include "lib/PSMfeatureExtractor.h"
#include "lib/initSIP.h"
#include <Rcpp.h>
#include <regex>
#include <omp.h>
using namespace Rcpp;

//' extractPSMfeatures extract featueres of top PSMs from multiple .Spe2Pep.txt files
//' @param Spe2PepFilePath a full path with .Spe2Pep.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @param ftFilepath a full path with .FT1 and .FT2 files in it
//' @param ThreadNumber read ThreadNumber of FT file at the same time, it will increase ram usage
//' @return the PSMs in a dataframe in a list
//' @examples
//' psm <- extractPSMfeatures(Spe2PepFilePath, topN, ftFilepath, 3)
//' @export
// [[Rcpp::export]]

List extractPSMfeatures(String Spe2PepFilePath, int topN,
                        String ftFilepath, int ThreadNumber = 1)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    PSMfeatureExtractor extractor;
    extractor.extractPSMfeatureParallel(Spe2PepFilePath, topN, ftFilepath, ThreadNumber);
    std::vector<sipPSM> &sipPSMs = extractor.mSpe2PepFileReader.sipPSMs;
    List psmList(sipPSMs.size());
    List isotopicPeaksList;
    std::vector<double> MZs;
    std::vector<int> charges;
    std::vector<double> intensities;
    DataFrame psmDf;
    std::vector<std::string> PSMnames;
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        psmDf = DataFrame::create(Named("fileNames") = sipPSMs[i].fileNames,
                                  _["scanNumbers"] = sipPSMs[i].scanNumbers,
                                  _["precurosrScanNumbers"] = sipPSMs[i].precursorScanNumbers,
                                  _["parentCharges"] = sipPSMs[i].parentCharges,
                                  _["isolationWindowCenterMZs"] = sipPSMs[i].isolationWindowCenterMZs,
                                  _["measuredParentMasses"] = sipPSMs[i].measuredParentMasses,
                                  _["calculatedParentMasses"] = sipPSMs[i].calculatedParentMasses,
                                  _["searchNames"] = sipPSMs[i].searchNames,
                                  _["retentionTimes"] = sipPSMs[i].retentionTimes,
                                  _["MVHscores"] = sipPSMs[i].MVHscores,
                                  _["XcorrScores"] = sipPSMs[i].XcorrScores,
                                  _["WDPscores"] = sipPSMs[i].WDPscores,
                                  _["ranks"] = sipPSMs[i].ranks,
                                  _["identifiedPeptides"] = sipPSMs[i].identifiedPeptides,
                                  _["originalPeptides"] = sipPSMs[i].originalPeptides,
                                  _["proteinNames"] = sipPSMs[i].proteinNames,
                                  _["istopicPeakNumbers"] = sipPSMs[i].istopicPeakNumbers,
                                  _["MS1IsotopicAbundances"] = sipPSMs[i].MS1IsotopicAbundances);
        isotopicPeaksList = List(sipPSMs[i].isotopicPeakss.size());
        PSMnames.clear();
        PSMnames.reserve(sipPSMs[i].isotopicPeakss.size());
        for (size_t j = 0; j < sipPSMs[i].isotopicPeakss.size(); j++)
        {
            MZs.clear();
            charges.clear();
            intensities.clear();
            for (size_t k = 0; k < sipPSMs[i].isotopicPeakss[j].size(); k++)
            {
                MZs.push_back(sipPSMs[i].isotopicPeakss[j][k].mz);
                charges.push_back(sipPSMs[i].isotopicPeakss[j][k].charge);
                intensities.push_back(sipPSMs[i].isotopicPeakss[j][k].intensity);
            }
            isotopicPeaksList[j] = DataFrame::create(Named("MZs") = MZs,
                                                     _("Charges") = charges,
                                                     _("Intensities") = intensities);
            PSMnames.push_back(to_string(sipPSMs[i].scanNumbers[j]) + "\t" +
                               to_string(sipPSMs[i].precursorScanNumbers[j]) + "\t" + sipPSMs[i].identifiedPeptides[j]);
        }
        isotopicPeaksList.names() = PSMnames;
        psmList[i] = List::create(Named("PSMdataframe") = psmDf,
                                  _("isotopicPeakList") = isotopicPeaksList);
    }
    psmList.names() = extractor.mSpe2PepFileReader.FT2s;
    return psmList;
}