#include "lib/PSMfeatureExtractor.h"
#include "lib/initSIP.h"
#include <Rcpp.h>
#include <regex>
#include <omp.h>
using namespace Rcpp;

// for create data.frame with more than 20 columns
class ListBuilder
{

public:
    ListBuilder(){};
    ~ListBuilder(){};

    inline ListBuilder &add(std::string const &name, SEXP x)
    {
        names.push_back(name);

        // NOTE: we need to protect the SEXPs we pass in; there is
        // probably a nicer way to handle this but ...
        elements.push_back(PROTECT(x));

        return *this;
    }

    inline operator List() const
    {
        List result(elements.size());
        for (size_t i = 0; i < elements.size(); ++i)
        {
            result[i] = elements[i];
        }
        result.attr("names") = wrap(names);
        UNPROTECT(elements.size());
        return result;
    }

    inline operator DataFrame() const
    {
        List result = static_cast<List>(*this);
        result.attr("class") = "data.frame";
        result.attr("row.names") = IntegerVector::create(NA_INTEGER, XLENGTH(elements[0]));
        return result;
    }

private:
    std::vector<std::string> names;
    std::vector<SEXP> elements;

    ListBuilder(ListBuilder const &){}; // not safe to copy
};

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
                        String ftFilepath, int ThreadNumber = 3)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    PSMfeatureExtractor extractor;
    extractor.extractPSMfeatureParallel(Spe2PepFilePath, topN, ftFilepath, ThreadNumber);
    std::vector<sipPSM> &sipPSMs = extractor.mSpe2PepFileReader.sipPSMs;
    List psmList(sipPSMs.size());
    List isotopicPeaksList;
    std::vector<std::string> isotopicMZs, isotopicCharges, isotopicIntensities;
    isotopicMZs.reserve(sipPSMs[0].fileNames.size());
    isotopicCharges.reserve(sipPSMs[0].fileNames.size());
    isotopicIntensities.reserve(sipPSMs[0].fileNames.size());
    std::string isotopicMZstr, isotopicChargeStr, isotopicIntensityStr;
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        isotopicMZs.clear();
        isotopicCharges.clear();
        isotopicIntensities.clear();
        for (size_t j = 0; j < sipPSMs[i].isotopicPeakss.size(); j++)
        {
            isotopicMZstr.clear();
            isotopicChargeStr.clear();
            isotopicIntensityStr.clear();
            for (size_t k = 0; k < sipPSMs[i].isotopicPeakss[j].size(); k++)
            {
                isotopicMZstr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].mz) + ",";
                isotopicChargeStr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].charge) + ",";
                isotopicIntensityStr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].intensity) + ",";
            }
            // remove last ","
            isotopicMZstr.pop_back();
            isotopicChargeStr.pop_back();
            isotopicIntensityStr.pop_back();
            isotopicMZs.push_back(isotopicMZstr);
            isotopicCharges.push_back(isotopicChargeStr);
            isotopicIntensities.push_back(isotopicIntensityStr);
        }
        DataFrame psmDf = ListBuilder()
                              .add("fileNames", wrap(sipPSMs[i].fileNames))
                              .add("scanNumbers", wrap(sipPSMs[i].scanNumbers))
                              .add("precursorScanNumbers", wrap(sipPSMs[i].precursorScanNumbers))
                              .add("retentionTimes", wrap(sipPSMs[i].retentionTimes))
                              .add("isDecoys", wrap(sipPSMs[i].isDecoys))
                              .add("parentCharges", wrap(sipPSMs[i].parentCharges))
                              .add("isolationWindowCenterMZs", wrap(sipPSMs[i].isolationWindowCenterMZs))
                              .add("measuredParentMasses", wrap(sipPSMs[i].measuredParentMasses))
                              .add("calculatedParentMasses", wrap(sipPSMs[i].calculatedParentMasses))
                              .add("mzShiftFromisolationWindowCenters", wrap(sipPSMs[i].mzShiftFromisolationWindowCenters))
                              .add("isotopicMassWindowShifts", wrap(sipPSMs[i].isotopicMassWindowShifts))
                              .add("massErrors", wrap(sipPSMs[i].massErrors))
                              .add("peptideLengths", wrap(sipPSMs[i].peptideLengths))
                              .add("missCleavageSiteNumbers", wrap(sipPSMs[i].missCleavageSiteNumbers))
                              .add("PTMnumbers", wrap(sipPSMs[i].PTMnumbers))
                              .add("searchNames", wrap(sipPSMs[i].searchNames))
                              .add("WDPscores", wrap(sipPSMs[i].WDPscores))
                              .add("XcorrScores", wrap(sipPSMs[i].XcorrScores))
                              .add("MVHscores", wrap(sipPSMs[i].MVHscores))                              
                              .add("diffScores", wrap(sipPSMs[i].diffScores))
                              .add("ranks", wrap(sipPSMs[i].ranks))
                              .add("identifiedPeptides", wrap(sipPSMs[i].identifiedPeptides))
                              .add("originalPeptides", wrap(sipPSMs[i].originalPeptides))
                              .add("nakePeptides", wrap(sipPSMs[i].nakePeptides))
                              .add("proteinNames", wrap(sipPSMs[i].proteinNames))
                              .add("istopicPeakNumbers", wrap(sipPSMs[i].isotopicPeakNumbers))
                              .add("MS1IsotopicAbundances", wrap(sipPSMs[i].MS1IsotopicAbundances))
                              .add("MS2IsotopicAbundances", wrap(sipPSMs[i].MS2IsotopicAbundances))
                              .add("isotopicMZs", wrap(isotopicMZs))
                              .add("isotopicCharges", wrap(isotopicCharges))
                              .add("isotopicIntensities", wrap(isotopicIntensities))
                              .add("isotopicAbundanceDiffs", wrap(sipPSMs[i].isotopicAbundanceDiffs));
        psmList[i] = psmDf;
    }
    psmList.names() = extractor.mSpe2PepFileReader.FT2s;
    return psmList;
}

//' extractPSMfeaturesTargetAndDecoy extract featueres of top PSMs from target and decoy .Spe2Pep.txt files
//' @param targetPath a full path with target .Spe2PepFile.txt files in it
//' @param decoyPath a full path with decoy .Spe2PepFile.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @param ftFilepath a full path with .FT1 and .FT2 files in it
//' @param ThreadNumber read ThreadNumber of FT file at the same time, it will increase ram usage
//' @return the PSMs in a dataframe in a list
//' @examples
//' psm <- extractPSMfeaturesTargetAndDecoy("targetPath", "decoyPath", topN, "ftFilepath", 3)
//' @export
// [[Rcpp::export]]

List extractPSMfeaturesTargetAndDecoy(String targetPath, String decoyPath, int topN,
                                      String ftFilepath, int ThreadNumber = 3)
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    PSMfeatureExtractor extractor;
    extractor.extractPSMfeatureParallel(targetPath, decoyPath, topN, ftFilepath, ThreadNumber);
    std::vector<sipPSM> &sipPSMs = extractor.mSpe2PepFileReader.sipPSMs;
    List psmList(sipPSMs.size());
    List isotopicPeaksList;
    std::vector<std::string> isotopicMZs, isotopicCharges, isotopicIntensities;
    isotopicMZs.reserve(sipPSMs[0].fileNames.size());
    isotopicCharges.reserve(sipPSMs[0].fileNames.size());
    isotopicIntensities.reserve(sipPSMs[0].fileNames.size());
    std::string isotopicMZstr, isotopicChargeStr, isotopicIntensityStr;
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        isotopicMZs.clear();
        isotopicCharges.clear();
        isotopicIntensities.clear();
        for (size_t j = 0; j < sipPSMs[i].isotopicPeakss.size(); j++)
        {
            isotopicMZstr.clear();
            isotopicChargeStr.clear();
            isotopicIntensityStr.clear();
            for (size_t k = 0; k < sipPSMs[i].isotopicPeakss[j].size(); k++)
            {
                isotopicMZstr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].mz) + ",";
                isotopicChargeStr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].charge) + ",";
                isotopicIntensityStr += std::to_string(sipPSMs[i].isotopicPeakss[j][k].intensity) + ",";
            }
            // remove last ","
            if (!isotopicMZstr.empty())
            {
                isotopicMZstr.pop_back();
                isotopicChargeStr.pop_back();
                isotopicIntensityStr.pop_back();
            }
            isotopicMZs.push_back(isotopicMZstr);
            isotopicCharges.push_back(isotopicChargeStr);
            isotopicIntensities.push_back(isotopicIntensityStr);
        }
        DataFrame psmDf = ListBuilder()
                              .add("fileNames", wrap(sipPSMs[i].fileNames))
                              .add("scanNumbers", wrap(sipPSMs[i].scanNumbers))
                              .add("precursorScanNumbers", wrap(sipPSMs[i].precursorScanNumbers))
                              .add("retentionTimes", wrap(sipPSMs[i].retentionTimes))
                              .add("isDecoys", wrap(sipPSMs[i].isDecoys))
                              .add("parentCharges", wrap(sipPSMs[i].parentCharges))
                              .add("isolationWindowCenterMZs", wrap(sipPSMs[i].isolationWindowCenterMZs))
                              .add("measuredParentMasses", wrap(sipPSMs[i].measuredParentMasses))
                              .add("calculatedParentMasses", wrap(sipPSMs[i].calculatedParentMasses))
                              .add("mzShiftFromisolationWindowCenters", wrap(sipPSMs[i].mzShiftFromisolationWindowCenters))
                              .add("isotopicMassWindowShifts", wrap(sipPSMs[i].isotopicMassWindowShifts))
                              .add("massErrors", wrap(sipPSMs[i].massErrors))
                              .add("peptideLengths", wrap(sipPSMs[i].peptideLengths))
                              .add("missCleavageSiteNumbers", wrap(sipPSMs[i].missCleavageSiteNumbers))
                              .add("PTMnumbers", wrap(sipPSMs[i].PTMnumbers))
                              .add("searchNames", wrap(sipPSMs[i].searchNames))
                              .add("WDPscores", wrap(sipPSMs[i].WDPscores))
                              .add("XcorrScores", wrap(sipPSMs[i].XcorrScores))
                              .add("MVHscores", wrap(sipPSMs[i].MVHscores))
                              .add("diffScores", wrap(sipPSMs[i].diffScores))
                              .add("ranks", wrap(sipPSMs[i].ranks))
                              .add("identifiedPeptides", wrap(sipPSMs[i].identifiedPeptides))
                              .add("originalPeptides", wrap(sipPSMs[i].originalPeptides))
                              .add("nakePeptides", wrap(sipPSMs[i].nakePeptides))
                              .add("proteinNames", wrap(sipPSMs[i].proteinNames))
                              .add("isotopicPeakNumbers", wrap(sipPSMs[i].isotopicPeakNumbers))
                              .add("MS1IsotopicAbundances", wrap(sipPSMs[i].MS1IsotopicAbundances))
                              .add("MS2IsotopicAbundances", wrap(sipPSMs[i].MS2IsotopicAbundances))
                              .add("isotopicMZs", wrap(isotopicMZs))
                              .add("isotopicCharges", wrap(isotopicCharges))
                              .add("isotopicIntensities", wrap(isotopicIntensities))
                              .add("isotopicAbundanceDiffs", wrap(sipPSMs[i].isotopicAbundanceDiffs));
        psmList[i] = psmDf;
    }
    psmList.names() = extractor.mSpe2PepFileReader.FT2s;
    return psmList;
}

//' extractPSMfeaturesTargetAndDecoytoPercolatorPin extract featueres of top PSMs from target and decoy .Spe2Pep.txt files
//' to pecorlator pin format
//' @param targetPath a full path with target .Spe2PepFile.txt files in it
//' @param decoyPath a full path with decoy .Spe2PepFile.txt files in it
//' @param topN store top N PSMs of each scan of one .FT2 file
//' @param ftFilepath a full path with .FT1 and .FT2 files in it
//' @param ThreadNumber read ThreadNumber of FT file at the same time, it will increase ram usage
//' @param doProteinInference out put protein inference format or only PSM format
//' @param fileName output path of the percolator tsv file
//' @examples
//' extractPSMfeaturesTargetAndDecoytoPercolatorPin("targetPath", "decoyPath", topN, "ftFilepath", 3, "a.pin")
//' @export
// [[Rcpp::export]]

void extractPSMfeaturesTargetAndDecoytoPercolatorPin(String targetPath, String decoyPath, int topN,
                                                     String ftFilepath, int ThreadNumber = 3,
                                                     bool doProteinInference = false,
                                                     String fileName = "a.pin")
{
    // read default config
    string config = get_extdata();
    ProNovoConfig::setFilename(config);
    PSMfeatureExtractor extractor;
    extractor.extractPSMfeatureParallel(targetPath, decoyPath, topN, ftFilepath, ThreadNumber);
    extractor.writePecorlatorPin(fileName, doProteinInference);
}