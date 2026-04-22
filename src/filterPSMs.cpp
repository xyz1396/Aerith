#include "lib/PSMsFiltrator.h"
#include "lib/sipFileReader.h"
#include <Rcpp.h>
#include <utility>
using namespace Rcpp;

//' getUnfilteredPSMs
//' @param sipPath a full path with .sip files in it
//' @param ftPath a full path with .ft files in it
//' @param topN store top N PSMs of each scan of one .FT file
//' @return data.frame of PSMs
//' @examples
//' demo_dir <- system.file("extdata", package = "Aerith")
//' head(getUnfilteredPSMs(demo_dir, demo_dir, 10))
//' @export
// [[Rcpp::export]]
DataFrame getUnfilteredPSMs(String sipPath, String ftPath, size_t topN)
{
    sipFileReader reader(sipPath);
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    PSMsFiltrator filtrator(sipPath, ftPath);
    sipPSMinfo msipPSMinfo = filtrator.convertFilesScansTopPSMs(reader.filesScansTopPSMs);
    return DataFrame::create(Named("psmIDs") = std::move(msipPSMinfo.psmIDs),
                             _("ftFileNames") = std::move(msipPSMinfo.fileNames),
                             _["scanNumbers"] = std::move(msipPSMinfo.scanNumbers),
                             _["retentionTimes"] = std::move(msipPSMinfo.retentionTimes),
                             _["scores"] = std::move(msipPSMinfo.scores),
                             _["ranks"] = std::move(msipPSMinfo.ranks),
                             _["parentCharges"] = std::move(msipPSMinfo.parentCharges),
                             _["pcts"] = std::move(msipPSMinfo.pcts),
                             _["searchNames"] = std::move(msipPSMinfo.searchNames),
                             _["isDecoys"] = std::move(msipPSMinfo.isDecoys),
                             _["measuredParentMasses"] = std::move(msipPSMinfo.measuredParentMasses),
                             _["calculatedParentMasses"] = std::move(msipPSMinfo.calculatedParentMasses),
                             _["identifiedPepSeqs"] = std::move(msipPSMinfo.identifiedPeptides),
                             _["originalPepSeqs"] = std::move(msipPSMinfo.originalPeptides),
                             _["realPepSeqs"] = std::move(msipPSMinfo.realPepSeqs),
                             _["formatedPepSeqs"] = std::move(msipPSMinfo.formatedPepSeqs),
                             _["pepLengths"] = std::move(msipPSMinfo.pepLengths),
                             _["proNames"] = std::move(msipPSMinfo.proteinNames),
                             _["trimedProteinNames"] = std::move(msipPSMinfo.trimedProteinNames),
                             _["proCounts"] = std::move(msipPSMinfo.proCounts));
}
