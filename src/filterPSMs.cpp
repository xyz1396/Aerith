#include "lib/PSMsFiltrator.h"
#include <Rcpp.h>
using namespace Rcpp;

//' getUnfilteredPSMs
//' @param sipPath a full path with .sip files in it
//' @param ftPath a full path with .ft files in it
//' @param topN store top N PSMs of each scan of one .FT file
//' @return a dataframe of unique PSMs and whether it is decoy sequence
//' @export
// [[Rcpp::export]]
DataFrame getUnfilteredPSMs(String sipPath, String ftPath, size_t topN)
{
    sipFileReader reader(sipPath);
    reader.topN = topN;
    reader.readAllFilesTopPSMs();
    PSMsFiltrator filtrator(sipPath, ftPath);
    sipPSMinfo msipPSMinfo = filtrator.convertFilesScansTopPSMs(reader.filesScansTopPSMs);
    return DataFrame::create(Named("psmIDs") = move(msipPSMinfo.psmIDs),
                             _("ftFileNames") = move(msipPSMinfo.fileNames),
                             _["scanNumbers"] = move(msipPSMinfo.scanNumbers),
                             _["retentionTimes"] = move(msipPSMinfo.retentionTimes),
                             _["scores"] = move(msipPSMinfo.scores),
                             _["ranks"] = move(msipPSMinfo.ranks),
                             _["parentCharges"] = move(msipPSMinfo.parentCharges),
                             _["pcts"] = move(msipPSMinfo.pcts),
                             _["searchNames"] = move(msipPSMinfo.searchNames),
                             _["isDecoys"] = move(msipPSMinfo.isDecoys),
                             _["measuredParentMasses"] = move(msipPSMinfo.measuredParentMasses),
                             _["calculatedParentMasses"] = move(msipPSMinfo.calculatedParentMasses),
                             _["identifiedPepSeqs"] = move(msipPSMinfo.identifiedPeptides),
                             _["originalPepSeqs"] = move(msipPSMinfo.originalPeptides),
                             _["realPepSeqs"] = move(msipPSMinfo.realPepSeqs),
                             _["formatedPepSeqs"] = move(msipPSMinfo.formatedPepSeqs),
                             _["pepLengths"] = move(msipPSMinfo.pepLengths),
                             _["proNames"] = move(msipPSMinfo.proteinNames),
                             _["trimedProteinNames"] = move(msipPSMinfo.trimedProteinNames),
                             _["proCounts"] = move(msipPSMinfo.proCounts));
}