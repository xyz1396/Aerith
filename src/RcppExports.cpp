// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// extractPSMfeatures
List extractPSMfeatures(String Spe2PepFilePath, int topN, String ftFilepath, int ThreadNumber);
RcppExport SEXP _Aerith_extractPSMfeatures(SEXP Spe2PepFilePathSEXP, SEXP topNSEXP, SEXP ftFilepathSEXP, SEXP ThreadNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type Spe2PepFilePath(Spe2PepFilePathSEXP);
    Rcpp::traits::input_parameter< int >::type topN(topNSEXP);
    Rcpp::traits::input_parameter< String >::type ftFilepath(ftFilepathSEXP);
    Rcpp::traits::input_parameter< int >::type ThreadNumber(ThreadNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(extractPSMfeatures(Spe2PepFilePath, topN, ftFilepath, ThreadNumber));
    return rcpp_result_gen;
END_RCPP
}
// extractPSMfeaturesTargetAndDecoy
List extractPSMfeaturesTargetAndDecoy(String targetPath, String decoyPath, int topN, String ftFilepath, int ThreadNumber);
RcppExport SEXP _Aerith_extractPSMfeaturesTargetAndDecoy(SEXP targetPathSEXP, SEXP decoyPathSEXP, SEXP topNSEXP, SEXP ftFilepathSEXP, SEXP ThreadNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type targetPath(targetPathSEXP);
    Rcpp::traits::input_parameter< String >::type decoyPath(decoyPathSEXP);
    Rcpp::traits::input_parameter< int >::type topN(topNSEXP);
    Rcpp::traits::input_parameter< String >::type ftFilepath(ftFilepathSEXP);
    Rcpp::traits::input_parameter< int >::type ThreadNumber(ThreadNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(extractPSMfeaturesTargetAndDecoy(targetPath, decoyPath, topN, ftFilepath, ThreadNumber));
    return rcpp_result_gen;
END_RCPP
}
// extractPSMfeaturesTargetAndDecoytoPercolatorPin
void extractPSMfeaturesTargetAndDecoytoPercolatorPin(String targetPath, String decoyPath, int topN, String ftFilepath, int ThreadNumber, bool doProteinInference, String fileName);
RcppExport SEXP _Aerith_extractPSMfeaturesTargetAndDecoytoPercolatorPin(SEXP targetPathSEXP, SEXP decoyPathSEXP, SEXP topNSEXP, SEXP ftFilepathSEXP, SEXP ThreadNumberSEXP, SEXP doProteinInferenceSEXP, SEXP fileNameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type targetPath(targetPathSEXP);
    Rcpp::traits::input_parameter< String >::type decoyPath(decoyPathSEXP);
    Rcpp::traits::input_parameter< int >::type topN(topNSEXP);
    Rcpp::traits::input_parameter< String >::type ftFilepath(ftFilepathSEXP);
    Rcpp::traits::input_parameter< int >::type ThreadNumber(ThreadNumberSEXP);
    Rcpp::traits::input_parameter< bool >::type doProteinInference(doProteinInferenceSEXP);
    Rcpp::traits::input_parameter< String >::type fileName(fileNameSEXP);
    extractPSMfeaturesTargetAndDecoytoPercolatorPin(targetPath, decoyPath, topN, ftFilepath, ThreadNumber, doProteinInference, fileName);
    return R_NilValue;
END_RCPP
}
// getUnfilteredPSMs
DataFrame getUnfilteredPSMs(String sipPath, String ftPath, size_t topN);
RcppExport SEXP _Aerith_getUnfilteredPSMs(SEXP sipPathSEXP, SEXP ftPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type sipPath(sipPathSEXP);
    Rcpp::traits::input_parameter< String >::type ftPath(ftPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(getUnfilteredPSMs(sipPath, ftPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// getUnfilteredPeptides
DataFrame getUnfilteredPeptides(CharacterVector workingPath);
RcppExport SEXP _Aerith_getUnfilteredPeptides(SEXP workingPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    rcpp_result_gen = Rcpp::wrap(getUnfilteredPeptides(workingPath));
    return rcpp_result_gen;
END_RCPP
}
// getFilterThreshold
DataFrame getFilterThreshold(CharacterVector workingPath, NumericVector OverallThreshold);
RcppExport SEXP _Aerith_getFilterThreshold(SEXP workingPathSEXP, SEXP OverallThresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type OverallThreshold(OverallThresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(getFilterThreshold(workingPath, OverallThreshold));
    return rcpp_result_gen;
END_RCPP
}
// getFilterThresholdTopPSMs
List getFilterThresholdTopPSMs(CharacterVector workingPath, NumericVector OverallThreshold, size_t topN);
RcppExport SEXP _Aerith_getFilterThresholdTopPSMs(SEXP workingPathSEXP, SEXP OverallThresholdSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type OverallThreshold(OverallThresholdSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(getFilterThresholdTopPSMs(workingPath, OverallThreshold, topN));
    return rcpp_result_gen;
END_RCPP
}
// getFilterThresholdTopPSMsSpe2Pep
List getFilterThresholdTopPSMsSpe2Pep(String workingPath, float OverallThreshold, size_t topN);
RcppExport SEXP _Aerith_getFilterThresholdTopPSMsSpe2Pep(SEXP workingPathSEXP, SEXP OverallThresholdSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< float >::type OverallThreshold(OverallThresholdSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(getFilterThresholdTopPSMsSpe2Pep(workingPath, OverallThreshold, topN));
    return rcpp_result_gen;
END_RCPP
}
// generateOneCFG
bool generateOneCFG(String cfgPath, String outPath, String element, int pct, int center, int width);
RcppExport SEXP _Aerith_generateOneCFG(SEXP cfgPathSEXP, SEXP outPathSEXP, SEXP elementSEXP, SEXP pctSEXP, SEXP centerSEXP, SEXP widthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type cfgPath(cfgPathSEXP);
    Rcpp::traits::input_parameter< String >::type outPath(outPathSEXP);
    Rcpp::traits::input_parameter< String >::type element(elementSEXP);
    Rcpp::traits::input_parameter< int >::type pct(pctSEXP);
    Rcpp::traits::input_parameter< int >::type center(centerSEXP);
    Rcpp::traits::input_parameter< int >::type width(widthSEXP);
    rcpp_result_gen = Rcpp::wrap(generateOneCFG(cfgPath, outPath, element, pct, center, width));
    return rcpp_result_gen;
END_RCPP
}
// generateCFGs
bool generateCFGs(String cfgPath, String outPath, String element);
RcppExport SEXP _Aerith_generateCFGs(SEXP cfgPathSEXP, SEXP outPathSEXP, SEXP elementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type cfgPath(cfgPathSEXP);
    Rcpp::traits::input_parameter< String >::type outPath(outPathSEXP);
    Rcpp::traits::input_parameter< String >::type element(elementSEXP);
    rcpp_result_gen = Rcpp::wrap(generateCFGs(cfgPath, outPath, element));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator
DataFrame precursor_peak_calculator(CharacterVector AAstr);
RcppExport SEXP _Aerith_precursor_peak_calculator(SEXP AAstrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator(AAstr));
    return rcpp_result_gen;
END_RCPP
}
// residue_peak_calculator_DIY
DataFrame residue_peak_calculator_DIY(String residue, String Atom, double Prob);
RcppExport SEXP _Aerith_residue_peak_calculator_DIY(SEXP residueSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type residue(residueSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(residue_peak_calculator_DIY(residue, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator_DIY
DataFrame precursor_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom, NumericVector Prob);
RcppExport SEXP _Aerith_precursor_peak_calculator_DIY(SEXP AAstrSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator_DIY(AAstr, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// calPepAtomCount
DataFrame calPepAtomCount(StringVector AAstrs);
RcppExport SEXP _Aerith_calPepAtomCount(SEXP AAstrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    rcpp_result_gen = Rcpp::wrap(calPepAtomCount(AAstrs));
    return rcpp_result_gen;
END_RCPP
}
// calBYAtomCountAndBaseMass
List calBYAtomCountAndBaseMass(StringVector AAstrs);
RcppExport SEXP _Aerith_calBYAtomCountAndBaseMass(SEXP AAstrsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    rcpp_result_gen = Rcpp::wrap(calBYAtomCountAndBaseMass(AAstrs));
    return rcpp_result_gen;
END_RCPP
}
// calPepPrecursorMass
NumericVector calPepPrecursorMass(StringVector AAstrs, String Atom, NumericVector Probs);
RcppExport SEXP _Aerith_calPepPrecursorMass(SEXP AAstrsSEXP, SEXP AtomSEXP, SEXP ProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Probs(ProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(calPepPrecursorMass(AAstrs, Atom, Probs));
    return rcpp_result_gen;
END_RCPP
}
// calPepNeutronMass
NumericVector calPepNeutronMass(StringVector AAstrs, String Atom, NumericVector Probs);
RcppExport SEXP _Aerith_calPepNeutronMass(SEXP AAstrsSEXP, SEXP AtomSEXP, SEXP ProbsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Probs(ProbsSEXP);
    rcpp_result_gen = Rcpp::wrap(calPepNeutronMass(AAstrs, Atom, Probs));
    return rcpp_result_gen;
END_RCPP
}
// precursor_peak_calculator_DIY_averagine
List precursor_peak_calculator_DIY_averagine(StringVector AAstrs, String Atom, double Prob);
RcppExport SEXP _Aerith_precursor_peak_calculator_DIY_averagine(SEXP AAstrsSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< StringVector >::type AAstrs(AAstrsSEXP);
    Rcpp::traits::input_parameter< String >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(precursor_peak_calculator_DIY_averagine(AAstrs, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// BYion_peak_calculator_DIY
DataFrame BYion_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom, NumericVector Prob);
RcppExport SEXP _Aerith_BYion_peak_calculator_DIY(SEXP AAstrSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type AAstr(AAstrSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(BYion_peak_calculator_DIY(AAstr, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// readOneScanMS2
List readOneScanMS2(const String& ftFile, const size_t scanNumber);
RcppExport SEXP _Aerith_readOneScanMS2(SEXP ftFileSEXP, SEXP scanNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String& >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type scanNumber(scanNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(readOneScanMS2(ftFile, scanNumber));
    return rcpp_result_gen;
END_RCPP
}
// readOneScanMS1
List readOneScanMS1(const String& ftFile, const size_t scanNumber);
RcppExport SEXP _Aerith_readOneScanMS1(SEXP ftFileSEXP, SEXP scanNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String& >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type scanNumber(scanNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(readOneScanMS1(ftFile, scanNumber));
    return rcpp_result_gen;
END_RCPP
}
// readFTheader
List readFTheader(String ftFile);
RcppExport SEXP _Aerith_readFTheader(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readFTheader(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS1
List readScansMS1(const String ftFile, const size_t startScanNumber, const size_t endScanNumber);
RcppExport SEXP _Aerith_readScansMS1(SEXP ftFileSEXP, SEXP startScanNumberSEXP, SEXP endScanNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type startScanNumber(startScanNumberSEXP);
    Rcpp::traits::input_parameter< const size_t >::type endScanNumber(endScanNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS1(ftFile, startScanNumber, endScanNumber));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS1Vector
List readScansMS1Vector(const String ftFile, const NumericVector scanNumbersVector);
RcppExport SEXP _Aerith_readScansMS1Vector(SEXP ftFileSEXP, SEXP scanNumbersVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type scanNumbersVector(scanNumbersVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS1Vector(ftFile, scanNumbersVector));
    return rcpp_result_gen;
END_RCPP
}
// readAllScanMS1
List readAllScanMS1(const String ftFile);
RcppExport SEXP _Aerith_readAllScanMS1(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readAllScanMS1(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS2
List readScansMS2(const String ftFile, const size_t startScanNumber, const size_t endScanNumber);
RcppExport SEXP _Aerith_readScansMS2(SEXP ftFileSEXP, SEXP startScanNumberSEXP, SEXP endScanNumberSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type startScanNumber(startScanNumberSEXP);
    Rcpp::traits::input_parameter< const size_t >::type endScanNumber(endScanNumberSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS2(ftFile, startScanNumber, endScanNumber));
    return rcpp_result_gen;
END_RCPP
}
// readScansMS2Vector
List readScansMS2Vector(const String ftFile, const NumericVector scanNumbersVector);
RcppExport SEXP _Aerith_readScansMS2Vector(SEXP ftFileSEXP, SEXP scanNumbersVectorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type scanNumbersVector(scanNumbersVectorSEXP);
    rcpp_result_gen = Rcpp::wrap(readScansMS2Vector(ftFile, scanNumbersVector));
    return rcpp_result_gen;
END_RCPP
}
// readAllScanMS2
List readAllScanMS2(const String ftFile);
RcppExport SEXP _Aerith_readAllScanMS2(SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const String >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readAllScanMS2(ftFile));
    return rcpp_result_gen;
END_RCPP
}
// readSip
List readSip(CharacterVector sipFile);
RcppExport SEXP _Aerith_readSip(SEXP sipFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type sipFile(sipFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readSip(sipFile));
    return rcpp_result_gen;
END_RCPP
}
// readSips
List readSips(CharacterVector workingPath);
RcppExport SEXP _Aerith_readSips(SEXP workingPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    rcpp_result_gen = Rcpp::wrap(readSips(workingPath));
    return rcpp_result_gen;
END_RCPP
}
// readFilesScansTopPSMs
DataFrame readFilesScansTopPSMs(CharacterVector workingPath, size_t topN);
RcppExport SEXP _Aerith_readFilesScansTopPSMs(SEXP workingPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< CharacterVector >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readFilesScansTopPSMs(workingPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// readFilesScansTopPSMsFromOneFT2
DataFrame readFilesScansTopPSMsFromOneFT2(String workingPath, String pattern, size_t topN);
RcppExport SEXP _Aerith_readFilesScansTopPSMsFromOneFT2(SEXP workingPathSEXP, SEXP patternSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< String >::type pattern(patternSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readFilesScansTopPSMsFromOneFT2(workingPath, pattern, topN));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2Pep
List readSpe2Pep(String Spe2PepFile);
RcppExport SEXP _Aerith_readSpe2Pep(SEXP Spe2PepFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type Spe2PepFile(Spe2PepFileSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2Pep(Spe2PepFile));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2Peps
List readSpe2Peps(String workingPath);
RcppExport SEXP _Aerith_readSpe2Peps(SEXP workingPathSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2Peps(workingPath));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2PepFilesScansTopPSMs
DataFrame readSpe2PepFilesScansTopPSMs(String workingPath, size_t topN);
RcppExport SEXP _Aerith_readSpe2PepFilesScansTopPSMs(SEXP workingPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2PepFilesScansTopPSMs(workingPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2PepFilesScansTopPSMsFromOneFT2
DataFrame readSpe2PepFilesScansTopPSMsFromOneFT2(String workingPath, String pattern, size_t topN);
RcppExport SEXP _Aerith_readSpe2PepFilesScansTopPSMsFromOneFT2(SEXP workingPathSEXP, SEXP patternSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< String >::type pattern(patternSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2PepFilesScansTopPSMsFromOneFT2(workingPath, pattern, topN));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2PepFilesScansTopPSMsFromEachFT2Parallel
List readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(String workingPath, size_t topN);
RcppExport SEXP _Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(SEXP workingPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(workingPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel
List readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(String targetPath, String decoyPath, size_t topN);
RcppExport SEXP _Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(SEXP targetPathSEXP, SEXP decoyPathSEXP, SEXP topNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type targetPath(targetPathSEXP);
    Rcpp::traits::input_parameter< String >::type decoyPath(decoyPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    rcpp_result_gen = Rcpp::wrap(readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(targetPath, decoyPath, topN));
    return rcpp_result_gen;
END_RCPP
}
// writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel
void writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(String workingPath, size_t topN, String fileName);
RcppExport SEXP _Aerith_writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(SEXP workingPathSEXP, SEXP topNSEXP, SEXP fileNameSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< String >::type workingPath(workingPathSEXP);
    Rcpp::traits::input_parameter< size_t >::type topN(topNSEXP);
    Rcpp::traits::input_parameter< String >::type fileName(fileNameSEXP);
    writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel(workingPath, topN, fileName);
    return R_NilValue;
END_RCPP
}
// scoreIntensity
double scoreIntensity(const bool observed, const double realIntensity, const double expectedIntensity, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scoreIntensity(SEXP observedSEXP, SEXP realIntensitySEXP, SEXP expectedIntensitySEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const bool >::type observed(observedSEXP);
    Rcpp::traits::input_parameter< const double >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const double >::type expectedIntensity(expectedIntensitySEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scoreIntensity(observed, realIntensity, expectedIntensity, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// scoreIntensityByCE
double scoreIntensityByCE(const NumericVector& expectedIntensity, const NumericVector& observedIntensity);
RcppExport SEXP _Aerith_scoreIntensityByCE(SEXP expectedIntensitySEXP, SEXP observedIntensitySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type expectedIntensity(expectedIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type observedIntensity(observedIntensitySEXP);
    rcpp_result_gen = Rcpp::wrap(scoreIntensityByCE(expectedIntensity, observedIntensity));
    return rcpp_result_gen;
END_RCPP
}
// scorePSM
double scorePSM(const NumericVector& realMZ, const NumericVector& realIntensity, const NumericVector& realCharge, int parentCharge, const String& pepSeq, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scorePSM(SEXP realMZSEXP, SEXP realIntensitySEXP, SEXP realChargeSEXP, SEXP parentChargeSEXP, SEXP pepSeqSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type realMZ(realMZSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realCharge(realChargeSEXP);
    Rcpp::traits::input_parameter< int >::type parentCharge(parentChargeSEXP);
    Rcpp::traits::input_parameter< const String& >::type pepSeq(pepSeqSEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scorePSM(realMZ, realIntensity, realCharge, parentCharge, pepSeq, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// annotatePSM
List annotatePSM(const NumericVector& realMZ, const NumericVector& realIntensity, const NumericVector& realCharge, const String& pepSeq, const NumericVector charges, const String& Atom, double Prob, const double isoCenter, const double isoWidth);
RcppExport SEXP _Aerith_annotatePSM(SEXP realMZSEXP, SEXP realIntensitySEXP, SEXP realChargeSEXP, SEXP pepSeqSEXP, SEXP chargesSEXP, SEXP AtomSEXP, SEXP ProbSEXP, SEXP isoCenterSEXP, SEXP isoWidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type realMZ(realMZSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realCharge(realChargeSEXP);
    Rcpp::traits::input_parameter< const String& >::type pepSeq(pepSeqSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type charges(chargesSEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    Rcpp::traits::input_parameter< const double >::type isoCenter(isoCenterSEXP);
    Rcpp::traits::input_parameter< const double >::type isoWidth(isoWidthSEXP);
    rcpp_result_gen = Rcpp::wrap(annotatePSM(realMZ, realIntensity, realCharge, pepSeq, charges, Atom, Prob, isoCenter, isoWidth));
    return rcpp_result_gen;
END_RCPP
}
// scorePSMold
double scorePSMold(const NumericVector& realMZ, const NumericVector& realIntensity, const NumericVector& realCharge, const String& pepSeq, const String& Atom, double Prob);
RcppExport SEXP _Aerith_scorePSMold(SEXP realMZSEXP, SEXP realIntensitySEXP, SEXP realChargeSEXP, SEXP pepSeqSEXP, SEXP AtomSEXP, SEXP ProbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type realMZ(realMZSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realIntensity(realIntensitySEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type realCharge(realChargeSEXP);
    Rcpp::traits::input_parameter< const String& >::type pepSeq(pepSeqSEXP);
    Rcpp::traits::input_parameter< const String& >::type Atom(AtomSEXP);
    Rcpp::traits::input_parameter< double >::type Prob(ProbSEXP);
    rcpp_result_gen = Rcpp::wrap(scorePSMold(realMZ, realIntensity, realCharge, pepSeq, Atom, Prob));
    return rcpp_result_gen;
END_RCPP
}
// rankyfify
NumericVector rankyfify(NumericVector a);
RcppExport SEXP _Aerith_rankyfify(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(rankyfify(a));
    return rcpp_result_gen;
END_RCPP
}
// denoiseOneMS2ScanHasCharge
List denoiseOneMS2ScanHasCharge(List scanList, float window, float step, float threshold);
RcppExport SEXP _Aerith_denoiseOneMS2ScanHasCharge(SEXP scanListSEXP, SEXP windowSEXP, SEXP stepSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type scanList(scanListSEXP);
    Rcpp::traits::input_parameter< float >::type window(windowSEXP);
    Rcpp::traits::input_parameter< float >::type step(stepSEXP);
    Rcpp::traits::input_parameter< float >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(denoiseOneMS2ScanHasCharge(scanList, window, step, threshold));
    return rcpp_result_gen;
END_RCPP
}
// writeAllScanMS2
bool writeAllScanMS2(List header, List scansList, const String& ftFile);
RcppExport SEXP _Aerith_writeAllScanMS2(SEXP headerSEXP, SEXP scansListSEXP, SEXP ftFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type header(headerSEXP);
    Rcpp::traits::input_parameter< List >::type scansList(scansListSEXP);
    Rcpp::traits::input_parameter< const String& >::type ftFile(ftFileSEXP);
    rcpp_result_gen = Rcpp::wrap(writeAllScanMS2(header, scansList, ftFile));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Aerith_extractPSMfeatures", (DL_FUNC) &_Aerith_extractPSMfeatures, 4},
    {"_Aerith_extractPSMfeaturesTargetAndDecoy", (DL_FUNC) &_Aerith_extractPSMfeaturesTargetAndDecoy, 5},
    {"_Aerith_extractPSMfeaturesTargetAndDecoytoPercolatorPin", (DL_FUNC) &_Aerith_extractPSMfeaturesTargetAndDecoytoPercolatorPin, 7},
    {"_Aerith_getUnfilteredPSMs", (DL_FUNC) &_Aerith_getUnfilteredPSMs, 3},
    {"_Aerith_getUnfilteredPeptides", (DL_FUNC) &_Aerith_getUnfilteredPeptides, 1},
    {"_Aerith_getFilterThreshold", (DL_FUNC) &_Aerith_getFilterThreshold, 2},
    {"_Aerith_getFilterThresholdTopPSMs", (DL_FUNC) &_Aerith_getFilterThresholdTopPSMs, 3},
    {"_Aerith_getFilterThresholdTopPSMsSpe2Pep", (DL_FUNC) &_Aerith_getFilterThresholdTopPSMsSpe2Pep, 3},
    {"_Aerith_generateOneCFG", (DL_FUNC) &_Aerith_generateOneCFG, 6},
    {"_Aerith_generateCFGs", (DL_FUNC) &_Aerith_generateCFGs, 3},
    {"_Aerith_precursor_peak_calculator", (DL_FUNC) &_Aerith_precursor_peak_calculator, 1},
    {"_Aerith_residue_peak_calculator_DIY", (DL_FUNC) &_Aerith_residue_peak_calculator_DIY, 3},
    {"_Aerith_precursor_peak_calculator_DIY", (DL_FUNC) &_Aerith_precursor_peak_calculator_DIY, 3},
    {"_Aerith_calPepAtomCount", (DL_FUNC) &_Aerith_calPepAtomCount, 1},
    {"_Aerith_calBYAtomCountAndBaseMass", (DL_FUNC) &_Aerith_calBYAtomCountAndBaseMass, 1},
    {"_Aerith_calPepPrecursorMass", (DL_FUNC) &_Aerith_calPepPrecursorMass, 3},
    {"_Aerith_calPepNeutronMass", (DL_FUNC) &_Aerith_calPepNeutronMass, 3},
    {"_Aerith_precursor_peak_calculator_DIY_averagine", (DL_FUNC) &_Aerith_precursor_peak_calculator_DIY_averagine, 3},
    {"_Aerith_BYion_peak_calculator_DIY", (DL_FUNC) &_Aerith_BYion_peak_calculator_DIY, 3},
    {"_Aerith_readOneScanMS2", (DL_FUNC) &_Aerith_readOneScanMS2, 2},
    {"_Aerith_readOneScanMS1", (DL_FUNC) &_Aerith_readOneScanMS1, 2},
    {"_Aerith_readFTheader", (DL_FUNC) &_Aerith_readFTheader, 1},
    {"_Aerith_readScansMS1", (DL_FUNC) &_Aerith_readScansMS1, 3},
    {"_Aerith_readScansMS1Vector", (DL_FUNC) &_Aerith_readScansMS1Vector, 2},
    {"_Aerith_readAllScanMS1", (DL_FUNC) &_Aerith_readAllScanMS1, 1},
    {"_Aerith_readScansMS2", (DL_FUNC) &_Aerith_readScansMS2, 3},
    {"_Aerith_readScansMS2Vector", (DL_FUNC) &_Aerith_readScansMS2Vector, 2},
    {"_Aerith_readAllScanMS2", (DL_FUNC) &_Aerith_readAllScanMS2, 1},
    {"_Aerith_readSip", (DL_FUNC) &_Aerith_readSip, 1},
    {"_Aerith_readSips", (DL_FUNC) &_Aerith_readSips, 1},
    {"_Aerith_readFilesScansTopPSMs", (DL_FUNC) &_Aerith_readFilesScansTopPSMs, 2},
    {"_Aerith_readFilesScansTopPSMsFromOneFT2", (DL_FUNC) &_Aerith_readFilesScansTopPSMsFromOneFT2, 3},
    {"_Aerith_readSpe2Pep", (DL_FUNC) &_Aerith_readSpe2Pep, 1},
    {"_Aerith_readSpe2Peps", (DL_FUNC) &_Aerith_readSpe2Peps, 1},
    {"_Aerith_readSpe2PepFilesScansTopPSMs", (DL_FUNC) &_Aerith_readSpe2PepFilesScansTopPSMs, 2},
    {"_Aerith_readSpe2PepFilesScansTopPSMsFromOneFT2", (DL_FUNC) &_Aerith_readSpe2PepFilesScansTopPSMsFromOneFT2, 3},
    {"_Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2Parallel", (DL_FUNC) &_Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2Parallel, 2},
    {"_Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel", (DL_FUNC) &_Aerith_readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel, 3},
    {"_Aerith_writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel", (DL_FUNC) &_Aerith_writeSpe2PepFilesScansTopPSMsFromEachFT2Parallel, 3},
    {"_Aerith_scoreIntensity", (DL_FUNC) &_Aerith_scoreIntensity, 5},
    {"_Aerith_scoreIntensityByCE", (DL_FUNC) &_Aerith_scoreIntensityByCE, 2},
    {"_Aerith_scorePSM", (DL_FUNC) &_Aerith_scorePSM, 7},
    {"_Aerith_annotatePSM", (DL_FUNC) &_Aerith_annotatePSM, 9},
    {"_Aerith_scorePSMold", (DL_FUNC) &_Aerith_scorePSMold, 6},
    {"_Aerith_rankyfify", (DL_FUNC) &_Aerith_rankyfify, 1},
    {"_Aerith_denoiseOneMS2ScanHasCharge", (DL_FUNC) &_Aerith_denoiseOneMS2ScanHasCharge, 4},
    {"_Aerith_writeAllScanMS2", (DL_FUNC) &_Aerith_writeAllScanMS2, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_Aerith(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
