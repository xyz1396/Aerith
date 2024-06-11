#include "sipPSM.h"

scanTopPSM::scanTopPSM(int parentCharge,
                       int precursorScanNumber,
                       double isolationWindowCenterMZ,
                       double measuredParentMass,
                       double calculatedParentMass,
                       std::string searchName,
                       std::string identifiedPeptide,
                       std::string originalPeptide,
                       std::string nakePeptide,
                       std::string proteinName,
                       float retentionTime,
                       float MVHscore,
                       float XcorrScore,
                       float WDPscore,
                       bool isDecoy) : parentCharge(parentCharge),
                                       precursorScanNumber(precursorScanNumber),
                                       isolationWindowCenterMZ(isolationWindowCenterMZ),
                                       measuredParentMass(measuredParentMass),
                                       calculatedParentMass(calculatedParentMass), searchName(searchName),
                                       identifiedPeptide(identifiedPeptide), originalPeptide(originalPeptide),
                                       nakePeptide(nakePeptide),
                                       proteinName(proteinName), retentionTime(retentionTime),
                                       MVHscore(MVHscore), XcorrScore(XcorrScore), WDPscore(WDPscore),isDecoy(isDecoy)
{
}

scanTopPSM::scanTopPSM(int parentCharge,
                       double measuredParentMass,
                       double calculatedParentMass,
                       std::string searchName,
                       float score,
                       std::string identifiedPeptide,
                       std::string originalPeptide,
                       std::string proteinName) : parentCharge(parentCharge), measuredParentMass(measuredParentMass),
                                                  calculatedParentMass(calculatedParentMass), searchName(searchName),
                                                  score(score),
                                                  identifiedPeptide(identifiedPeptide), originalPeptide(originalPeptide),
                                                  proteinName(proteinName)
{
}

scanTopPSM::~scanTopPSM() {}