#pragma once
#include "isotopicPeak.h"
#include <string>
#include <vector>

class alignas(64) sipPSM
{
public:
    std::string fileName;

    // for multiple .FT files
    std::vector<std::string> fileNames;

    std::vector<int> scanNumbers;
    std::vector<int> parentCharges;
    std::vector<double> measuredParentMasses;
    std::vector<double> calculatedParentMasses;
    std::string scanType;
    std::string searchName;

    // for multiple PCT in SIP search
    std::vector<std::string> searchNames;

    std::string scoringFunction;
    std::vector<int> ranks;
    std::vector<float> scores;
    std::vector<std::string> identifiedPeptides;
    std::vector<std::string> originalPeptides;
    std::vector<std::string> nakePeptides;
    std::vector<std::string> proteinNames;

    // for .Spe2Pep.txt file
    std::vector<int> precursorScanNumbers;
    std::vector<double> isolationWindowCenterMZs;
    std::vector<float> retentionTimes;
    std::vector<float> MVHscores;
    std::vector<float> XcorrScores;
    std::vector<float> WDPscores;

    // new featuers
    std::vector<bool> isDecoys;
    std::vector<int> peptideLengths;
    std::vector<int> missCleavageSiteNumbers;
    std::vector<int> PTMnumbers;
    std::vector<float> MVHdiffScores;
    std::vector<double> mzShiftFromisolationWindowCenters;
    std::vector<int> isotopicMassWindowShifts;
    std::vector<double> massErrors;
    // isotopic envelope in MS1
    std::vector<std::vector<isotopicPeak>> isotopicPeakss;
    std::vector<int> isotopicPeakNumbers;
    std::vector<double> precursorIntensities;
    std::vector<double> MS1IsotopicAbundances;
    std::vector<double> MS2IsotopicAbundances;
    std::vector<double> isotopicAbundanceDiffs;
};

class alignas(64) scanTopPSM
{
public:
    int parentCharge;
    int precursorScanNumber; // for .Spe2Pep.txt file
    double isolationWindowCenterMZ;
    double measuredParentMass;
    double calculatedParentMass;
    std::string searchName;
    float score;
    std::string identifiedPeptide;
    std::string originalPeptide;
    std::string nakePeptide;
    std::string proteinName;
    // for .Spe2Pep.txt file
    float retentionTime;
    float MVHscore;
    float XcorrScore;
    float WDPscore;
    bool isDecoy;
    // for .sip file
    scanTopPSM(int parentCharge,
               double measuredParentMass,
               double calculatedParentMass,
               std::string searchName,
               float score,
               std::string identifiedPeptide,
               std::string originalPeptide,
               std::string proteinName);
    // for Spe2Pep.txt file
    scanTopPSM(int parentCharge,
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
               bool isDecoy);
    ~scanTopPSM();
};