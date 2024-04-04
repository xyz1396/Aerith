#pragma once
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
    std::vector<std::string> proteinNames;

    // for .Spe2Pep.txt file
    std::vector<float> retentionTimes;
    std::vector<float> MVHscores;
    std::vector<float> XcorrScores;
    std::vector<float> WDPscores;
    std::vector<double> isolationWindowCenterMZs;

    // new featuers
    std::vector<double> massErrors;
    std::vector<double> mzShiftFromisolationWindowCenters;
    std::vector<int> isotopicMassWindowShifts;
};

class alignas(64) scanTopPSM
{
public:
    int parentCharge;
    double isolationWindowCenterMZ;
    double measuredParentMass;
    double calculatedParentMass;
    std::string searchName;
    float score;
    std::string identifiedPeptide;
    std::string originalPeptide;
    std::string proteinName;
    // for .Spe2Pep.txt file
    float retentionTime;
    float MVHscore;
    float XcorrScore;
    float WDPscore;
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
               double isolationWindowCenterMZ,
               double measuredParentMass,
               double calculatedParentMass,
               std::string searchName,
               std::string identifiedPeptide,
               std::string originalPeptide,
               std::string proteinName,
               float retentionTime,
               float MVHscore,
               float XcorrScore,
               float WDPscore);
    ~scanTopPSM();
};