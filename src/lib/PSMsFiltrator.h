#pragma once
#include "sipPSM.h"
#include <tuple>
#include <unordered_map>

class alignas(64) sipPSMinfo : public sipPSM
{
public:
    std::vector<float> retentionTimes;
    std::vector<bool> isDecoys;
    // sip abundance
    std::vector<int> pcts;
    std::vector<int> pepLengths;
    std::vector<int> proCounts;
    std::vector<std::string> psmIDs;
    std::vector<std::string> trimedProteinNames;
    std::vector<std::string> realPepSeqs;
    std::vector<std::string> formatedPepSeqs;
    // for accurate isotopic abundance calculation and quantify
    std::vector<std::vector<double>> precursorIsotopicMasses;
    std::vector<std::vector<double>> precursorIsotopicIntensities;
};

class PSMsFiltrator
{
private:
public:
    float FDRthreshold;
    std::string sipPath, ftPath;
    std::vector<std::string> tokens;
    PSMsFiltrator(const std::string &msipPath, const std::string &mftPath);
    PSMsFiltrator(const std::string &workPath);
    ~PSMsFiltrator();
    void splitString(const std::string &mString);
    // <proCount,isDecoy,protein name without fake decoy>
    std::tuple<int, bool, std::string> detectProDecoy(const std::string &proteinName);
    float getRentionTime();
    int getPct(const std::string &searchName);
    // <pepLength, pepSeq in "[]", formated peptide seq for percolator>
    std::tuple<int, std::string, std::string> getRealPep(const std::string &pepSeq);
    // <precursorIsotopicMass, precursorIsotopicIntensity>
    std::pair<std::vector<double>, std::vector<double>> getPrecursorIsotopicPeak();
    void writePercolatorTSV();
    void readPercolatorTSV();
    sipPSMinfo convertFilesScansTopPSMs(const std::unordered_map<std::string, std::unordered_map<int, std::vector<scanTopPSM>>> &filesScansTopPSMs);
};
