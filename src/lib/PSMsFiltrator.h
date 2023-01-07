#pragma once
#include "sipFileReader.h"
#include <tuple>

struct alignas(64) sipPSMinfo : sipPSM
{
    vector<float> retentionTimes;
    vector<bool> isDecoys;
    // sip abundance
    vector<int> pcts;
    vector<int> pepLengths;
    vector<int> proCounts;
    vector<string> psmIDs;
    vector<string> trimedProteinNames;
    vector<string> realPepSeqs;
    vector<string> formatedPepSeqs;
    // for accurate isotopic abundance calculation and quantify
    vector<vector<double>> precursorIsotopicMasses;
    vector<vector<double>> precursorIsotopicIntensities;
};

class PSMsFiltrator
{
private:
public:
    float FDRthreshold;
    string sipPath, ftPath;
    vector<string> tokens;
    PSMsFiltrator(const string &msipPath, const string &mftPath);
    PSMsFiltrator(const string &workPath);
    ~PSMsFiltrator();
    void splitString(const string &mString);
    // <proCount,isDecoy,protein name without fake decoy>
    tuple<int, bool, string> detectProDecoy(const string &proteinName);
    float getRentionTime();
    int getPct(const string &searchName);
    // <pepLength, pepSeq in "[]", formated peptide seq for percolator>
    tuple<int, string, string> getRealPep(const string &pepSeq);
    // <precursorIsotopicMass, precursorIsotopicIntensity>
    pair<vector<double>, vector<double>> getPrecursorIsotopicPeak();
    void writePercolatorTSV();
    void readPercolatorTSV();
    sipPSMinfo convertFilesScansTopPSMs(const unordered_map<string, unordered_map<int, vector<topPSM>>> &filesScansTopPSMs);
};
