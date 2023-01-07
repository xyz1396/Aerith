#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>
namespace fs = std::filesystem;
using namespace std;

struct alignas(64) sipPSM
{
    string fileName;
    // for multiple .FT files
    vector<string> fileNames;
    vector<int> scanNumbers;
    vector<int> parentCharges;
    vector<double> measuredParentMasses;
    vector<double> calculatedParentMasses;
    string scanType;
    string searchName;
    // for multiple PCT in SIP search
    vector<string> searchNames;
    string scoringFunction;
    vector<int> ranks;
    vector<float> scores;
    vector<string> identifiedPeptides;
    vector<string> originalPeptides;
    vector<string> proteinNames;
};

struct alignas(64) topPSM
{
    int parentCharge;
    double measuredParentMass;
    double calculatedParentMass;
    string searchName;
    float score;
    string identifiedPeptide;
    string originalPeptide;
    string proteinName;
};

class sipFileReader
{
private:
public:
    string workingPath;
    vector<string> sipFileNames;
    vector<sipPSM> sipPSMs;
    // store top N PSMs of each scan of one .FT file
    size_t topN = 5;
    // unordered_map<.FT file name, unordered_map<scanNumber, vector<top N PSM>>>
    // store topN PSMs of each scan of each .FT file
    unordered_map<string, unordered_map<int, vector<topPSM>>> filesScansTopPSMs;
    vector<string> tokens;
    sipPSM currentSipPSM;
    string currentFilePath;
    fstream sipFileStream;
    // for single file
    sipFileReader();
    // for multi files;
    sipFileReader(string mWorkingPath);
    ~sipFileReader();
    void splitString(const string &mString);
    // fill vectors in currentSipPSM
    void fillVectors();
    void readOneFile(string sipFileName);
    void readAllFiles();
    void readAllFilesTopPSMs();
    void fillScanTopPSMs(vector<topPSM> &topPSMs, const int psmIX);
    // converted from filesScansTopPSMs
    // make it good for output
    sipPSM convertFilesScansTopPSMs();
};
