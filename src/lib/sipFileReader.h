#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include "sipPSM.h"
namespace fs = std::filesystem;

class sipFileReader
{
private:
public:
    std::string workingPath;
    std::vector<std::string> sipFileNames;
    std::vector<sipPSM> sipPSMs;
    // store top N PSMs of each scan of one .FT file
    size_t topN = 5;
    // unordered_map<.FT file name, unordered_map<scanNumber, std::vector<top N PSM>>>
    // store topN PSMs of each scan of each .FT file
    std::unordered_map<std::string, std::unordered_map<int, std::vector<scanTopPSM>>> filesScansTopPSMs;
    std::vector<std::string> tokens;
    sipPSM currentSipPSM;
    std::string currentFilePath;
    std::fstream sipFileStream;
    // for single file
    sipFileReader();
    // for multi files;
    sipFileReader(std::string mWorkingPath);
    ~sipFileReader();
    void splitString(const std::string &mString);
    // fill vectors in currentSipPSM
    void fillVectors();
    void readOneFile(std::string sipFileName);
    void readAllFiles();
    void readAllFilesTopPSMs();
    void fillScanTopPSMs(std::vector<scanTopPSM> &scanTopPSMs, const int psmIX);
    // converted from filesScansTopPSMs
    // make it good for output
    sipPSM convertFilesScansTopPSMs();
};
