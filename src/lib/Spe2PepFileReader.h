#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <omp.h>
#include <algorithm>
#include "sipPSM.h"
namespace fs = std::filesystem;

class Spe2PepFileReader
{
private:
public:
    std::string workingPath;
    std::vector<std::string> sipFileNames;
    std::vector<std::string> FT2s;
    std::vector<sipPSM> sipPSMs;
    // store top N PSMs of each scan of one .FT file
    size_t topN = 5;
    // std::unordered_map<.FT file name, std::unordered_map<scanNumber, std::vector<top N PSM>>>
    // store topN PSMs of each scan of each .FT file
    std::unordered_map<std::string, std::unordered_map<int, std::vector<scanTopPSM>>> filesScansTopPSMs;
    std::vector<std::string> tokens;
    std::vector<std::string> Sep2PepFileChunkLines;
    sipPSM currentSipPSM;
    std::string currentFilePath;
    std::string currentLine;
    std::fstream sipFileStream;
    // for single file
    Spe2PepFileReader();
    // for multi files;
    Spe2PepFileReader(const std::string mWorkingPath);
    ~Spe2PepFileReader();
    std::vector<std::string> getSpe2PepFiles(const std::string &mWorkingPath);
    void splitString(const std::string &mString);
    // fill std::vectors in currentSipPSM
    void fillVectors();
    void readOneFile(std::string sipFileName);
    void readFileChunk();
    void readAllFiles();
    void readAllFilesTopPSMs();
    void fillScanTopPSMs(std::vector<scanTopPSM> &topPSMs, const int psmIX);
    // converted from filesScansTopPSMs
    // make it good for output
    sipPSM convertFilesScansTopPSMs();
    std::unordered_map<std::string, std::vector<std::string>>
    getFT2Spe2pepFileMap(const std::string &workingPath);
    void readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(const std::string &workingPath, size_t topN);
    void mergeDecoyToTarget(Spe2PepFileReader &targetReader, Spe2PepFileReader &decoyReader);
    void readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(
        const std::string &targetPath, const std::string &decoyPath, size_t topN = 5);
    void writeTSV(const std::string fileName);
};
