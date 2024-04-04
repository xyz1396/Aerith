#include "sipFileReader.h"

sipFileReader::sipFileReader()
{
}

// get all .sip files' full path in workingPath
sipFileReader::sipFileReader(std::string mWorkingPath) : workingPath(mWorkingPath)
{
    std::string ext(".sip");
    // recursive_directory_iterator can get .sip file in all child path
    for (auto &p : fs::directory_iterator(workingPath))
    {
        if (p.path().extension() == ext)
            sipFileNames.push_back(p.path().string());
    }
    if (sipFileNames.size() == 0)
        std::cout << "no .sip file was found in " << workingPath << " !" << std::endl;
}

sipFileReader::~sipFileReader()
{
}

void sipFileReader::splitString(const std::string &mString)
{
    std::string sep = "\t";
    size_t start = 0;
    size_t end = mString.find(sep);
    tokens.clear();
    while (end != std::string::npos)
    {
        tokens.push_back(mString.substr(start, end - start));
        start = end + sep.length();
        end = mString.find(sep, start);
    }
    tokens.push_back(mString.substr(start));
}

void sipFileReader::fillVectors()
{
    currentSipPSM.scanNumbers.push_back(stoi(tokens[1]));
    currentSipPSM.parentCharges.push_back(stoi(tokens[2]));
    currentSipPSM.measuredParentMasses.push_back(stod(tokens[3]));
    currentSipPSM.calculatedParentMasses.push_back(stod(tokens[4]));
    currentSipPSM.ranks.push_back(stoi(tokens[8]));
    currentSipPSM.scores.push_back(stof(tokens[9]));
    currentSipPSM.identifiedPeptides.push_back(tokens[10]);
    currentSipPSM.originalPeptides.push_back(tokens[11]);
    currentSipPSM.proteinNames.push_back(tokens[12]);
}

void sipFileReader::readOneFile(std::string sipFileName)
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    if (fs::exists(sipFileName))
    {
        sipFileStream.open(sipFileName.c_str(), std::ios::in);
        if (!sipFileStream.is_open())
        {
            std::cout << "Cannot open " << sipFileName << std::endl;
        }
        std::string currentLine;
        // ignore annotation in sip result file
        while (!sipFileStream.eof())
        {
            getline(sipFileStream, currentLine);
            if (currentLine[0] != '#')
                break;
        }
        // read first line
        if (!sipFileStream.eof())
        {
            getline(sipFileStream, currentLine);
            splitString(currentLine);
            currentSipPSM.fileName = tokens[0];
            currentSipPSM.scanType = tokens[5];
            currentSipPSM.searchName = tokens[6];
            currentSipPSM.scoringFunction = tokens[7];
            fillVectors();
        }
        // read more lines
        // ignore the last /n
        while (sipFileStream.peek() != EOF)
        {
            getline(sipFileStream, currentLine);
            splitString(currentLine);
            fillVectors();
        }
        if (sipFileStream.is_open())
            sipFileStream.close();
    }
    else
    {
        std::cout << sipFileName << " does not exists" << std::endl;
    }
}

void sipFileReader::readAllFiles()
{
    for (std::string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        sipPSMs.push_back(currentSipPSM);
    }
}

void sipFileReader::readAllFilesTopPSMs()
{

    for (std::string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        auto fileIX = filesScansTopPSMs.find(currentSipPSM.fileName);
        // if not find the fileName, add fileName-fileScansTopPSMs pair
        if (fileIX == filesScansTopPSMs.end())
        {
            std::unordered_map<int, std::vector<scanTopPSM>> fileScansTopPSMs;
            fileIX = filesScansTopPSMs.insert({currentSipPSM.fileName, fileScansTopPSMs}).first;
        }
        for (size_t i = 0; i < currentSipPSM.scanNumbers.size(); i++)
        {
            auto scanIX = fileIX->second.find(currentSipPSM.scanNumbers[i]);
            if (scanIX == fileIX->second.end())
            {
                std::vector<scanTopPSM> scanTopPSMs;
                scanTopPSMs.reserve(topN);
                scanIX = fileIX->second.insert({currentSipPSM.scanNumbers[i], scanTopPSMs}).first;
            }
            fillScanTopPSMs(scanIX->second, i);
        }
    }
}

void sipFileReader::fillScanTopPSMs(std::vector<scanTopPSM> &scanTopPSMs, const int psmIX)
{
    // for insert position
    size_t i = 0;
    bool insertSucced = false;
    scanTopPSM mScanTopPSM = scanTopPSM(currentSipPSM.parentCharges[psmIX],
                                        currentSipPSM.measuredParentMasses[psmIX],
                                        currentSipPSM.calculatedParentMasses[psmIX],
                                        currentSipPSM.searchName,
                                        currentSipPSM.scores[psmIX],
                                        currentSipPSM.identifiedPeptides[psmIX],
                                        currentSipPSM.originalPeptides[psmIX],
                                        currentSipPSM.proteinNames[psmIX]);
    while (i < scanTopPSMs.size())
    {
        if (currentSipPSM.scores[psmIX] > scanTopPSMs[i].score)
        {
            scanTopPSMs.insert(scanTopPSMs.begin() + i, mScanTopPSM);
            insertSucced = true;
            break;
        }
        i++;
    }
    // if insert fail and topPSMs has less than topN PSMs
    if (i == scanTopPSMs.size() && i < topN)
    {
        scanTopPSMs.push_back(mScanTopPSM);
        insertSucced = true;
    }
    // if insert succeed
    if (insertSucced)
    {
        for (size_t j = 0; j < scanTopPSMs.size(); j++)
        {
            if (mScanTopPSM.identifiedPeptide == scanTopPSMs[j].identifiedPeptide)
            {
                // remove the inserted one if PSM with higher score and same squence has already existed
                if (j < i)
                {
                    scanTopPSMs.erase(scanTopPSMs.begin() + i);
                    break;
                }
                // remove the PSM with lower score and same squence as inserted one
                if (j > i)
                {
                    scanTopPSMs.erase(scanTopPSMs.begin() + j);
                    break;
                }
            }
        }
    }
    // if insert succeed and topPSMs already has topN PSMs
    if (scanTopPSMs.size() > topN)
        scanTopPSMs.pop_back();
}

sipPSM sipFileReader::convertFilesScansTopPSMs()
{
    // converted from filesScansTopPSMs
    sipPSM scanTopPSMs;
    for (auto &fileIX : filesScansTopPSMs)
    {
        for (auto &scanIX : fileIX.second)
        {
            for (size_t i = 0; i < scanIX.second.size(); i++)
            {
                scanTopPSMs.fileNames.push_back(fileIX.first);
                scanTopPSMs.scanNumbers.push_back(scanIX.first);
                scanTopPSMs.parentCharges.push_back(scanIX.second[i].parentCharge);
                scanTopPSMs.measuredParentMasses.push_back(scanIX.second[i].measuredParentMass);
                scanTopPSMs.calculatedParentMasses.push_back(scanIX.second[i].calculatedParentMass);
                scanTopPSMs.searchNames.push_back(scanIX.second[i].searchName);
                scanTopPSMs.ranks.push_back(i + 1);
                scanTopPSMs.scores.push_back(scanIX.second[i].score);
                scanTopPSMs.identifiedPeptides.push_back(scanIX.second[i].identifiedPeptide);
                scanTopPSMs.originalPeptides.push_back(scanIX.second[i].originalPeptide);
                scanTopPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
            }
        }
    }
    return scanTopPSMs;
}