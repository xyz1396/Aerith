#include "sipFileReader.h"

sipFileReader::sipFileReader()
{
}

// get all .sip files' full path in workingPath
sipFileReader::sipFileReader(string mWorkingPath) : workingPath(mWorkingPath)
{
    string ext(".sip");
    // recursive_directory_iterator can get .sip file in all child path
    for (auto &p : fs::directory_iterator(workingPath))
    {
        if (p.path().extension() == ext)
            sipFileNames.push_back(p.path().string());
    }
    if (sipFileNames.size() == 0)
        cout << "no .sip file was found in " << workingPath << " !" << endl;
}

sipFileReader::~sipFileReader()
{
}

void sipFileReader::splitString(const string &mString)
{
    string sep = "\t";
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

void sipFileReader::readOneFile(string sipFileName)
{
    setlocale(LC_ALL, "C");
    ios_base::sync_with_stdio(false);
    if (fs::exists(sipFileName))
    {
        sipFileStream.open(sipFileName.c_str(), ios::in);
        if (!sipFileStream.is_open())
        {
            cout << "Cannot open " << sipFileName << endl;
        }
        string currentLine;
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
        cout << sipFileName << " does not exists" << endl;
    }
}

void sipFileReader::readAllFiles()
{
    for (string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        sipPSMs.push_back(currentSipPSM);
    }
}

void sipFileReader::readAllFilesTopPSMs()
{

    for (string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        auto fileIX = filesScansTopPSMs.find(currentSipPSM.fileName);
        // if not find the fileName, add fileName-fileScansTopPSMs pair
        if (fileIX == filesScansTopPSMs.end())
        {
            unordered_map<int, vector<topPSM>> fileScansTopPSMs;
            fileIX = filesScansTopPSMs.insert({currentSipPSM.fileName, fileScansTopPSMs}).first;
        }
        for (size_t i = 0; i < currentSipPSM.scanNumbers.size(); i++)
        {
            auto scanIX = fileIX->second.find(currentSipPSM.scanNumbers[i]);
            if (scanIX == fileIX->second.end())
            {
                vector<topPSM> scanTopPSMs;
                scanTopPSMs.reserve(topN);
                scanIX = fileIX->second.insert({currentSipPSM.scanNumbers[i], scanTopPSMs}).first;
            }
            fillScanTopPSMs(scanIX->second, i);
        }
    }
}

void sipFileReader::fillScanTopPSMs(vector<topPSM> &topPSMs, const int psmIX)
{
    size_t i = 0;
    while (i < topPSMs.size())
    {
        if (currentSipPSM.scores[psmIX] > topPSMs[i].score)
        {
            topPSMs.insert(topPSMs.begin() + i, topPSM({currentSipPSM.parentCharges[psmIX],
                                                        currentSipPSM.measuredParentMasses[psmIX],
                                                        currentSipPSM.calculatedParentMasses[psmIX],
                                                        currentSipPSM.searchName,
                                                        currentSipPSM.scores[psmIX],
                                                        currentSipPSM.identifiedPeptides[psmIX],
                                                        currentSipPSM.originalPeptides[psmIX],
                                                        currentSipPSM.proteinNames[psmIX]}));
            break;
        }
        i++;
    }
    // if insert fail and topPSMs has less than topN PSMs
    if (i == topPSMs.size() && i < topN)
        topPSMs.push_back(topPSM({currentSipPSM.parentCharges[psmIX],
                                  currentSipPSM.measuredParentMasses[psmIX],
                                  currentSipPSM.calculatedParentMasses[psmIX],
                                  currentSipPSM.searchName,
                                  currentSipPSM.scores[psmIX],
                                  currentSipPSM.identifiedPeptides[psmIX],
                                  currentSipPSM.originalPeptides[psmIX],
                                  currentSipPSM.proteinNames[psmIX]}));
    // if insert succeed and topPSMs already has topN PSMs
    if (topPSMs.size() > topN)
        topPSMs.pop_back();
}

sipPSM sipFileReader::convertFilesScansTopPSMs()
{
    // converted from filesScansTopPSMs
    sipPSM topPSMs;
    for (auto &fileIX : filesScansTopPSMs)
    {
        for (auto &scanIX : fileIX.second)
        {
            for (size_t i = 0; i < scanIX.second.size(); i++)
            {
                topPSMs.fileNames.push_back(fileIX.first);
                topPSMs.scanNumbers.push_back(scanIX.first);
                topPSMs.parentCharges.push_back(scanIX.second[i].parentCharge);
                topPSMs.measuredParentMasses.push_back(scanIX.second[i].measuredParentMass);
                topPSMs.calculatedParentMasses.push_back(scanIX.second[i].calculatedParentMass);
                topPSMs.searchNames.push_back(scanIX.second[i].searchName);
                topPSMs.ranks.push_back(i + 1);
                topPSMs.scores.push_back(scanIX.second[i].score);
                topPSMs.identifiedPeptides.push_back(scanIX.second[i].identifiedPeptide);
                topPSMs.originalPeptides.push_back(scanIX.second[i].originalPeptide);
                topPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
            }
        }
    }
    return topPSMs;
}