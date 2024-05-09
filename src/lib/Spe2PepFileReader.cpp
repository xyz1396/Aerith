#include "Spe2PepFileReader.h"

Spe2PepFileReader::Spe2PepFileReader()
{
}

// get all .Spe2Pep files' full path in workingPath
Spe2PepFileReader::Spe2PepFileReader(std::string mWorkingPath) : workingPath(mWorkingPath)
{
    // recursive_directory_iterator can get .Spe2Pep file in all child path
    for (auto &p : fs::directory_iterator(workingPath))
    {
        std::string stem = p.path().stem().string();
        std::string ext = p.path().extension().string();
        // for ".Spe2Pep.txt"
        if (stem.size() > 8)
        {
            if (ext == ".txt" && stem.substr(stem.size() - 8) == ".Spe2Pep")
                sipFileNames.push_back(p.path().string());
        }
    }
    if (sipFileNames.size() == 0)
        std::cout << "no .Spe2Pep.txt file was found in " << workingPath << " !" << std::endl;
}

Spe2PepFileReader::~Spe2PepFileReader()
{
}

void Spe2PepFileReader::splitString(const std::string &mString)
{
    std::string sep = "\t";
    size_t start = 0;
    size_t end = mString.find(sep);
    tokens.clear();
    tokens.reserve(10);
    while (end != std::string::npos)
    {
        tokens.push_back(mString.substr(start, end - start));
        start = end + sep.length();
        end = mString.find(sep, start);
    }
    tokens.push_back(mString.substr(start));
}

void Spe2PepFileReader::fillVectors()
{
    size_t i = 0;
    splitString(Sep2PepFileChunkLines[i]);
    int scanNumber = stoi(tokens[2]);
    double isolationWindowCenter = stod(tokens[4]);
    float retentionTime = stof(tokens[7]);
    int precursorScanNumber = stoi(tokens[8]);
    i++;
    while (i < Sep2PepFileChunkLines.size())
    {
        splitString(Sep2PepFileChunkLines[i]);
        currentSipPSM.scanNumbers.push_back(scanNumber);
        currentSipPSM.isolationWindowCenterMZs.push_back(isolationWindowCenter);
        currentSipPSM.retentionTimes.push_back(retentionTime);
        currentSipPSM.precursorScanNumbers.push_back(precursorScanNumber);
        currentSipPSM.parentCharges.push_back(stoi(tokens[8]));
        currentSipPSM.measuredParentMasses.push_back(stod(tokens[9]));
        currentSipPSM.calculatedParentMasses.push_back(stod(tokens[3]));
        currentSipPSM.ranks.push_back(i);
        currentSipPSM.MVHscores.push_back(stof(tokens[4]));
        currentSipPSM.XcorrScores.push_back(stof(tokens[5]));
        currentSipPSM.WDPscores.push_back(stof(tokens[6]));
        currentSipPSM.identifiedPeptides.push_back(tokens[1]);
        currentSipPSM.originalPeptides.push_back(tokens[2]);
        currentSipPSM.proteinNames.push_back(tokens[7]);
        i++;
    }
}

void Spe2PepFileReader::readFileChunk()
{
    Sep2PepFileChunkLines.clear();
    Sep2PepFileChunkLines.reserve(topN + 1);
    // line begin with "+"
    Sep2PepFileChunkLines.push_back(currentLine);
    while (!sipFileStream.eof())
    {
        getline(sipFileStream, currentLine);
        if (currentLine[0] != '*')
            break;
        Sep2PepFileChunkLines.push_back(currentLine);
    }
}

void Spe2PepFileReader::readOneFile(std::string sipFileName)
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
        // ignore annotation in .Spe2Pep.txt file
        while (!sipFileStream.eof())
        {
            getline(sipFileStream, currentLine);
            if (currentLine[0] != '#')
                break;
        }
        // ignore column names
        readFileChunk();
        // read first chunck
        if (!sipFileStream.eof())
        {
            readFileChunk();
            splitString(Sep2PepFileChunkLines[0]);
            currentSipPSM.fileName = tokens[1];
            currentSipPSM.scanType = tokens[5];
            currentSipPSM.searchName = tokens[6];
            fillVectors();
        }
        // read more chunks
        // ignore the last /n
        while (sipFileStream.peek() != EOF)
        {
            readFileChunk();
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

void Spe2PepFileReader::readAllFiles()
{
    for (std::string sipFileName : sipFileNames)
    {
        currentSipPSM = sipPSM();
        readOneFile(sipFileName);
        sipPSMs.push_back(currentSipPSM);
    }
}

void Spe2PepFileReader::readAllFilesTopPSMs()
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

void Spe2PepFileReader::fillScanTopPSMs(std::vector<scanTopPSM> &scanTopPSMs, const int psmIX)
{
    // for insert position
    size_t i = 0;
    bool insertSucced = false;
    scanTopPSM mScanTopPSM = scanTopPSM(currentSipPSM.parentCharges[psmIX],
                                        currentSipPSM.precursorScanNumbers[psmIX],
                                        currentSipPSM.isolationWindowCenterMZs[psmIX],
                                        currentSipPSM.measuredParentMasses[psmIX],
                                        currentSipPSM.calculatedParentMasses[psmIX],
                                        currentSipPSM.searchName,
                                        currentSipPSM.identifiedPeptides[psmIX],
                                        currentSipPSM.originalPeptides[psmIX],
                                        currentSipPSM.proteinNames[psmIX],
                                        currentSipPSM.retentionTimes[psmIX],
                                        currentSipPSM.MVHscores[psmIX],
                                        currentSipPSM.XcorrScores[psmIX],
                                        currentSipPSM.WDPscores[psmIX]);
    while (i < scanTopPSMs.size())
    {
        // according to the first score
        if (currentSipPSM.MVHscores[psmIX] > scanTopPSMs[i].MVHscore)
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

sipPSM Spe2PepFileReader::convertFilesScansTopPSMs()
{
    // converted from filesScansTopPSMs
    sipPSM topPSMs;
    std::map<int, std::vector<scanTopPSM>> orderedMap;
    for (auto &fileIX : filesScansTopPSMs)
    {
        // sort by scan number
        orderedMap = std::map<int, std::vector<scanTopPSM>>(fileIX.second.begin(), fileIX.second.end());
        for (auto &scanIX : orderedMap)
        {
            for (size_t i = 0; i < scanIX.second.size(); i++)
            {
                topPSMs.fileNames.push_back(fileIX.first);
                topPSMs.scanNumbers.push_back(scanIX.first);
                topPSMs.precursorScanNumbers.push_back(scanIX.second[i].precursorScanNumber);
                topPSMs.parentCharges.push_back(scanIX.second[i].parentCharge);
                topPSMs.isolationWindowCenterMZs.push_back(scanIX.second[i].isolationWindowCenterMZ);
                topPSMs.measuredParentMasses.push_back(scanIX.second[i].measuredParentMass);
                topPSMs.calculatedParentMasses.push_back(scanIX.second[i].calculatedParentMass);
                topPSMs.searchNames.push_back(scanIX.second[i].searchName);
                topPSMs.ranks.push_back(i + 1);
                topPSMs.scores.push_back(scanIX.second[i].score);
                topPSMs.identifiedPeptides.push_back(scanIX.second[i].identifiedPeptide);
                topPSMs.originalPeptides.push_back(scanIX.second[i].originalPeptide);
                topPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
                topPSMs.retentionTimes.push_back(scanIX.second[i].retentionTime);
                topPSMs.MVHscores.push_back(scanIX.second[i].MVHscore);
                topPSMs.XcorrScores.push_back(scanIX.second[i].XcorrScore);
                topPSMs.WDPscores.push_back(scanIX.second[i].WDPscore);
            }
        }
    }
    return topPSMs;
}

void Spe2PepFileReader::readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(const std::string &workingPath, size_t topN = 5)
{
    Spe2PepFileReader reader(workingPath);
    std::unordered_map<std::string, std::vector<std::string>> FT2map;
    size_t pos;
    // group the .Spe2PepFile.txt files by FT2 name
    for (const auto &str : reader.sipFileNames)
    {
        // Extract the filename
        std::filesystem::path full_path(str);
        std::string filename = full_path.filename().string();
        // get FT2 file name by substring before the last third "." character
        pos = filename.find_last_of(".");
        pos = filename.find_last_of(".", pos - 1);
        pos = filename.find_last_of(".", pos - 1);
        if (pos == std::string::npos)
            std::cout << filename << " 's name style is not right" << std::endl;
        else
            FT2map[filename.substr(0, pos)].push_back(str);
    }
    FT2s.clear();
    for (auto const &pair : FT2map)
    {
        FT2s.push_back(pair.first);
    }
    sipPSMs = std::vector<sipPSM>(FT2s.size());
    int num_cores = omp_get_num_procs();
    // Set the number of threads to the lesser of num_cores and 10
    int num_threads = std::min(num_cores, 10);
    omp_set_num_threads(num_threads);
#pragma omp parallel for
    for (size_t i = 0; i < FT2s.size(); i++)
    {
        Spe2PepFileReader reader;
        reader.sipFileNames = FT2map[FT2s[i]];
        reader.topN = topN;
        reader.readAllFilesTopPSMs();
        sipPSMs[i] = reader.convertFilesScansTopPSMs();
    }
}

void Spe2PepFileReader::writeTSV(const std::string fileName = "a.tsv")
{
    setlocale(LC_ALL, "C");
    std::ios_base::sync_with_stdio(false);
    std::ofstream file(fileName);
    if (!file)
    {
        std::cerr << "Unable to open file for writing.\n";
    }
    file << std::fixed << std::setprecision(6);

    const size_t chunkSize = 10000;
    std::stringstream ss;
    ss << "fileNames"
       << "\t"
       << "scanNumbers"
       << "\t"
       << "precursorScanNumbers"
       << "\t"
       << "parentCharges"
       << "\t"
       << "isolationWindowCenterMZs"
       << "\t"
       << "measuredParentMasses"
       << "\t"
       << "calculatedParentMasses"
       << "\t"
       << "searchNames"
       << "\t"
       << "retentionTimes"
       << "\t"
       << "MVHscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "WDPscores"
       << "\t"
       << "ranks"
       << "\t"
       << "identifiedPeptides"
       << "\t"
       << "originalPeptides"
       << "\t"
       << "proteinNames"
       << "\n";
    for (size_t i = 0; i < sipPSMs.size(); i++)
    {
        for (size_t j = 0; j < sipPSMs[i].fileNames.size(); j += chunkSize)
        {
            for (size_t k = j; k < std::min(j + chunkSize, sipPSMs[i].fileNames.size()); k++)
            {
                ss << sipPSMs[i].fileNames[k] << "\t"
                   << sipPSMs[i].scanNumbers[k] << "\t"
                   << sipPSMs[i].precursorScanNumbers[k] << "\t"
                   << sipPSMs[i].parentCharges[k] << "\t"
                   << sipPSMs[i].isolationWindowCenterMZs[k] << "\t"
                   << sipPSMs[i].measuredParentMasses[k] << "\t"
                   << sipPSMs[i].calculatedParentMasses[k] << "\t"
                   << sipPSMs[i].searchNames[k] << "\t"
                   << sipPSMs[i].retentionTimes[k] << "\t"
                   << sipPSMs[i].MVHscores[k] << "\t"
                   << sipPSMs[i].XcorrScores[k] << "\t"
                   << sipPSMs[i].WDPscores[k] << "\t"
                   << sipPSMs[i].ranks[k] << "\t"
                   << sipPSMs[i].identifiedPeptides[k] << "\t"
                   << sipPSMs[i].originalPeptides[k] << "\t"
                   << sipPSMs[i].proteinNames[k] << "\n";
            }
            file << ss.str();
            // Clear the stringstream
            ss.str(std::string());
            ss.clear();
        }
    }
    file.close();
}