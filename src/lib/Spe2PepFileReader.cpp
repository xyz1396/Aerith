#include "Spe2PepFileReader.h"

Spe2PepFileReader::Spe2PepFileReader()
{
}

// get all .Spe2Pep files' full path in workingPath
Spe2PepFileReader::Spe2PepFileReader(const std::string mWorkingPath) : workingPath(mWorkingPath)
{
    sipFileNames = getSpe2PepFiles(mWorkingPath);
}

Spe2PepFileReader::~Spe2PepFileReader()
{
}

std::vector<std::string> Spe2PepFileReader::getSpe2PepFiles(const std::string &mWorkingPath)
{
    std::vector<std::string> spe2PepFiles;
    // recursive_directory_iterator can get .Spe2Pep file in all child path
    for (auto &p : fs::directory_iterator(mWorkingPath))
    {
        std::string stem = p.path().stem().string();
        std::string ext = p.path().extension().string();
        // for ".Spe2Pep.txt"
        if (stem.size() > 8)
        {
            if (ext == ".txt" && stem.substr(stem.size() - 8) == ".Spe2Pep")
                spe2PepFiles.push_back(p.path().string());
        }
    }
    if (spe2PepFiles.size() == 0)
        std::cerr << "no .Spe2Pep.txt file was found in " << mWorkingPath << " !" << std::endl;
    return spe2PepFiles;
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

void Spe2PepFileReader::splitStringView(std::string_view str, char delimiter)
{
    tokens.clear();
    size_t start = 0;
    size_t end;
    while ((end = str.find(delimiter, start)) != std::string_view::npos)
    {
        tokens.emplace_back(str.substr(start, end - start));
        start = end + 1;
    }
    tokens.emplace_back(str.substr(start));
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
    std::size_t start, end;
    std::string nakeSeq;
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
        currentSipPSM.WDPscores.push_back(stof(tokens[4]));
        currentSipPSM.XcorrScores.push_back(stof(tokens[5]));
        currentSipPSM.MVHscores.push_back(stof(tokens[6]));
        currentSipPSM.identifiedPeptides.push_back(tokens[1]);
        currentSipPSM.originalPeptides.push_back(tokens[2]);
        // init isDecoys
        currentSipPSM.isDecoys.push_back(false);

        start = tokens[1].find_first_of('[');
        end = tokens[1].find_last_of(']');
        nakeSeq = tokens[1].substr(start + 1, end - start - 1);
        currentSipPSM.nakePeptides.push_back(nakeSeq);

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

void Spe2PepFileReader::readOneEntireFile(const std::string &sipFileName)
{
    // setlocale(LC_ALL, "C");
    // std::ios_base::sync_with_stdio(false);
    std::ifstream file(sipFileName, std::ios::in | std::ios::binary | std::ios::ate);
    if (!file)
    {
        std::cout << "Cannot open " << sipFileName << std::endl;
        return;
    }
    std::streamsize fileSize = file.tellg();
    std::string fileContent(fileSize, '\0');
    file.seekg(0, std::ios::beg);
    // Read the entire file into one string
    if (!file.read(&fileContent[0], fileSize))
    {
        std::cout << "Cannot open " << sipFileName << std::endl;
        return;
    }
    size_t pos = 0;
    size_t fileEnd = fileContent.size();

    // Skip annotations starting with '#'
    while (pos < fileEnd && fileContent[pos] == '#')
    {
        pos = fileContent.find('\n', pos);
        if (pos == std::string::npos)
            break;
        pos++;
    }

    // Skip column names
    pos = fileContent.find('\n', pos);
    if (pos != std::string::npos)
        pos++;
    pos = fileContent.find('\n', pos);
    if (pos != std::string::npos)
        pos++;

    int scanNumber;
    double isolationWindowCenter;
    float retentionTime;
    int precursorScanNumber;
    std::string nakeSeq;
    size_t rank = 1, start, end;
    std::string_view line;

    // read frist line begin with +
    size_t lineEnd = fileContent.find('\n', pos);
    if (pos != std::string::npos)
    {
        line = std::string_view(&fileContent[pos], lineEnd - pos);
        splitStringView(line, '\t');
        currentSipPSM.fileName = tokens[1];
        currentSipPSM.scanType = tokens[5];
        currentSipPSM.searchName = tokens[6];
        scanNumber = stoi(tokens[2]);
        isolationWindowCenter = stod(tokens[4]);
        retentionTime = stof(tokens[7]);
        precursorScanNumber = stoi(tokens[8]);
        pos++;
    }

    // Process lines
    while (pos < fileEnd)
    {
        // Find the end of the current line
        size_t lineEnd = fileContent.find('\n', pos);
        if (lineEnd == std::string::npos)
            lineEnd = fileEnd;

        // Create a string_view for the line
        line = std::string_view(&fileContent[pos], lineEnd - pos);

        // Process the line based on its prefix
        if (!line.empty())
        {
            splitStringView(line, '\t');
            if (line[0] == '+')
            {
                scanNumber = stoi(tokens[2]);
                isolationWindowCenter = stod(tokens[4]);
                retentionTime = stof(tokens[7]);
                precursorScanNumber = stoi(tokens[8]);
                rank = 1;
            }
            else if (line[0] == '*')
            {
                currentSipPSM.scanNumbers.push_back(scanNumber);
                currentSipPSM.isolationWindowCenterMZs.push_back(isolationWindowCenter);
                currentSipPSM.retentionTimes.push_back(retentionTime);
                currentSipPSM.precursorScanNumbers.push_back(precursorScanNumber);
                currentSipPSM.parentCharges.push_back(stoi(tokens[8]));
                currentSipPSM.measuredParentMasses.push_back(stod(tokens[9]));
                currentSipPSM.calculatedParentMasses.push_back(stod(tokens[3]));
                currentSipPSM.ranks.push_back(rank);
                rank++;
                currentSipPSM.WDPscores.push_back(stof(tokens[4]));
                currentSipPSM.XcorrScores.push_back(stof(tokens[5]));
                currentSipPSM.MVHscores.push_back(stof(tokens[6]));
                currentSipPSM.identifiedPeptides.push_back(tokens[1]);
                currentSipPSM.originalPeptides.push_back(tokens[2]);
                // init isDecoys
                currentSipPSM.isDecoys.push_back(false);

                start = tokens[1].find_first_of('[');
                end = tokens[1].find_last_of(']');
                nakeSeq = tokens[1].substr(start + 1, end - start - 1);
                currentSipPSM.nakePeptides.push_back(nakeSeq);

                currentSipPSM.proteinNames.push_back(tokens[7]);
            }
        }

        // Move to the next line
        pos = lineEnd + 1;
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
        readOneEntireFile(sipFileName);
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
                                        currentSipPSM.nakePeptides[psmIX],
                                        currentSipPSM.proteinNames[psmIX],
                                        currentSipPSM.retentionTimes[psmIX],
                                        currentSipPSM.MVHscores[psmIX],
                                        currentSipPSM.XcorrScores[psmIX],
                                        currentSipPSM.WDPscores[psmIX],
                                        currentSipPSM.isDecoys[psmIX]);
    while (i < scanTopPSMs.size())
    {
        // // according to the first score
        // if (currentSipPSM.MVHscores[psmIX] > scanTopPSMs[i].MVHscore)
        // according to the WDP score
        if (currentSipPSM.WDPscores[psmIX] > scanTopPSMs[i].WDPscore)
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
            // if (mScanTopPSM.originalPeptide == scanTopPSMs[j].originalPeptide)
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
                topPSMs.nakePeptides.push_back(scanIX.second[i].nakePeptide);
                topPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
                topPSMs.retentionTimes.push_back(scanIX.second[i].retentionTime);
                topPSMs.MVHscores.push_back(scanIX.second[i].MVHscore);
                topPSMs.XcorrScores.push_back(scanIX.second[i].XcorrScore);
                topPSMs.WDPscores.push_back(scanIX.second[i].WDPscore);
                topPSMs.isDecoys.push_back(scanIX.second[i].isDecoy);
            }
        }
    }
    return topPSMs;
}

std::unordered_map<std::string, std::vector<std::string>> Spe2PepFileReader::
    getFT2Spe2pepFileMap(const std::string &workingPath)
{
    std::unordered_map<std::string, std::vector<std::string>> FT2map;
    size_t pos;
    std::vector<std::string> spe2PepFiles = getSpe2PepFiles(workingPath);
    // group the .Spe2PepFile.txt files by FT2 name
    for (const auto &str : spe2PepFiles)
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
    return FT2map;
}

void Spe2PepFileReader::readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(const std::string &mWorkingPath,
                                                                        size_t topN = 5)
{
    std::unordered_map<std::string, std::vector<std::string>> FT2map = getFT2Spe2pepFileMap(mWorkingPath);
    FT2s.reserve(FT2map.size());
    for (auto const &pair : FT2map)
    {
        FT2s.push_back(pair.first);
    }
    sipPSMs = std::vector<sipPSM>(FT2s.size());
#ifdef _OPENMP
    int num_cores = omp_get_num_procs();
    // Set the number of threads to the lesse of num_cores and 10
    int num_threads = std::min(num_cores, 10);
    omp_set_num_threads(num_threads);
#endif
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

void Spe2PepFileReader::mergeDecoyToTarget(Spe2PepFileReader &targetReader,
                                           Spe2PepFileReader &decoyReader)
{
    std::unordered_set<std::string> targetSeqs;
    targetSeqs.reserve(10000);
    for (auto &targetFileIX : targetReader.filesScansTopPSMs)
    {
        for (auto &targetScanIX : targetFileIX.second)
        {
            for (auto &psm : targetScanIX.second)
            {
                psm.isDecoy = false;
                targetSeqs.insert(psm.nakePeptide);
            }
        }
    }
    for (auto &decoyFileIX : decoyReader.filesScansTopPSMs)
    {
        for (auto &decoyScanIX : decoyFileIX.second)
        {
            for (size_t i = 0; i < decoyScanIX.second.size(); i++)
            {
                if (targetSeqs.find(decoyScanIX.second[i].nakePeptide) == targetSeqs.end())
                {
                    decoyScanIX.second[i].isDecoy = true;
                }
                else
                    // remove decoy with the same pep seq as target
                    decoyScanIX.second.erase(decoyScanIX.second.begin() + i);
            }
        }
    }
    for (auto &decoyFileIX : decoyReader.filesScansTopPSMs)
    {
        auto targetFileIX = targetReader.filesScansTopPSMs.find(decoyFileIX.first);
        if (targetFileIX != targetReader.filesScansTopPSMs.end())
        {
            for (auto &decoyScanIX : decoyFileIX.second)
            {
                auto targetScanIX = targetFileIX->second.find(decoyScanIX.first);
                if (targetScanIX != targetFileIX->second.end())
                {
                    targetScanIX->second.insert(targetScanIX->second.end(),
                                                decoyScanIX.second.begin(), decoyScanIX.second.end());
                }
                else
                {
                    targetFileIX->second.insert(decoyScanIX);
                }
            }
        }
        else
        {
            targetReader.filesScansTopPSMs.insert(decoyFileIX);
            std::cerr << "target decoy files not match, " << decoyFileIX.first << " only exists in decoy" << std::endl;
        }
    }
    for (auto &targetFileIX : targetReader.filesScansTopPSMs)
    {
        for (auto &targetScanIX : targetFileIX.second)
        {
            std::sort(targetScanIX.second.begin(), targetScanIX.second.end(), [](const scanTopPSM &psm1, const scanTopPSM &psm2)
                    //   { return psm1.MVHscore > psm2.MVHscore; });
              { return psm1.WDPscore > psm2.WDPscore; });
        }
    }
}

void Spe2PepFileReader::readSpe2PepFilesScansTopPSMsFromEachFT2TargetAndDecoyParallel(
    const std::string &targetPath, const std::string &decoyPath, size_t topN)
{
    std::unordered_map<std::string, std::vector<std::string>> targetFT2map =
        getFT2Spe2pepFileMap(targetPath);
    std::unordered_map<std::string, std::vector<std::string>> decoyFT2map =
        getFT2Spe2pepFileMap(decoyPath);
    FT2s.reserve(targetFT2map.size());
    for (auto const &pair : targetFT2map)
    {
        FT2s.push_back(pair.first);
    }
    sipPSMs = std::vector<sipPSM>(FT2s.size());
#ifdef _OPENMP
    int num_cores = omp_get_num_procs();
    // Set the number of threads to the lesse of num_cores and 10
    int num_threads = std::min(num_cores, 10);
    omp_set_num_threads(num_threads);
#endif
#pragma omp parallel for
    for (size_t i = 0; i < FT2s.size(); i++)
    {
        Spe2PepFileReader targetReader;
        targetReader.sipFileNames = targetFT2map[FT2s[i]];
        targetReader.topN = topN;
        targetReader.readAllFilesTopPSMs();
        Spe2PepFileReader decoyReader;
        decoyReader.sipFileNames = decoyFT2map[FT2s[i]];
        decoyReader.topN = topN;
        decoyReader.readAllFilesTopPSMs();
        mergeDecoyToTarget(targetReader, decoyReader);
        sipPSMs[i] = targetReader.convertFilesScansTopPSMs();
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
       << "WDPscores"
       << "\t"
       << "XcorrScores"
       << "\t"
       << "MVHscores"              
       << "\t"
       << "ranks"
       << "\t"
       << "identifiedPeptides"
       << "\t"
       << "originalPeptides"
       << "\t"
       << "nakePeptides"
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
                   << sipPSMs[i].WDPscores[k] << "\t"
                   << sipPSMs[i].XcorrScores[k] << "\t"
                   << sipPSMs[i].MVHscores[k] << "\t"
                   << sipPSMs[i].ranks[k] << "\t"
                   << sipPSMs[i].identifiedPeptides[k] << "\t"
                   << sipPSMs[i].originalPeptides[k] << "\t"
                   << sipPSMs[i].nakePeptides[k] << "\t"
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