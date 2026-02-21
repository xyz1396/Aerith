#include "ftFileReader.h"

Scan::Scan() : scanNumber(0),
			   retentionTime(0),
			   precursorScanNumber(0),
			   precursorCharge(0),
			   isolationWindowCenterMZ(0) {}

Scan::Scan(int mScanNumber, float mRetentionTime, double mTIC, std::string mScanType, std::string mScanFilter) : scanNumber(mScanNumber),
																												 retentionTime(mRetentionTime),
																												 TIC(mTIC),
																												 scanType(mScanType),
																												 scanFilter(mScanFilter),
																												 precursorScanNumber(0),
																												 precursorCharge(0),
																												 isolationWindowCenterMZ(0) {}

Scan::Scan(int mScanNumber, float mRetentionTime, double mTIC,
		   std::string mScanType, std::string mScanFilter,
		   int mPrecursorScanNumber, int mPrecursorCharge, double mIsolationWindowCenterMZ,
		   std::vector<int> mPrecursorCharges,
		   std::vector<double> mPrecursorMZs) : scanNumber(mScanNumber),
												retentionTime(mRetentionTime),
												TIC(mTIC),
												scanType(mScanType),
												scanFilter(mScanFilter),
												precursorScanNumber(mPrecursorScanNumber),
												precursorCharge(mPrecursorCharge),
												isolationWindowCenterMZ(mIsolationWindowCenterMZ),
												precursorCharges(mPrecursorCharges),
												precursorMZs(mPrecursorMZs) {}

ftFileReader::ftFileReader()
{
}

ftFileReader::ftFileReader(std::string file) : ftFileName(file)
{
	setlocale(LC_ALL, "C");
	std::ios_base::sync_with_stdio(false);
	if (fs::exists(ftFileName))
	{
		ftFileStream.open(ftFileName.c_str(), std::ios::in);
		if (!ftFileStream.is_open())
		{
			isEmpty = true;
			std::cout << "Cannot open " << ftFileName << std::endl;
		}
		else
		{ // Set the buffer size to 2mb
			const int bufsize = 1024 * 1024 * 2;
			streamBuffer.resize(bufsize);
			ftFileStream.rdbuf()->pubsetbuf(streamBuffer.data(), streamBuffer.size());
		}
	}
	else
	{
		isEmpty = true;
		std::cout << ftFileName << " does not exists" << std::endl;
	}
}

ftFileReader::~ftFileReader()
{
	if (ftFileStream.is_open())
		ftFileStream.close();
}

void ftFileReader::splitString(const std::string &mString)
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

bool ftFileReader::detectPrecursorAndCharge()
{
	hasPrecursor = false;
	if (!getline(ftFileStream, currentLine))
		return false;
	while (!currentLine.empty() && currentLine[0] >= 'A' && currentLine[0] <= 'Z')
	{
		splitString(currentLine);
		if (tokens.size() > 2 && tokens[1] == "Instrument Model")
			instrument = tokens[2];
		else if (tokens.size() > 2 && tokens[1] == "ScanType")
		{
			scanType = tokens[2];
		}
		else if (tokens.size() > 2 && tokens[1] == "ScanFilter")
		{
			scanFilter = tokens[2];
		}
		else if (tokens.size() > 1 && tokens[1] == "ParentScanNumber")
			hasPrecursor = true;
		if (!getline(ftFileStream, currentLine))
		{
			currentLine.clear();
			break;
		}
	}

	// Check first peak format from the line where header parsing stopped.
	splitString(currentLine);
	size_t firstPeakTokenCount = tokens.size();

	// Return to the beginning of file for normal scan parsing.
	ftFileStream.clear();
	ftFileStream.seekg(0, std::ios::beg);
	if (!getline(ftFileStream, currentLine))
		return false;
	while (currentLine.empty() && getline(ftFileStream, currentLine))
	{
	}
	if (currentLine.empty())
		return false;

	if (firstPeakTokenCount == 2)
	{
		hasCharge = false;
		return true;
	}
	else if (firstPeakTokenCount == 6)
	{
		hasCharge = true;
		return true;
	}
	else
		return false;
}

bool ftFileReader::hasNextScan()
{
	return !ftFileStream.eof();
}

bool ftFileReader::isNumber(const std::string &s)
{
	if (s.empty())
		return false;
	return (s[0] >= '0' && s[0] <= '9');
}

Scan ftFileReader::readScanNumberRentionTime()
{
	continueRead = true;
	int mScanNumber = 0;
	float mRetentionTime = 0;
	double mTIC = 0;
	std::string mScanType = "";
	std::string mScanFilter = "";
	while (continueRead)
	{
		splitString(currentLine);
		// Check if tokens[0] is a number (spectrum data), if so, stop reading metadata
		if (isNumber(tokens[0]))
		{
			continueRead = false;
		}
		else
		{
			// Read metadata
			if (tokens[0] == "S")
			{
				mScanNumber = stoi(tokens[1]);
				mTIC = stod(tokens[2]);
			}
			else if (tokens.size() > 2 && tokens[1] == "RetentionTime")
			{
				mRetentionTime = stof(tokens[2]);
			}
			else if (tokens.size() > 2 && tokens[1] == "ScanType")
			{
				mScanType = tokens[2];
			}
			else if (tokens.size() > 2 && tokens[1] == "ScanFilter")
			{
				mScanFilter = tokens[2];
			}
		}
		getline(ftFileStream, currentLine);
	}
	return Scan(mScanNumber, mRetentionTime, mTIC, mScanType, mScanFilter);
}

int safe_stoi(const std::string &s)
{
	try
	{
		return std::stoi(s);
	}
	catch (const std::invalid_argument &e)
	{
		std::cerr << s << " is not a valid integer" << '\n';
		return 0;
	}
}

Scan ftFileReader::readScanNumberRentionTimePrecursor()
{
	continueRead = true;
	int mScanNumber = 0;
	double mIsolationWindowCenterMz = 0;
	int mPrecursorCharge = 0;
	float mRetentionTime = 0;
	double mTIC = 0;
	std::string mScanType = "";
	std::string mScanFilter = "";
	int mPrecursorScanNumber = 0;
	std::vector<int> mPrecursorCharges;
	std::vector<double> mPrecursorMZs;
	while (continueRead)
	{
		splitString(currentLine);
		// Stop when peak lins starts.
		if (isNumber(tokens[0]))
		{
			continueRead = false;
		}
		else if (tokens[0] == "D" && tokens.size() > 2)
		{
			mPrecursorScanNumber = stoi(tokens[2]);
		}
		else if (tokens[0] == "S" && tokens.size() > 3)
		{
			mScanNumber = stoi(tokens[1]);
			mIsolationWindowCenterMz = stod(tokens[2]);
			mTIC = stod(tokens[3]);
		}
		else if (tokens[0] == "Z" && tokens.size() > 1)
		{
			mPrecursorCharge = stoi(tokens[1]);
			// For DIA precursor and wide-window DDA, read isolation-window entries.
			for (size_t i = 3; i + 1 < tokens.size(); i += 2)
			{
				mPrecursorCharges.push_back(stoi(tokens[i]));
				mPrecursorMZs.push_back(stod(tokens[i + 1]));
			}
		}
		else if (tokens.size() > 2 && tokens[1] == "RetentionTime")
		{
			mRetentionTime = stof(tokens[2]);
		}
		else if (tokens.size() > 2 && tokens[1] == "ScanType")
		{
			mScanType = tokens[2];
		}
		else if (tokens.size() > 2 && tokens[1] == "ScanFilter")
		{
			mScanFilter = tokens[2];
		}
		getline(ftFileStream, currentLine);
	}

	return Scan(mScanNumber, mRetentionTime, mTIC,
				mScanType, mScanFilter, mPrecursorScanNumber,
				mPrecursorCharge, mIsolationWindowCenterMz,
				mPrecursorCharges, mPrecursorMZs);
}

void ftFileReader::readPeakCharge()
{
	while (hasNextScan())
	{
		if ((currentLine.at(0) >= '0') && (currentLine.at(0) <= '9'))
		{
			splitString(currentLine);
			if (hasCharge)
			{
				currentScan.mz.push_back(stod(tokens[0]));
				currentScan.intensity.push_back(stod(tokens[1]));
				// some times resolution is ∞
				currentScan.resolution.push_back(safe_stoi(tokens[2]));
				currentScan.baseLine.push_back(stof(tokens[3]));
				currentScan.signalToNoise.push_back(stof(tokens[4]));
				currentScan.charge.push_back(stoi(tokens[5]));
			}
			else
			{
				currentScan.mz.push_back(stod(tokens[0]));
				currentScan.intensity.push_back(stod(tokens[1]));
				currentScan.resolution.push_back(0);
				currentScan.baseLine.push_back(0);
				currentScan.signalToNoise.push_back(0);
				currentScan.charge.push_back(0);
			}
		}
		else
			break;
		getline(ftFileStream, currentLine);
	}
}

void ftFileReader::skipScans(const size_t scanNumber)
{
	size_t currentScanNumber = 0;
	while (hasNextScan())
	{
		if (currentLine[0] == 'S')
		{
			splitString(currentLine);
			currentScanNumber = stoi(tokens[1]);
			if (currentScanNumber >= scanNumber)
			{
				break;
			}
		}
		getline(ftFileStream, currentLine);
	}
}

void ftFileReader::readNextScan()
{
	if (hasPrecursor)
		currentScan = readScanNumberRentionTimePrecursor();
	else
	{
		currentScan = readScanNumberRentionTime();
	}
	readPeakCharge();
}

Scan ftFileReader::readOneScan(const size_t scanNumber)
{
	Scan emptyScan(0, 0.0f, 0.0, "", "", 0, 0, 0.0, std::vector<int>(), std::vector<double>());
	if (!detectPrecursorAndCharge())
	{
		std::cout << "Cannot read the first scan" << std::endl;
		return (emptyScan);
	}
	skipScans(scanNumber);
	if (hasNextScan())
	{
		readNextScan();
		if (currentScan.scanNumber == scanNumber)
			return currentScan;
		else
		{
			std::cout << scanNumber << " not exist" << std::endl;
			return emptyScan;
		}
	}
	else
	{
		std::cout << scanNumber << " is larger than the last scan number" << std::endl;
		return emptyScan;
	}
}

void ftFileReader::readScans(const size_t startScanNumber, const size_t endScanNumber)
{
	if (detectPrecursorAndCharge())
	{
		Scans.clear();
		skipScans(startScanNumber);
		while (hasNextScan())
		{
			readNextScan();
			Scans.push_back(currentScan);
			if (currentScan.scanNumber == endScanNumber)
			{
				break;
			}
			else if (currentScan.scanNumber > endScanNumber)
			{
				Scans.pop_back();
				break;
			}
		}
	}
	else
		std::cout << "Cannot read the first scan" << std::endl;
}

void ftFileReader::readScans(std::vector<size_t> &scanNumbers)
{
	// sort ascending
	std::sort(scanNumbers.begin(), scanNumbers.end());
	// remove duplicates
	auto new_end = std::unique(scanNumbers.begin(), scanNumbers.end());
	scanNumbers.erase(new_end, scanNumbers.end());
	// convert to set
	std::set<size_t> scanNumbersSet(scanNumbers.begin(), scanNumbers.end());
	if (detectPrecursorAndCharge())
	{
		size_t i = 0;
		Scans.clear();
		while (hasNextScan() && i < scanNumbers.size())
		{
			skipScans(scanNumbers[i]);
			readNextScan();
			if (scanNumbersSet.find(currentScan.scanNumber) != scanNumbersSet.end())
				Scans.push_back(currentScan);
			// skip scan not exist
			while (currentScan.scanNumber > scanNumbers[i] && i < scanNumbers.size())
			{
				std::cout << "Scan " << scanNumbers[i] << " not exist!" << std::endl;
				i++;
			}
			if (currentScan.scanNumber == scanNumbers[i] && i < scanNumbers.size())
				i++;
			if (currentScan.scanNumber >= scanNumbers.back())
				break;
		}
	}
	else
		std::cout << "Cannot read the first scan" << std::endl;
}

void ftFileReader::readAllScan()
{
	if (detectPrecursorAndCharge())
	{
		Scans.clear();
		while (hasNextScan())
		{
			readNextScan();
			Scans.push_back(currentScan);
		}
	}
	else
		std::cout << "Cannot read the first scan" << std::endl;
}

void ftFileReader::printFileInfo()
{
	std::cout << "ftFileName: " << ftFileName << std::endl;
	std::cout << "has precursor: " << hasPrecursor << std::endl;
	std::cout << "has charge: " << hasCharge << std::endl;
	std::cout << "scanType: " << scanType << std::endl;
	std::cout << "instrument: " << instrument << std::endl;
}

void ftFileReader::printScan(const Scan &mScan)
{
}
