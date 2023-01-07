#include "ftFileReader.h"

Scan::Scan() : scanNumber(0),
			   retentionTime(0),
			   precursorScanNumber(0),
			   precursorMz(0),
			   precursorCharge(0) {}

Scan::Scan(int mScanNumber, float mRetentionTime, double mTIC) : scanNumber(mScanNumber),
																 retentionTime(mRetentionTime),
																 TIC(mTIC),
																 precursorScanNumber(0),
																 precursorMz(0),
																 precursorCharge(0) {}

Scan::Scan(int mScanNumber, float mRetentionTime, double mTIC, int mPrecusorScanNumber,
		   double mPrecusorMz, int mPrecusorCharge) : scanNumber(mScanNumber),
													  retentionTime(mRetentionTime),
													  TIC(mTIC),
													  precursorScanNumber(mPrecusorScanNumber),
													  precursorMz(mPrecusorMz),
													  precursorCharge(mPrecusorCharge) {}

ftFileReader::ftFileReader(string file) : ftFileName(file)
{
	setlocale(LC_ALL, "C");
	ios_base::sync_with_stdio(false);
	if (fs::exists(ftFileName))
	{
		ftFileStream.open(ftFileName.c_str(), ios::in);
		if (!ftFileStream.is_open())
		{
			isEmpty = true;
			cout << "Cannot open " << ftFileName << endl;
		}
	}
	else
	{
		isEmpty = true;
		cout << ftFileName << " does not exists" << endl;
	}
}

ftFileReader::~ftFileReader()
{
	if (ftFileStream.is_open())
		ftFileStream.close();
}

void ftFileReader::splitString(const string &mString)
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

bool ftFileReader::detectPrecursorAndCharge()
{
	hasPrecursor = false;
	getline(ftFileStream, currentLine);
	while (currentLine[0] >= 'A' && currentLine[0] <= 'Z')
	{
		splitString(currentLine);
		if (tokens[1] == "Instrument Model")
			instrument = tokens[2];
		else if (tokens[1] == "ScanType")
		{
			scanType = tokens[2];
		}
		else if (tokens[1] == "ScanFilter")
		{
			scanFilter = tokens[2];
		}
		else if (tokens[1] == "ParentScanNumber")
			hasPrecursor = true;
		getline(ftFileStream, currentLine);
	}
	// return to beginning of line of first peak
	ftFileStream.seekg(0);
	splitString(currentLine);
	if (tokens.size() == 2)
	{
		hasCharge = false;
		return true;
	}
	else if (tokens.size() == 6)
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

Scan ftFileReader::readScanNumberRentionTime()
{
	continueRead = true;
	int mScanNumber = 0;
	float mRetentionTime = 0;
	double mTIC = 0;
	while (continueRead)
	{
		splitString(currentLine);
		if (tokens[0] != "I")
		{
			if (tokens[0] == "S")
			{
				mScanNumber = stoi(tokens[1]);
				mTIC = stod(tokens[2]);
			}
		}
		else
		{
			continueRead = false;
			mRetentionTime = stof(tokens[2]);
			// read 3 more lines to the first peak
			getline(ftFileStream, currentLine);
			getline(ftFileStream, currentLine);
		}
		getline(ftFileStream, currentLine);
	}
	return Scan(mScanNumber, mRetentionTime, mTIC);
}

Scan ftFileReader::readScanNumberRentionTimePrecursor()
{
	continueRead = true;
	int mScanNumber = 0;
	double mPrecusorMz = 0;
	double mPrecusorCharge = 0;
	float mRetentionTime = 0;
	double mTIC = 0;
	int mPrecusorScanNumber = 0;
	while (continueRead)
	{
		splitString(currentLine);
		if (tokens[0] != "D")
		{
			if (tokens[0] == "S")
			{
				mScanNumber = stoi(tokens[1]);
				mPrecusorMz = stod(tokens[2]);
				mTIC = stod(tokens[3]);
			}
			else if (tokens[0] == "Z")
				mPrecusorCharge = stoi(tokens[1]);
			else if (tokens[1] == "RetentionTime")
				mRetentionTime = stof(tokens[2]);
		}
		else
		{
			continueRead = false;
			mPrecusorScanNumber = stoi(tokens[2]);
		}
		getline(ftFileStream, currentLine);
	}
	return Scan(mScanNumber, mRetentionTime, mTIC, mPrecusorScanNumber,
				mPrecusorMz, mPrecusorCharge);
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
				currentScan.resolution.push_back(stoi(tokens[2]));
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

Scan ftFileReader::readOneScan(const int scanCount)
{
	if (detectPrecursorAndCharge())
	{
		for (int i = 0; i < scanCount; i++)
		{
			if (hasNextScan())
				readNextScan();
			else
			{
				cout << ftFileName << " has " << i << " scans < " << scanCount << endl;
				Scan emptyScan(0, 0, 0, 0, 0, 0);
				return emptyScan;
			}
		}
	}
	else
		cout << "Cannot read the first scan" << endl;
	return currentScan;
}

void ftFileReader::readScans(const int scanCount)
{
	if (detectPrecursorAndCharge())
	{
		Scans.clear();
		for (int i = 0; i < scanCount; i++)
		{
			if (hasNextScan())
			{
				readNextScan();
				Scans.push_back(currentScan);
			}
			else
			{
				cout << ftFileName << " has " << i << " scans < " << scanCount << endl;
				break;
			}
		}
	}
	else
		cout << "Cannot read the first scan" << endl;
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
		cout << "Cannot read the first scan" << endl;
}

void ftFileReader::printFileInfo()
{
	cout << "ftFileName: " << ftFileName << endl;
	cout << "has precursor: " << hasPrecursor << endl;
	cout << "has charge: " << hasCharge << endl;
	cout << "scanType: " << scanType << endl;
	cout << "instrument: " << instrument << endl;
}

void ftFileReader::printScan(const Scan &mScan)
{
}
