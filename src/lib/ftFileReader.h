#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <filesystem>
#include <algorithm>
namespace fs = std::filesystem;

struct alignas(64) Scan
{
	size_t scanNumber;
	float retentionTime;
	double TIC;
	std::vector<double> mz;
	std::vector<double> mass;
	std::vector<double> intensity;

	// only MS2 scan has follows
	int precursorScanNumber;
	double precursorMz;
	int precursorCharge;

	// only orbitrap scan has follows
	std::vector<int> resolution;
	std::vector<float> baseLine;
	std::vector<float> signalToNoise;
	std::vector<int> charge;

	Scan();
	// for MS1 scans
	Scan(int mScanNumber, float mRetentionTime, double mTIC);
	// for MS2 scans
	Scan(int mScanNumber, float mRetentionTime, double mTIC, int mPrecusorScanNumber,
		 double mPrecusorMz, int mPrecusorCharge);
};

class ftFileReader
{
private:
public:
	std::string ftFileName;
	std::ifstream ftFileStream;
	std::string currentLine;
	bool continueRead;
	bool hasPrecursor;
	bool hasCharge;
	// avoid empty file
	bool isEmpty = false;
	std::string instrument;
	std::string scanType;
	std::string scanFilter;
	Scan currentScan;
	std::vector<Scan> Scans;
	std::vector<std::string> tokens;
	ftFileReader(std::string file);
	~ftFileReader();
	void splitString(const std::string &mString);
	bool detectPrecursorAndCharge();
	bool hasNextScan();
	Scan readScanNumberRentionTime();
	Scan readScanNumberRentionTimePrecursor();
	void readPeakCharge();
	// ignore scans before scanNumber
	void skipScans(const size_t scanNumber);
	void readNextScan();
	Scan readOneScan(const size_t scanNumber);
	// read scans in a range
	void readScans(const size_t startScanNumber, const size_t endScanNumber);
	// read scans of scanNumbers
	void readScans(std::vector<size_t> &scanNumbers);
	void readAllScan();

	// for test
	void printFileInfo();
	void printScan(const Scan &mScan);
};