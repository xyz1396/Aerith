#pragma once
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;
using namespace std;

struct alignas(64) Scan
{
	int scanNumber;
	float retentionTime;
	double TIC;
	vector<double> mz;
	vector<double> mass;
	vector<double> intensity;

	// only MS2 scan has follows
	int precursorScanNumber;
	double precursorMz;
	int precursorCharge;

	// only orbitrap scan has follows
	vector<int> resolution;
	vector<float> baseLine;
	vector<float> signalToNoise;
	vector<int> charge;

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
	string ftFileName;
	ifstream ftFileStream;
	string currentLine;
	bool continueRead;
	bool hasPrecursor;
	bool hasCharge;
	// avoid empty file
	bool isEmpty = false;
	string instrument;
	string scanType;
	string scanFilter;
	Scan currentScan;
	vector<Scan> Scans;
	vector<string> tokens;
	ftFileReader(string file);
	~ftFileReader();
	void splitString(const string &mString);
	bool detectPrecursorAndCharge();
	bool hasNextScan();
	Scan readScanNumberRentionTime();
	Scan readScanNumberRentionTimePrecursor();
	void readPeakCharge();
	void readNextScan();
	Scan readOneScan(const int scanNumber);
	void readScans(const int scanCount);
	void readAllScan();

	// for test
	void printFileInfo();
	void printScan(const Scan &mScan);
};