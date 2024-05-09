#pragma once
#include "Spe2PepFileReader.h"
#include "ftFileReader.h"
#include "averagine.h"
#include "isotopicPeak.h"
#include <unordered_map>
#include <omp.h>

class PSMfeatureExtractor
{
public:
    PSMfeatureExtractor();
    Spe2PepFileReader mSpe2PepFileReader;
    averagine mAveragine = averagine();
    std::vector<Scan> FT1Scans;
    std::vector<Scan> FT2Scans;
    sipPSM *mSipPSM;
    std::unordered_map<size_t, Scan *> scanNumerFT1ScanMap;
    std::unordered_map<size_t, Scan *> scanNumerFT2ScanMap;
    size_t binarySearchPeak(const Scan *mScan, double Mz, int charge);
    // N isotopic peaks on each side to consider
    const static int NisotopicPeak = 20;
    std::vector<int> vertexIXs;
    // filter isotopic peaks by isotopic envolope shape
    void filterIsotopicPeaks(std::vector<isotopicPeak> &isotopicPeaks, const double calculatedPrecursorMZ);
    // return precursor scan number and isotopic peaks
    std::vector<isotopicPeak> findIsotopicPeaks(const size_t MS1ScanNumber,
                                                const int precursorCharge,
                                                const double observedPrecursorMass,
                                                const double calculatedPrecursorMass);
    double getSIPelementAbundanceFromMS1(const std::string &peptideSeq,
                                         const std::vector<isotopicPeak> &isotopicPeaks, const int precursorCharge);
    void loadFT1(const std::string &FTfileBasePath);
    void loadFT1FT2fileParallel(const std::string &FTfileBasePath);
    void extractFeaturesOfEachPSM();
    void extractPSMfeature(const std::string &Spe2PepFilePath, const int topN,
                           const std::string &ftFilepath);
    void extractPSMfeatureParallel(const std::string &Spe2PepFilePath, const int topN,
                                   const std::string &ftFilepath, const int threadNumber = 3);
    void writeTSV(const std::string &fileName);

private:
};
