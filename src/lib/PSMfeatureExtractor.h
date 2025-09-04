#pragma once
#include "Spe2PepFileReader.h"
#include "ftFileReader.h"
#include "averagine.h"
#include "isotopicPeak.h"
#include <unordered_map>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <array>

class PSMfeatureExtractor
{
public:
    PSMfeatureExtractor();
    Spe2PepFileReader mSpe2PepFileReader;
    averagine mAveragine = averagine();
    std::vector<Scan> FT1Scans;
    std::vector<Scan> FT2Scans;
    sipPSM *mSipPSM;
    std::unordered_map<size_t, Scan *> scanNumerFT1ScanMap, scanNumerFT2ScanMap;
    // N isotopic peaks on each side to consider
    const static int NisotopicPeak = 20;
    std::vector<int> vertexIXs;
    std::array<char, 2> cleavageSites = {'K', 'R'};

    void loadFT1(const std::string &FTfileBasePath);
    void loadFT2(const std::string &FTfileBasePath);
    void loadFT1FT2fileParallel(const std::string &FTfileBasePath);
    size_t binarySearchPeak(const Scan *mScan, double Mz, int charge);
    // filter isotopic peaks by isotopic envolope shape
    void filterIsotopicPeaks(std::vector<isotopicPeak> &isotopicPeaks, const double calculatedPrecursorMZ);
    void filterIsotopicPeaksTopN(std::vector<isotopicPeak> &isotopicPeaks, const double observedPrecursorMZ,
                                 const size_t topN);
    // return precursor scan number and isotopic peaks
    std::vector<isotopicPeak> findIsotopicPeaks(int &MS1ScanNumber,
                                                const int precursorCharge,
                                                const double observedPrecursorMass);
    double getSIPelementAbundanceFromMS1(const std::string &peptideSeq,
                                         const std::vector<isotopicPeak> &isotopicPeaks, const int precursorCharge);
    std::pair<int, int> getSeqLengthAndMissCleavageSiteNumber(const std::string &peptideSeq);
    int getPTMnumber(const std::string &peptideSeq);
    std::pair<int, double> getMassWindowShiftAndError(const double observedPrecursorMass,
                                                      const double calculatedPrecursorMass);
    double getMS2IsotopicAbundance(const std::string &searchName);
    void extractFeaturesOfEachPSM();
    void extractPSMfeature(const std::string &Spe2PepFilePath, const int topN,
                           const std::string &ftFilepath);
    void extractPSMfeatureParallel(const std::string &ftFilepath, const int threadNumber = 3);
    void extractPSMfeatureParallel(const std::string &Spe2PepFilePath, const int topN,
                                   const std::string &ftFilepath, const int threadNumber = 3);
    void extractPSMfeatureParallel(const std::string &targetPath, const std::string &decoyPath, const int topN,
                                   const std::string &ftFilepath, const int threadNumber = 3);
    void writeTSV(const std::string &fileName);
    void writePecorlatorPin(const std::string &fileName, bool doProteinInference);

private:
};
