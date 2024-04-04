#pragma once
#include "Spe2PepFileReader.h"
#include "ftFileReader.h"
#include "averagine.h"

struct isotopicPeak
{
    double mass;
    int charge;
    double intensity;
};

class PSMfeatureExtractor
{
public:
    PSMfeatureExtractor();
    Spe2PepFileReader mSpe2PepFileReader;
    std::vector<Scan> FT1Scans;
    std::vector<Scan> FT2Scans;
    std::vector<isotopicPeak> findIsotopicPeaks(const double mz, const int charge,
                                                const int scanNumber,
                                                const std::vector<Scan> &ft1Scans,
                                                const std::vector<Scan> &ft2Scans);
    double getSIPelementAbundance(const std::string &peptideSeq,
                                  const std::vector<isotopicPeak> isotopicPeaks);
    void extractPSMfeature(const std::string &Spe2PepFilePath, const int topN,
                           const std::string &ftFilepath);

private:
};
