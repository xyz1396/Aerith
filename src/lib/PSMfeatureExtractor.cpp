#include "PSMfeatureExtractor.h"

PSMfeatureExtractor::PSMfeatureExtractor()
{
}

std::vector<isotopicPeak> PSMfeatureExtractor::findIsotopicPeaks(const double mz, const int charge,
                                                                 const int scanNumber,
                                                                 const std::vector<Scan> &ft1Scans,
                                                                 const std::vector<Scan> &ft2Scans)
{
}

double PSMfeatureExtractor::getSIPelementAbundance(const std::string &peptideSeq,
                                                   const std::vector<isotopicPeak> isotopicPeaks)
{
}

void PSMfeatureExtractor::extractPSMfeature(const std::string &Spe2PepFilePath, const int topN,
                                            const std::string &ftFilepath)
{
    mSpe2PepFileReader.readSpe2PepFilesScansTopPSMsFromEachFT2Parallel(Spe2PepFilePath, topN);
    std::string FTfileBaseName;
    for (size_t i = 0; i < mSpe2PepFileReader.FT2s.size(); i++)
    {
        ftFileReader FT2FileReader(ftFilepath + "/" + mSpe2PepFileReader.FT2s[i]);
        FT2FileReader.readAllScan();
        FT2Scans = FT2FileReader.Scans;
        FTfileBaseName = mSpe2PepFileReader.FT2s[i]
                             .substr(0, mSpe2PepFileReader.FT2s[i].size() - 4);
        ftFileReader FT1FileReader(ftFilepath + "/" + FTfileBaseName + ".FT1");
        FT2FileReader.readAllScan();
        FT2Scans = FT2FileReader.Scans;
    }
}