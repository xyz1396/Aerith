#pragma once
#include "sipFileReader.h"
#include <map>
#include <unordered_map>
#include <algorithm>
#include <unordered_set>

struct peptideInfo
{
    // for filtering
    bool isDecoy;
    float bestScore;
    peptideInfo(bool mIsDecoy, float mBestScore);
    // for output
    double parentMass;
    std::string originalPepSeq;
    std::vector<std::string> proteinNames;
    // follows are in std::pairs
    std::vector<float> scores;
    std::vector<std::string> searchNames;
    std::vector<std::string> ftFileNames;
    std::vector<int> scanNumbers;
    std::vector<int> parentCharges;
};

class peptidesFiltrator
{
private:
public:
    float FDRthreshold;
    std::vector<std::string> tokens;
    // identified peptides with different charge states
    // we always don't select singly charged precursor in expriment
    // charge 2: 3 : >3 is nearly 100:10:1
    // score of precursor with high charge is low due to resolution limit
    // peptideMapCharge2 contains peptides with charge <=2
    std::unordered_map<std::string, peptideInfo> peptideMapCharge2;
    std::unordered_map<std::string, peptideInfo> peptideMapCharge3;
    std::unordered_map<std::string, peptideInfo> peptideMapChargeLargerThan3;
    int decoyCountCharge2, decoyCountCharge3, decoyCountChargeLargerThan3;
    int pepCountCharge2, pepCountCharge3, pepCountChargeLargerThan3;
    float scoreThresholdCharge2, scoreThresholdCharge3, scoreThresholdChargeLargerThan3;
    // for output filtered peptides
    std::unordered_map<std::string, peptideInfo> peptideMap;
    peptidesFiltrator(const std::vector<sipPSM> &sipPSMs, float mFDRthreshold);
    ~peptidesFiltrator();
    void splitString(const std::string mString);
    bool detectDecoy(const std::string &proteinName);
    void initFillPeptideMap(const std::string &pepSeq, const float score, const std::string &proName,
                            std::unordered_map<std::string, peptideInfo> &mMap);
    void sortPeptideBestScore(const std::unordered_map<std::string, peptideInfo> &mPepMap, std::vector<std::pair<float, bool>> &bestScoreDecoyPairs);
    void sortPeptideAllScore(const std::unordered_map<std::string, peptideInfo> &mPepMap, std::vector<std::pair<float, std::pair<std::string, bool>>> &scoreDecoyPairs);
    // return <decoyCount,pepCount,scoreThreshold>
    std::tuple<size_t, size_t, float> getDecoyCountScoreThreshold(std::vector<std::pair<float, bool>> &bestScoreDecoyPairs);
    std::tuple<size_t, size_t, float> getDecoyCountScoreThreshold(std::vector<std::pair<float, std::pair<std::string, bool>>> &scoreDecoyPairs);
    void filterPeptideMap();
    void fillPeptideMapExtraInfo();
};
