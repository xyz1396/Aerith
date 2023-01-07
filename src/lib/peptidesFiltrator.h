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
    string originalPepSeq;
    vector<string> proteinNames;
    // follows are in pairs
    vector<float> scores;
    vector<string> searchNames;
    vector<string> ftFileNames;
    vector<int> scanNumbers;
    vector<int> parentCharges;
};

class peptidesFiltrator
{
private:
public:
    float FDRthreshold;
    vector<string> tokens;
    // identified peptides with different charge states
    // we always don't select singly charged precursor in expriment
    // charge 2: 3 : >3 is nearly 100:10:1
    // score of precursor with high charge is low due to resolution limit
    // peptideMapCharge2 contains peptides with charge <=2
    unordered_map<string, peptideInfo> peptideMapCharge2;
    unordered_map<string, peptideInfo> peptideMapCharge3;
    unordered_map<string, peptideInfo> peptideMapChargeLargerThan3;
    int decoyCountCharge2, decoyCountCharge3, decoyCountChargeLargerThan3;
    int pepCountCharge2, pepCountCharge3, pepCountChargeLargerThan3;
    float scoreThresholdCharge2, scoreThresholdCharge3, scoreThresholdChargeLargerThan3;
    // for output filtered peptides
    unordered_map<string, peptideInfo> peptideMap;
    peptidesFiltrator(const vector<sipPSM> &sipPSMs, float mFDRthreshold);
    ~peptidesFiltrator();
    void splitString(const string mString);
    bool detectDecoy(const string &proteinName);
    void initFillPeptideMap(const string &pepSeq, const float score, const string &proName,
                            unordered_map<string, peptideInfo> &mMap);
    void sortPeptideBestScore(const unordered_map<string, peptideInfo> &mPepMap, vector<pair<float, bool>> &bestScoreDecoyPairs);
    void sortPeptideAllScore(const unordered_map<string, peptideInfo> &mPepMap, vector<pair<float, pair<string, bool>>> &scoreDecoyPairs);
    // return <decoyCount,pepCount,scoreThreshold>
    tuple<size_t, size_t, float> getDecoyCountScoreThreshold(vector<pair<float, bool>> &bestScoreDecoyPairs);
    tuple<size_t, size_t, float> getDecoyCountScoreThreshold(vector<pair<float, pair<string, bool>>> &scoreDecoyPairs);
    void filterPeptideMap();
    void fillPeptideMapExtraInfo();
};
