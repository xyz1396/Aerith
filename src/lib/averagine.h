#pragma once
#include "proNovoConfig.h"

class averagine
{
private:
    /* data */
public:
    const string averagineResidue = "a";
    const string SIPatoms = "CHONPS";
    // C,H,O,N,P,S count of averagine residue
    vector<int> averagineAtomCounts;
    // C,H,O,N,P,S count of averagine peptides at different length
    vector<vector<int>> averaginePepAtomCountss;
    // Atom count difference bettween averagine and peptide
    vector<int> diffAtomCounts;
    IsotopeDistribution averagineSIPdistribution;
    vector<IsotopeDistribution> averaginePepSIPdistributions;
    int minPepLen, maxPepLen, pepLenRange, SIPatomIX;
    averagine(const int minPepLen, const int maxPepLen);
    ~averagine();
    void changeAtomSIPabundance(const char SIPatom, const double pct);
    void calAveraginePepAtomCounts();
    vector<int>* getAveraginePepAtomCounts(const int pepLen);
    void calAveraginePepSIPdistributions();
    IsotopeDistribution* getAveraginePepSIPdistribution(const int pepLen);
    void calDiffAtomCounts(const string &pepSeq);
    void calPrecusorIsotopeDistribution(const string &pepSeq, IsotopeDistribution &tempSIPdistribution);
};
