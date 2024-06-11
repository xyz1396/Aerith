#pragma once
#include "proNovoConfig.h"
#include <numeric>
#include <functional>

class averagine
{
private:
    /* data */
public:
    double *C12Mass, *C13Mass, *C13Abundance, *H1Mass, *H2Mass, *H2Abundance;
    double *O16Mass, *O17Mass, *O17Abundance, *O18Mass, *O18Abundance;
    double *N14Mass, *N15Mass, *N15Abundance, *PfakeLowMass, *PfakeMass, *PfakeAbundance;
    double *S32Mass, *S33Mass, *S33Abundance, *S34Mass, *S34Abundance, *S36Mass, *S36Abundance;
    const string averagineResidue = "a";
    const string SIPatoms = "CHONPS";
    // C,H,O,N,P,S count of averagine residue
    vector<int> averagineAtomCounts;
    // C,H,O,N,P,S count of averagine peptides at different length
    vector<vector<int>> averaginePepAtomCountss;
    // C,H,O,N,P,S Atom count of peptides
    std::array<int, 6> pepAtomCounts = {0};
    // Atom count difference bettween averagine and peptide
    vector<int> diffAtomCounts;
    IsotopeDistribution averagineSIPdistribution;
    vector<IsotopeDistribution> averaginePepSIPdistributions;
    int minPepLen, maxPepLen, pepLenRange, SIPatomIX;
    averagine(const int minPepLen, const int maxPepLen);
    averagine();
    ~averagine();
    void changeAtomSIPabundance(const char SIPatom, const double pct);
    double weighted_mean(const std::vector<double> &values, const std::vector<double> &weights);
    double calNetronMass(const string &pepSeq);
    void calAveraginePepAtomCounts();
    vector<int> *getAveraginePepAtomCounts(const int pepLen);
    void calAveraginePepSIPdistributions();
    IsotopeDistribution *getAveraginePepSIPdistribution(const int pepLen);
    void calPepAtomCounts(const string &pepSeq);
    // init it in init function and changeAtomSIPabundance
    void adjustEstimatePrecursorMassbyNP(); 
    std::function<double(double, int, double)> estimatePrecursorMassbyNP;
    // for peptide base mass without isotope
    double calPrecursorBaseMass(const string &pepSeq);
    double calPrecursorMass(const string &pepSeq);
    void calDiffAtomCounts(const string &pepSeq);
    void calPrecursorIsotopeDistribution(const string &pepSeq, IsotopeDistribution &tempSIPdistribution);
};
