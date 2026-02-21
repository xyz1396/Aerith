#pragma once
#include <vector>
#include "ftFileReader.h"
#include "averagine.h"
#include "ms2scan.h"
#include <unordered_map>

using SparseBins = std::unordered_map<int, double>;

class PSMpeakAnnotator
{
public:
    enum ionKind
    {
        // for fragment ions
        B,
        BisotopicPeak,
        Y,
        YisotopicPeak,
        // for precursor ions
        P,
        PisotopicPeak
    };
    struct ionMatch
    {
        double expectedMZ;
        double observedMZ;
        double expectedIntensity;
        double observedIntensity;
        int expectedCharge;
        int observedCharge;
        ionKind kind;
        int residuePosition;
        size_t matchedIndex;
    };
    PSMpeakAnnotator(const double tolerancePPM);
    ~PSMpeakAnnotator();
    void analyzePSM(const std::string &peptide, Scan *realScan,
                    // charges of BY ions in consideration
                    const std::vector<int> &charges, const double isoCenter = 0, 
                    const double isoWidth = 0, const bool calScores = false);
    // analyze precursor isotopic envelope
    // realScan should be MS1 scan 
    void analyzePrecursor(const std::string &peptide, Scan *realScan,
                    const int charge, const double isoCenter = 0, 
                    const double isoWidth = 0, const bool calScores = false
    );
    double getScore();
    double getMatchedSpectraEntropyScore();
    std::vector<int> getMatchedIndices();
    std::vector<ionKind> getIonKinds();
    // get SIP abundance from BY ion isotopic envelope by binomial distribution estimation
    std::vector<double> getSIPabundances();
    std::vector<int> getResiduePositions();
    std::vector<double> getExpectedMZs();
    std::vector<double> getExpectedIntensities();
    std::vector<int> getExpectedCharges();
    double getMVHscore();
    double getWDPscore();
    double getXcorrScore();

private:
    std::string peptide;
    averagine mAveragine = averagine();
    Scan *realScan;
    double precursorBaseMass = 0.0;
    int precursorSIPatomCount = 0;
    // unit is ppm
    double tolerancePPM = 10;
    double matchedSpectraEntropyScore = 0, Score = 0, MVHscore = 0, WDPscore = 0, XcorrScore = 0;
    vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;

    // theoretical spectra info of same length vectors
    std::vector<double> expectedMZs, expectedIntensities;
    std::vector<int> expectedCharges;
    std::vector<ionKind> ionKinds;
    // residue position of BY/P ions
    std::vector<int> residuePositions, matchedIndices;
    // SIP abundance from BY ion isotopic envelope by binomial distribution estimation
    std::vector<double> SIPabundances;
    // matched B/Y ion isotopic peaks of envolope for scoring
    // key: ion kind and residue position, value: matched indices, first ion is the highest intensity ion
    std::map<std::pair<ionKind, int>, std::vector<ionMatch>> matchedIonIsotopicEnvelopes;

    void generateTheoreticalSpectra(const std::string &peptide);
    void removePeaksInIsolationWindow(Scan *mScan, const double isoCenter, const double isoWidth);
    void keepPeaksInIsolationWindow(Scan *mScan, const double isoCenter, const double isoWidth);
    size_t binarySearchPeak(const Scan *mScan, double Mz);
    void findIsotopicPeaks(const std::vector<double> &ionMasses,
                           const std::vector<double> &ionIntensities,
                           const Scan *mScan,
                           const int residuePosition,
                           const int charge,
                           const PSMpeakAnnotator::ionKind BYkind);
    // calculate SIP abundance from BY ion isotopic envelope by binomial distribution estimation
    double calSIPabundancesOfBYion(const double baseMass,
                                   const std::vector<int> &matchedIXs,
                                   const Scan *mScan, const int SIPelementCount,
                                   const int charge);
    void matchIsotopicEnvelopes(Scan *mRealScan, const int charge);
    void calMatchedSpectraEntropyScore();
    void scorePSM();
    double logBinom(int n, int k);
    double hypergeomProbability(int n_theo, int n_obs, int j, int N);
    void calMVHscore();
    void calWDPscore();
    // select top N peaks of each isotopic envelope
    void selectTopPeaks(std::vector<double> &theo_mz,
        std::vector<double> &theo_intensity, int topN, double intensityThreshold);
    SparseBins sparseBinSpectrum(const std::vector<double> &mz,
        const vector<double> &intensity, double binWidth,
        double minMz, double maxMz);
    double sparseCorrelationAtShift(const SparseBins &theoBins,
        const SparseBins &obsBins, int shift);
    void calXcorrScore();
};