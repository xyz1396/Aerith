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
        B,
        BisotopicPeak,
        Y,
        YisotopicPeak
    };
    PSMpeakAnnotator(const double tolerancePPM);
    ~PSMpeakAnnotator();
    void analyzePSM(const std::string &peptide, Scan *realScan,
                    // charges of BY ions in consideration
                    const std::vector<int> &charges, const double isoCenter = 0, 
                    const double isoWidth = 0, const bool calScores = false);
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
    // unit is ppm
    double tolerancePPM = 10;
    double matchedSpectraEntropyScore = 0, Score = 0, MVHscore = 0, WDPscore = 0, XcorrScore = 0;
    vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;

    // theoretical spectra info of same length vectors
    std::vector<double> expectedMZs, expectedIntensities;
    std::vector<int> expectedCharges;
    std::vector<ionKind> ionKinds;
    // residue position of BY ions
    std::vector<int> residuePositions, matchedIndices;
    // SIP abundance from BY ion isotopic envelope by binomial distribution estimation
    std::vector<double> SIPabundances;

    void generateTheoreticalSpectra(const std::string &peptide);
    void removePeaksInIsolationWindow(Scan *mScan, const double isoCenter, const double isoWidth);
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