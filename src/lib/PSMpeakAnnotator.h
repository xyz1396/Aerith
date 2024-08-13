#pragma once
#include <vector>
#include "ftFileReader.h"
#include "averagine.h"
#include "ms2scan.h"

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
                    const std::vector<int> &charges, const double isoCenter = 0, const double isoWidth = 0);
    double getScore();
    double getMatchedSpectraEntropy();
    std::vector<int> getMatchedIndices();
    std::vector<ionKind> getIonKinds();
    std::vector<int> getResiduePositions();
    std::vector<double> getExpectedMZs();
    std::vector<double> getExpectedIntensities();
    std::vector<int> getExpectedCharges();

private:
    std::string peptide;
    Scan *realScan;
    // unit is ppm
    double tolerancePPM = 10;
    double matchedSpectraEntropy = 0;
    double Score = 0;
    vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;

    // theoretical spectra info of same length vectors
    std::vector<double> expectedMZs, expectedIntensities;
    std::vector<int> expectedCharges;
    std::vector<ionKind> ionKinds;
    // residue position of BY ions
    std::vector<int> residuePositions, matchedIndices;

    void generateTheoreticalSpectra(const std::string &peptide);
    void removePeaksInIsolationWindow(Scan *mScan, const double isoCenter, const double isoWidth);
    size_t binarySearchPeak(const Scan *mScan, double Mz);
    void findIsotopicPeaks(const std::vector<double> &ionMasses,
                           const std::vector<double> &ionIntensities,
                           const Scan *mScan,
                           const int residuePosition,
                           const int charge,
                           const PSMpeakAnnotator::ionKind BYkind);
    void matchIsotopicEnvelopes(Scan *mRealScan, const int charge);
    void calMatchedSpectraEntropy();
    double scorePSM();
};