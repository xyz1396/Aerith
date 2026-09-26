#pragma once
#include "parameters.h"

// Calculations use each peptide's actual sourced atom composition.
class PeptideIsotopeCalculator
{
public:
    const string SIPatoms = "CHONPS";
    Composition pepComposition;
    AtomCounts pepAtomCounts{};
    // C,H,O,N,P,S Atom count of BYions
    std::vector<Composition> BionsCompositions, YionsCompositions;
    std::vector<AtomCounts> BionsAtomCounts;
    std::vector<std::array<int, 6>> YionsAtomCounts;
    std::vector<double> BionsBaseMasses;
    std::vector<double> YionsBaseMasses;
    int SIPatomIX;
    PeptideIsotopeCalculator();
    static bool changeAtomProbability(std::vector<double> &probs, char atom, const double pct);
    void changeAtomSIPabundance(const char SIPatom, const double pct);
    double calNetronMass(const string &pepSeq);
    void calPepAtomCounts(const string &pepSeq);
    void calBYionsAtomCounts(const string &pepSeq);
    // Synchronize the target element with the current parameters.
    void updateSIPelement();

    // for peptide base mass without isotope
    double calPrecursorBaseMass(const string &pepSeq);
    void calBYionBaseMasses(const string &pepSeq);
    double calPrecursorMass(const string &pepSeq);
};
