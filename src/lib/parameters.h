#pragma once
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>
#include <map>
#include <limits>
#include <numeric>
#include <sstream>
#include <algorithm>
#include "isotopologue.h"

// Runtime SIP settings are reset for each public calculation;
// no filesystem configuration is involved.
class AerithParameters
{
public:
    AerithParameters();
    static AerithParameters &current();
    static void reset();
    static char isotopeElement(const std::string &isotope);
    static void validateAbundance(double abundance);

    std::string chemistryProfile = "sipros5/source-aware-cam-tryptic-water/v1";
    std::string searchType = "SIP";
    int minPeptideLength = 7;
    int maxPeptideLength = 60;
    double fragmentToleranceDa = 0.01;
    double parentToleranceDa = 0.01;
    std::string cleavageAfter = "KR";
    std::string cleavageBefore = "ACDEFGHIJKLMNPQRSTVWY";
    int maxMissedCleavages = 2;
    std::vector<std::string> fixedPtms{"carbamidomethyl"};
    std::map<std::string, std::string> ptmSites;
    Isotopologue isotopologue;

    std::string sipElement = "C";
    double deductionMinValue = 0.005;
    double deductionFold = 4.0;
    double neutronMass = 1.003355;
    double deductionCoefficient = 0.0;
    void setDeductionCoefficient();
    double getDeductionCoefficient() const { return deductionCoefficient; }
    double getMassAccuracyFragmentIon() const { return fragmentToleranceDa; }
    int getMinPeptideLength() const { return minPeptideLength; }
    const std::string &getSearchType() const { return searchType; }
    double getProtonMass() const { return 1.00727646688; }
    double getNeutronMass() const { return neutronMass; }
    double getTerminusMassN() const { return 1.007825; }
    double getTerminusMassC() const { return 17.002740; }
    std::vector<std::pair<std::string, std::string>> getNeutralLossList() const
    { return {{">", "1"}, {"<", "2"}}; }
    double scoreError(double error) const
    { return std::erfc(std::abs(error) / (fragmentToleranceDa / 2 * std::sqrt(2.0))); }
};
