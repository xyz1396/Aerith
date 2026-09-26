#include "parameters.h"
#include <algorithm>

AerithParameters &AerithParameters::current()
{
    static AerithParameters parameters;
    return parameters;
}

void AerithParameters::reset()
{
    // Rebuild the object explicitly; a failed initialization is an error.
    current() = AerithParameters();
}

AerithParameters::AerithParameters()
{
    std::vector<IsotopeDistribution> atoms;
    std::map<std::string, Composition> residues;
    atoms = {
        IsotopeDistribution({12.000000, 13.003355},
                            {0.9893, 0.0107}),
        IsotopeDistribution({1.007825, 2.014102},
                            {0.999885, 0.000115}),
        IsotopeDistribution({15.994915, 16.999132, 17.999160},
                            {0.99757, 0.00038, 0.00205}),
        IsotopeDistribution({14.003074, 15.000109},
                            {0.99632, 0.00368}),
        IsotopeDistribution({30.973762}, {1.0}),
        IsotopeDistribution(
            {31.972071, 32.971459, 33.967867, 34.967867, 35.967081},
            {0.9493, 0.0076, 0.0429, 0.0, 0.0002})};

    const auto biosynthetic = [](const sipros::AtomCounts &counts)
    {
        return sipros::compositionFrom(
            sipros::IsotopeSource::Biosynthetic, counts);
    };
    const auto reagentNatural = [](const sipros::AtomCounts &counts)
    {
        return sipros::compositionFrom(
            sipros::IsotopeSource::ReagentNatural, counts);
    };
    const auto digestionSolvent = [](const sipros::AtomCounts &counts)
    {
        return sipros::compositionFrom(
            sipros::IsotopeSource::DigestionSolvent, counts);
    };

    residues = {
        {"Nterm", digestionSolvent({0, 1, 0, 0, 0, 0})},
        {"Cterm", digestionSolvent({0, 1, 1, 0, 0, 0})},
        {"J", biosynthetic({6, 11, 1, 1, 0, 0})},
        {"I", biosynthetic({6, 11, 1, 1, 0, 0})},
        {"L", biosynthetic({6, 11, 1, 1, 0, 0})},
        {"A", biosynthetic({3, 5, 1, 1, 0, 0})},
        {"S", biosynthetic({3, 5, 2, 1, 0, 0})},
        {"G", biosynthetic({2, 3, 1, 1, 0, 0})},
        {"V", biosynthetic({5, 9, 1, 1, 0, 0})},
        {"E", biosynthetic({5, 7, 3, 1, 0, 0})},
        {"K", biosynthetic({6, 12, 1, 2, 0, 0})},
        {"T", biosynthetic({4, 7, 2, 1, 0, 0})},
        {"D", biosynthetic({4, 5, 3, 1, 0, 0})},
        {"R", biosynthetic({6, 12, 1, 4, 0, 0})},
        {"P", biosynthetic({5, 7, 1, 1, 0, 0})},
        {"N", biosynthetic({4, 6, 2, 2, 0, 0})},
        {"F", biosynthetic({9, 9, 1, 1, 0, 0})},
        {"Q", biosynthetic({5, 8, 2, 2, 0, 0})},
        {"Y", biosynthetic({9, 9, 2, 1, 0, 0})},
        {"M", biosynthetic({5, 9, 1, 1, 0, 1})},
        {"H", biosynthetic({6, 7, 1, 3, 0, 0})},
        {"C", biosynthetic({3, 5, 1, 1, 0, 1})},
        {"W", biosynthetic({11, 10, 1, 2, 0, 0})},
        // Oxidation adds oxygen from an exogenous reagent/air pool.
        {"~", reagentNatural({0, 0, 1, 0, 0, 0})},
        // Deamidation removes peptide H/N and incorporates natural O.
        {"!", biosynthetic({0, -1, 0, -1, 0, 0}) +
                  reagentNatural({0, 0, 1, 0, 0, 0})},
        {"@", biosynthetic({0, 1, 3, 0, 1, 0})},
        {">", biosynthetic({0, 1, 3, 0, 1, 0})},
        {"<", biosynthetic({0, 1, 3, 0, 1, 0})},
        // Chemistry-only fragment replacements for phospho neutral loss.
        {"1", biosynthetic({0, 0, 0, 0, 0, 0})},
        {"2", biosynthetic({0, -2, -1, 0, 0, 0})},
        {"%", biosynthetic({2, 2, 1, 0, 0, 0})},
        {"^", biosynthetic({1, 2, 0, 0, 0, 0})},
        {"&", biosynthetic({2, 4, 0, 0, 0, 0})},
        {"*", biosynthetic({3, 6, 0, 0, 0, 0})},
        {")", biosynthetic({0, -1, 0, 0, 0, 0}) +
              reagentNatural({0, 0, 2, 1, 0, 0})},
        {"$", biosynthetic({1, 2, 0, 0, 0, 1})}};

    const auto cam = reagentNatural({2, 3, 1, 1, 0, 0});
    residues.at("C") += cam;
    residues["/"] = {};
    residues["("] = biosynthetic({0, -1, 0, 0, 0, 0}) +
        reagentNatural({0, 0, 1, 1, 0, 0}) - cam;
    ptmSites = {{"~", "M"}, {"!", "NQ"}, {"@", "STYHD"},
        {">", "STYHD"}, {"<", "ST"}, {"%", "K"}, {"^", "KRED"},
        {"&", "KR"}, {"*", "K"}, {"(", "C"}, {")", "Y"},
        {"/", "C"}, {"$", "D"}};
    isotopologue.setupIsotopologue(residues, atoms);
    setDeductionCoefficient();
}

char AerithParameters::isotopeElement(const std::string &isotope)
{
    if (isotope == "C13") return 'C';
    if (isotope == "H2") return 'H';
    if (isotope == "O18") return 'O';
    if (isotope == "N15") return 'N';
    if (isotope == "S34") return 'S';
    throw std::invalid_argument("Unsupported SIP isotope; use C13, H2, O18, N15 or S34.");
}

void AerithParameters::validateAbundance(double abundance)
{
    if (!std::isfinite(abundance) || abundance < 0 || abundance > 1)
        throw std::invalid_argument("SIP abundance must be finite and within [0, 1].");
}

void AerithParameters::setDeductionCoefficient()
{
    const auto element = std::string("CHONPS").find(sipElement);
    if (sipElement.size() != 1 || element == std::string::npos || element == 4)
        throw std::invalid_argument("Unsupported SIP element.");
    const size_t isotope = element == 2 || element == 5 ? 2 : 1;
    const auto &atom = isotopologue.vAtomIsotopicDistribution.at(element);
    neutronMass = (atom.vMass.at(isotope) - atom.vMass.at(0)) /
        ((element == 2 || element == 5) ? 2.0 : 1.0);
    deductionCoefficient = -(deductionMinValue + deductionFold *
        std::pow(atom.vProb.at(isotope) - 0.5, 8));
}
