#include "PeptideIsotopeCalculator.h"

PeptideIsotopeCalculator::PeptideIsotopeCalculator()
{
    updateSIPelement();
}

void PeptideIsotopeCalculator::updateSIPelement()
{
    const auto index = SIPatoms.find(AerithParameters::current().sipElement);
    if (index == string::npos || index == 4)
        throw std::invalid_argument("Unsupported SIP element.");
    SIPatomIX = static_cast<int>(index);
}

void PeptideIsotopeCalculator::changeAtomSIPabundance(const char SIPatom, const double pct)
{
    const size_t atomPos = SIPatoms.find(SIPatom);
    if (atomPos == string::npos || atomPos == 4)
        throw std::invalid_argument("Unsupported SIP element.");
    AerithParameters::validateAbundance(pct);
    const int atomIndex = static_cast<int>(atomPos);
    SIPatomIX = atomIndex;
    AerithParameters::current().sipElement = SIPatom;
    updateSIPelement();
    // Each call starts from the natural profile, including minor O/S isotopes.
    AerithParameters::current().isotopologue.vAtomIsotopicDistribution[atomIndex].vProb =
        AerithParameters::current().isotopologue.naturalAtomIsotopicDistribution[atomIndex].vProb;
    changeAtomProbability(
        AerithParameters::current().isotopologue.vAtomIsotopicDistribution[atomIndex].vProb, SIPatom, pct);
    AerithParameters::current().setDeductionCoefficient();
}

bool PeptideIsotopeCalculator::changeAtomProbability(std::vector<double> &probs, const char atom, const double pct)
{
    AerithParameters::validateAbundance(pct);
    const size_t target = (atom == 'O' || atom == 'S') ? 2u : 1u;
    if (target >= probs.size())
        throw std::invalid_argument("SIP isotope is missing from the parameters.");
    double others = 0.0;
    for (size_t i = 0; i < probs.size(); ++i)
        if (i != target) others += probs[i];
    if (!(others > 0))
        throw std::invalid_argument("Natural isotope probabilities are invalid.");
    for (size_t i = 0; i < probs.size(); ++i)
        if (i != target) probs[i] *= (1.0 - pct) / others;
    probs[target] = pct;
    return true;
}

void PeptideIsotopeCalculator::calPepAtomCounts(const string &pepSeq)
{
    pepComposition = AerithParameters::current().isotopologue.peptideComposition(pepSeq);
    pepAtomCounts = pepComposition.total();
}

void PeptideIsotopeCalculator::calBYionsAtomCounts(const string &pepSeq)
{
    const auto fragments = AerithParameters::current().isotopologue.fragmentCompositions(pepSeq);
    BionsCompositions = fragments.first;
    YionsCompositions = fragments.second;
    BionsAtomCounts.clear(); YionsAtomCounts.clear();
    for (const auto &ion : BionsCompositions) BionsAtomCounts.push_back(ion.total());
    for (const auto &ion : YionsCompositions) YionsAtomCounts.push_back(ion.total());
}

namespace {
double baseMass(const Composition &composition)
{
    const auto counts = composition.total();
    const auto &atoms = AerithParameters::current().isotopologue.naturalAtomIsotopicDistribution;
    double mass = 0.0;
    for (size_t i = 0; i < counts.size(); ++i)
        mass += counts[i] * atoms[i].vMass[0];
    return mass;
}

// Expected isotope mass excess and nominal excess-neutron count for both pools.
std::pair<double, double> isotopeExcess(const Composition &composition)
{
    const auto &isotopologue = AerithParameters::current().isotopologue;
    double mass = 0.0, neutrons = 0.0;
    for (size_t i = 0; i < composition[IsotopeSource::Biosynthetic].size(); ++i)
    {
        for (int pool = 0; pool < 2; ++pool)
        {
            const int count = pool == 0 ? composition[IsotopeSource::Biosynthetic][i] : composition.naturalSourceTotal()[i];
            const auto &atom = pool == 0 ? isotopologue.vAtomIsotopicDistribution[i]
                : isotopologue.naturalAtomIsotopicDistribution[i];
            for (size_t j = 1; j < atom.vMass.size(); ++j)
            {
                const double delta = atom.vMass[j] - atom.vMass[0];
                mass += count * atom.vProb[j] * delta;
                neutrons += count * atom.vProb[j] * std::round(delta);
            }
        }
    }
    return {mass, neutrons};
}
}

double PeptideIsotopeCalculator::calNetronMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    const auto excess = isotopeExcess(pepComposition);
    if (!(excess.second > 0.0) || !std::isfinite(excess.first))
        throw std::runtime_error("Cannot estimate isotope spacing from this composition.");
    return excess.first / excess.second;
}

double PeptideIsotopeCalculator::calPrecursorBaseMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    return baseMass(pepComposition);
}

void PeptideIsotopeCalculator::calBYionBaseMasses(const string &pepSeq)
{
    calBYionsAtomCounts(pepSeq);
    BionsBaseMasses.clear();
    YionsBaseMasses.clear();
    for (const auto &ion : BionsCompositions) BionsBaseMasses.push_back(baseMass(ion));
    for (const auto &ion : YionsCompositions) YionsBaseMasses.push_back(baseMass(ion));
}

double PeptideIsotopeCalculator::calPrecursorMass(const string &pepSeq)
{
    calPepAtomCounts(pepSeq);
    const auto excess = isotopeExcess(pepComposition);
    if (!(excess.second > 0.0) || !std::isfinite(excess.first))
        throw std::runtime_error("Cannot estimate precursor mass from this composition.");
    IsotopeDistribution envelope;
    AerithParameters::current().isotopologue.computeIsotopicDistribution(pepComposition, envelope);
    if (envelope.vMass.empty() || envelope.vMass.size() != envelope.vProb.size())
        throw std::runtime_error("Cannot estimate precursor mass from an invalid isotope envelope.");
    const size_t apex = std::distance(envelope.vProb.begin(),
        std::max_element(envelope.vProb.begin(), envelope.vProb.end()));
    const double lightest = baseMass(pepComposition);
    const double nominalShift = std::round(envelope.vMass[apex] - lightest);
    return lightest + nominalShift * excess.first / excess.second;
}

