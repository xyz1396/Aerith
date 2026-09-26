#include "isotopologue.h"
#include <Rcpp.h>
#include <stdexcept>

IsotopeDistribution::IsotopeDistribution()
{
}

IsotopeDistribution::IsotopeDistribution(vector<double> vItsMass, vector<double> vItsProb)
{
	vMass = vItsMass;
	vProb = vItsProb;
}

IsotopeDistribution::~IsotopeDistribution()
{
	// destructor
}

void IsotopeDistribution::print()
{
	Rcpp::Rcout << "Mass " << '\t' << "Inten" << endl;
	for (unsigned int i = 0; i < vMass.size(); i++)
	{
		Rcpp::Rcout << setprecision(8) << vMass[i] << '\t' << vProb[i] << endl;
	}
}

double IsotopeDistribution::getMostAbundantMass()
{
	double dMaxProb = 0;
	double dMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		if (dMaxProb < vProb[i])
		{
			dMaxProb = vProb[i];
			dMass = vMass[i];
		}
	}
	return dMass;
}

double IsotopeDistribution::getAverageMass()
{
	double dSumProb = 0;
	double dSumMass = 0;
	for (unsigned int i = 0; i < vMass.size(); ++i)
	{
		dSumProb = dSumProb + vProb[i];
		dSumMass = dSumMass + vProb[i] * vMass[i];
	}

	if (dSumProb <= 0)
		return 1.0;

	return (dSumMass / dSumProb);
}
void IsotopeDistribution::filterProbCutoff(double dProbCutoff)
{
	vector<double> vMassCopy = vMass;
	vector<double> vProbCopy = vProb;
	vMass.clear();
	vProb.clear();
	for (unsigned int i = 0; i < vProbCopy.size(); ++i)
	{
		if (vProbCopy[i] >= dProbCutoff)
		{
			vMass.push_back(vMassCopy[i]);
			vProb.push_back(vProbCopy[i]);
		}
	}
}
double IsotopeDistribution::getLowestMass()
{
	return *min_element(vMass.begin(), vMass.end());
}

Isotopologue::Isotopologue() : MassPrecision(0.01), ProbabilityCutoff(0.000000001)
{
}

Isotopologue::~Isotopologue()
{
	// destructor
}

void Isotopologue::setupIsotopologue(const map<string, Composition> &residues,
    const vector<IsotopeDistribution> &atoms)
{
    if (atoms.size() != sipros::ElementCount || !residues.count("Nterm") ||
        !residues.count("Cterm"))
        throw std::invalid_argument("Incomplete compiled chemistry parameters.");
    mResidueCompositions = residues;
    vAtomIsotopicDistribution = atoms;
    for (auto &atom : vAtomIsotopicDistribution)
    {
        if (atom.vMass.empty() || atom.vMass.size() != atom.vProb.size())
            throw std::invalid_argument("Invalid compiled isotope distribution.");
        double total = 0;
        for (size_t i = 0; i < atom.vMass.size(); ++i)
        {
            if (!std::isfinite(atom.vMass[i]) || !std::isfinite(atom.vProb[i]) || atom.vProb[i] < 0)
                throw std::invalid_argument("Invalid compiled isotope probability.");
            total += atom.vProb[i];
        }
        if (std::abs(total - 1) > 1e-8 || !CheckMass(atom.vMass, atom.vProb))
            throw std::invalid_argument("Invalid compiled isotope distribution.");
    }
    naturalAtomIsotopicDistribution = vAtomIsotopicDistribution;
    refreshResidueDistributions();
}

double Isotopologue::computeMostAbundantMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getMostAbundantMass();
}

double Isotopologue::computeAverageMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getAverageMass();
}

double Isotopologue::computeMonoisotopicMass(string sSequence)
{
	IsotopeDistribution tempIsotopeDistribution;
	if (!computeIsotopicDistribution(sSequence, tempIsotopeDistribution))
		return 0;
	else
		return tempIsotopeDistribution.getLowestMass();
}

bool Isotopologue::getSingleResidueMostAbundantMasses(vector<string> &vsResidues, vector<double> &vdMostAbundantMasses, double &dTerminusMassN,
													  double &dTerminusMassC)
{
	vsResidues.clear();
	vdMostAbundantMasses.clear();

	map<string, IsotopeDistribution>::iterator ResidueIter;
	string sCurrentResidue;
	IsotopeDistribution currentDistribution;
	double dCurrentMostAbundantMasses;

	// for single amino acid

	for (ResidueIter = vResidueIsotopicDistribution.begin(); ResidueIter != vResidueIsotopicDistribution.end(); ResidueIter++)
	{
		sCurrentResidue = ResidueIter->first;
		currentDistribution = ResidueIter->second;
		dCurrentMostAbundantMasses = currentDistribution.getMostAbundantMass();
		if (sCurrentResidue.size() == 1)
		{
			vsResidues.push_back(sCurrentResidue);
			vdMostAbundantMasses.push_back(dCurrentMostAbundantMasses);
			//			cout << sCurrentResidue << "  " << dCurrentMostAbundantMasses << endl;
		}
		else if (sCurrentResidue == "NTerm" || sCurrentResidue == "Nterm")
		{
			dTerminusMassN = dCurrentMostAbundantMasses;
		}
		else if (sCurrentResidue == "CTerm" || sCurrentResidue == "Cterm")
		{
			dTerminusMassC = dCurrentMostAbundantMasses;
		}
		else
		{
		Rcpp::Rcerr << "ERROR: Cannot recognize the configuration for residue " << sCurrentResidue << endl;
		}
	}

	unsigned int i;

	// bubble sort the list by mass
	unsigned int n = vsResidues.size();
	unsigned int pass;
	double dCurrentMass;
	for (pass = 1; pass < n; pass++)
	{ // count how many times
		// This next loop becomes shorter and shorter
		for (i = 0; i < n - pass; i++)
		{
			if (vdMostAbundantMasses.at(i) > vdMostAbundantMasses.at(i + 1))
			{
				// exchange
				dCurrentMass = vdMostAbundantMasses.at(i);
				sCurrentResidue = vsResidues.at(i);

				vdMostAbundantMasses.at(i) = vdMostAbundantMasses.at(i + 1);
				vsResidues.at(i) = vsResidues.at(i + 1);

				vdMostAbundantMasses.at(i + 1) = dCurrentMass;
				vsResidues.at(i + 1) = sCurrentResidue;
			}
		}
	}

	return true;
}

bool Isotopologue::computeIsotopicDistribution(string sequence, IsotopeDistribution &distribution)
{
    return computeIsotopicDistribution(peptideComposition(sequence), distribution);
}

bool Isotopologue::computeProductIon(string sequence,
    vector<vector<double>> &yMass, vector<vector<double>> &yProb,
    vector<vector<double>> &bMass, vector<vector<double>> &bProb)
{
    const auto fragments = fragmentCompositions(sequence);
    yMass.clear(); yProb.clear(); bMass.clear(); bProb.clear();
    for (const auto &composition : fragments.first)
    {
        IsotopeDistribution distribution;
        computeIsotopicDistribution(composition, distribution);
        bMass.push_back(distribution.vMass);
        bProb.push_back(distribution.vProb);
    }
    for (const auto &composition : fragments.second)
    {
        IsotopeDistribution distribution;
        computeIsotopicDistribution(composition, distribution);
        yMass.push_back(distribution.vMass);
        yProb.push_back(distribution.vProb);
    }
    return true;
}

void Isotopologue::refreshResidueDistributions()
{
    vResidueIsotopicDistribution.clear();
    for (const auto &residue : mResidueCompositions)
        computeIsotopicDistribution(residue.second,
            vResidueIsotopicDistribution[residue.first]);
}

bool Isotopologue::computeIsotopicDistribution(const Composition &composition,
    IsotopeDistribution &distribution)
{
    distribution = IsotopeDistribution({0.0}, {1.0});
    const auto naturalCounts = composition.naturalSourceTotal();
    for (size_t i = 0; i < sipros::ElementCount; ++i)
    {
        const int biological = composition[IsotopeSource::Biosynthetic][i];
        const auto &active = vAtomIsotopicDistribution.at(i);
        const auto &natural = naturalAtomIsotopicDistribution.at(i);
        if (active.vProb == natural.vProb)
        {
            if (biological + naturalCounts[i] != 0)
                distribution = sum(multiply(natural, biological + naturalCounts[i]), distribution);
        }
        else
        {
            if (biological != 0)
                distribution = sum(multiply(active, biological), distribution);
            if (naturalCounts[i] != 0)
                distribution = sum(multiply(natural, naturalCounts[i]), distribution);
        }
    }
    return true;
}

namespace {
void validateComposition(const Composition &composition)
{
    for (const auto &source : composition.atoms)
        for (int count : source)
            if (count < 0)
                throw std::invalid_argument("PTM removes more atoms than the residue contains.");
}

struct ParsedPeptide
{
    vector<Composition> residues;
    Composition nTermPtm, cTermPtm;
};

ParsedPeptide parsePeptide(const Isotopologue &iso, string sequence, bool fragments)
{
    ParsedPeptide result;
    string suffix;
    if (!sequence.empty() && sequence.front() == '[')
    {
        const auto end = sequence.find(']');
        if (end == string::npos)
            throw std::invalid_argument("Missing closing peptide bracket.");
        suffix = sequence.substr(end + 1);
        sequence = sequence.substr(1, end - 1);
    }
    auto lookup = [&](char symbol) -> Composition {
        if (fragments && symbol == '>') symbol = '1';
        if (fragments && symbol == '<') symbol = '2';
        const auto entry = iso.mResidueCompositions.find(string(1, symbol));
        if (entry == iso.mResidueCompositions.end())
            throw std::invalid_argument("Unknown residue or PTM: " + string(1, symbol));
        return entry->second;
    };
    auto isResidue = [](char symbol) {
        return std::isalpha(static_cast<unsigned char>(symbol)) != 0;
    };
    size_t i = 0;
    if (i < sequence.size() && !isResidue(sequence[i]))
        result.nTermPtm = lookup(sequence[i++]);
    while (i < sequence.size())
    {
        if (!isResidue(sequence[i]))
            throw std::invalid_argument("Only one PTM per residue is supported.");
        Composition residue = lookup(sequence[i++]);
        if (i < sequence.size() && !isResidue(sequence[i]))
            residue += lookup(sequence[i++]);
        validateComposition(residue);
        result.residues.push_back(residue);
    }
    if (result.residues.empty())
        throw std::invalid_argument("Peptide sequence must contain a residue.");
    if (!suffix.empty())
    {
        if (suffix.size() != 1 || isResidue(suffix[0]))
            throw std::invalid_argument("Invalid C-terminal PTM.");
        result.cTermPtm = lookup(suffix[0]);
    }
    validateComposition(result.nTermPtm);
    validateComposition(iso.mResidueCompositions.at("Nterm") +
        iso.mResidueCompositions.at("Cterm") + result.cTermPtm);
    return result;
}
}

Composition Isotopologue::peptideComposition(const string &sequence) const
{
    const auto peptide = parsePeptide(*this, sequence, false);
    Composition result = mResidueCompositions.at("Nterm") + mResidueCompositions.at("Cterm") +
        peptide.nTermPtm + peptide.cTermPtm;
    for (const auto &residue : peptide.residues) result += residue;
    return result;
}

std::pair<vector<Composition>, vector<Composition>>
Isotopologue::fragmentCompositions(const string &sequence) const
{
    const auto peptide = parsePeptide(*this, sequence, true);
    vector<Composition> bIons, yIons;
    Composition b = peptide.nTermPtm;
    Composition y = mResidueCompositions.at("Nterm") + mResidueCompositions.at("Cterm") + peptide.cTermPtm;
    for (size_t i = 0; i + 1 < peptide.residues.size(); ++i)
    {
        b += peptide.residues[i];
        y += peptide.residues[peptide.residues.size() - 1 - i];
        bIons.push_back(b); yIons.push_back(y);
    }
    return {bIons, yIons};
}

bool Isotopologue::computeAtomicComposition(string sequence, vector<int> &counts)
{
    const auto total = peptideComposition(sequence).total();
    counts.assign(total.begin(), total.end());
    return true;
}

double Isotopologue::estimateSIPAbundance(const Composition &composition,
    size_t element, double meanMassExcess) const
{
    const int count = composition[IsotopeSource::Biosynthetic].at(element);
    if (count <= 0) return 0.0;
    const auto &target = naturalAtomIsotopicDistribution.at(element);
    const size_t isotope = (element == 2 || element == 5) ? 2 : 1;
    const double delta = target.vMass.at(isotope) - target.vMass[0];
    double background = 0.0, minorShift = 0.0;
    for (size_t i = 0; i < composition[IsotopeSource::Biosynthetic].size(); ++i)
    {
        const auto &atom = naturalAtomIsotopicDistribution[i];
        double meanShift = 0.0;
        for (size_t j = 1; j < atom.vMass.size(); ++j)
        {
            const double shift = (atom.vMass[j] - atom.vMass[0]) * atom.vProb[j];
            meanShift += shift;
            if (i == element && j != isotope)
            {
                minorShift += shift;
            }
        }
        background += composition.naturalSourceTotal()[i] * meanShift;
        if (i != element) background += composition[IsotopeSource::Biosynthetic][i] * meanShift;
    }
    const double conditionalMinorShift = minorShift / (1.0 - target.vProb[isotope]);
    const double abundance = (meanMassExcess - background - count * conditionalMinorShift) /
        (count * (delta - conditionalMinorShift));
    return 100.0 * std::max(0.0, std::min(1.0, abundance));
}

IsotopeDistribution Isotopologue::sum(const IsotopeDistribution &distribution0, const IsotopeDistribution &distribution1)
{
	double ProbabilityCutoff_local = 0.000001;

	IsotopeDistribution sumDistribution;
	double currentMass;
	double currentProb;
	int iSizeDistribution0 = distribution0.vMass.size();
	int iSizeDistribution1 = distribution1.vMass.size();
	int iCount = 0;
	double dSum = 0;
	int newSize = iSizeDistribution0 + iSizeDistribution1 - 1;
	sumDistribution.vMass.reserve(newSize);
	sumDistribution.vProb.reserve(newSize);
#pragma omp simd
	for (int k = 0; k < newSize; k++)
	{
		double sumweight = 0, summass = 0;
		int start = k < (iSizeDistribution1 - 1) ? 0 : k - iSizeDistribution1 + 1; // max(0, k-f_n+1)
		int end = k < (iSizeDistribution0 - 1) ? k : iSizeDistribution0 - 1;	   // min(g_n - 1, k)
		iCount = 0;
		dSum = 0;
		for (int i = start; i <= end; i++)
		{
			double weight = distribution0.vProb[i] * distribution1.vProb[k - i];
			double mass = distribution0.vMass[i] + distribution1.vMass[k - i];
			sumweight += weight;
			summass += weight * mass;
			iCount += 1;
			dSum += mass;
		}
		if (sumweight == 0)
		{
			currentMass = dSum / ((double)iCount);
		}
		else
		{
			currentMass = summass / sumweight;
		}
		currentProb = sumweight;
		sumDistribution.vMass.push_back(currentMass);
		sumDistribution.vProb.push_back(currentProb);
	}

	// prune small probabilities
	vector<double>::iterator iteProb = sumDistribution.vProb.begin();
	vector<double>::iterator iteMass = sumDistribution.vMass.begin();
	while (iteProb != sumDistribution.vProb.end())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb++;
		iteMass++;
	}
	if (iteProb != sumDistribution.vProb.begin())
	{
		sumDistribution.vProb.erase(sumDistribution.vProb.begin(), iteProb);
		sumDistribution.vMass.erase(sumDistribution.vMass.begin(), iteMass);
	}

	iteProb = sumDistribution.vProb.end() - 1;
	iteMass = sumDistribution.vMass.end() - 1;
	while (iteProb != sumDistribution.vProb.begin())
	{
		if ((*iteProb) > ProbabilityCutoff_local)
		{
			break;
		}
		iteProb--;
		iteMass--;
	}
	if (iteProb != sumDistribution.vProb.end() - 1)
	{
		sumDistribution.vProb.erase(iteProb + 1, sumDistribution.vProb.end());
		sumDistribution.vMass.erase(iteMass + 1, sumDistribution.vMass.end());
	}

	// normalize the probability space to 1
	int iSizeSumDistribution;
	int i;
	double sumProb = 0;
	iSizeSumDistribution = sumDistribution.vMass.size();
	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumProb += sumDistribution.vProb[i];
	}

	if (sumProb <= 0)
	{
		throw std::runtime_error("Error: Sum of distribution is zero");
	}

	for (i = 0; i < iSizeSumDistribution; ++i)
	{
		sumDistribution.vProb[i] = sumDistribution.vProb[i] / sumProb;
	}

	return sumDistribution;
}

IsotopeDistribution Isotopologue::multiply(const IsotopeDistribution &distribution0, int count)
{
	if (count == 1)
		return distribution0;
	IsotopeDistribution productDistribution;
	productDistribution.vMass.push_back(0.0);
	productDistribution.vProb.push_back(1.0);

	if (count < 0)
	{
		IsotopeDistribution negativeDistribution = distribution0;
#pragma omp simd
		for (unsigned int n = 0; n < negativeDistribution.vMass.size(); n++)
		{
			negativeDistribution.vMass[n] = -negativeDistribution.vMass[n];
		}

		for (int i = 0; i < abs(count); ++i)
		{
			productDistribution = sum(productDistribution, negativeDistribution);
		}
	}
	else
	{
		for (int i = 0; i < count; ++i)
		{
			productDistribution = sum(productDistribution, distribution0);
		}
	}
	return productDistribution;
}

void Isotopologue::shiftMass(IsotopeDistribution &distribution0, double dMass)
{
	for (unsigned int i = 0; i < distribution0.vMass.size(); ++i)
	{
		distribution0.vMass[i] = distribution0.vMass[i] + dMass;
	}
}

bool Isotopologue::CheckMass(vector<double> &vdMass, vector<double> &vdNaturalCompositionTemp)
{
	double dDiff = 0;

	int iMassCount = vdMass.size();
	if (iMassCount < 2)
	{
		return true;
	}
	int pass, i;
	double dCurrentMass;
	double dCurrentComposition;
	// count how many times
	// This next loop becomes shorter and shorter
	for (pass = 1; pass < iMassCount; pass++)
	{
		for (i = 0; i < iMassCount - pass; i++)
		{
			if (vdMass.at(i) > vdMass.at(i + 1))
			{
				// exchange
				dCurrentMass = vdMass.at(i);
				dCurrentComposition = vdNaturalCompositionTemp.at(i);

				vdMass.at(i) = vdMass.at(i + 1);
				vdNaturalCompositionTemp.at(i) = vdNaturalCompositionTemp.at(i + 1);

				vdMass.at(i + 1) = dCurrentMass;
				vdNaturalCompositionTemp.at(i + 1) = dCurrentComposition;
			}
		}
	}

	for (i = 1; i < (int)vdMass.size(); ++i)
	{
		dDiff = round(vdMass.at(i) - vdMass.at(i - 1));
		if (dDiff != 1.0)
		{
			return false;
		}
	}
	return true;
}
