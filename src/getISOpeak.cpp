#include "lib/initSIP.h"
#include "lib/PeptideIsotopeCalculator.h"
#include <algorithm>
#include <Rcpp.h>
#include <utility>

using namespace Rcpp;

//' @title Precursor Peak Calculator
//' @description This function calculates the isotopic distribution of a given amino acid string and returns a DataFrame containing the mass and probability of each isotope.
//' @param AAstr A string representing the amino acid sequence.
//' @return A DataFrame with two columns: "Mass" containing the mass of each isotope and "Prob" containing the probability of each isotope.
//' @examples
//' a <- precursor_peak_calculator("PEPTIDE")
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator(String AAstr)
{ 
	AerithParameters::reset();
	IsotopeDistribution myIso;
	AerithParameters::current().isotopologue.computeIsotopicDistribution(AAstr, myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' Simple residue peak calculator of user defined isotopic distribution of one residue
//' @param residue residue name
//' @param Atom isotopes of "C13", "N15", "H2", "O18", "S34"
//' @param Prob its SIP abundance (0.0~1.0)
//' @return A DataFrame with two columns: "Mass" containing the mass of each isotope and "Prob" containing the probability of each isotope.
//' @examples
//' df <- residue_peak_calculator_DIY("A", "C13", 0.2)
//' @export
// [[Rcpp::export]]
DataFrame residue_peak_calculator_DIY(String residue, String Atom,
									  double Prob)
{
	AerithParameters::validateAbundance(Prob);
	// Reset the compiled parameter object.
	AerithParameters::reset();
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
	IsotopeDistribution myIso;
	auto residueIter = AerithParameters::current().isotopologue.vResidueIsotopicDistribution.find(residue);
	if (residueIter != AerithParameters::current().isotopologue.vResidueIsotopicDistribution.end())
		myIso = residueIter->second;
	else
		stop("Unknown residue or PTM: %s", residue.get_cstring());
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' @title Precursor Peak Calculator with User-Defined Isotopic Distribution
//' @description This function calculates the isotopic distribution of a given amino acid string with a user-defined isotopic distribution and returns a DataFrame containing the mass and probability of each isotope.
//' @param AAstr A string representing the amino acid sequence.
//' @param Atom A string representing the isotope ("C13", "N15", "H2", "O18", "S34").
//' @param Prob A double representing the abundance of the specified isotope (0.0 to 1.0).
//' @return A DataFrame with two columns: "Mass" containing the mass of each isotope and "Prob" containing the probability of each isotope.
//' @examples
//' # Example usage
//' df <- precursor_peak_calculator_DIY("PEPTIDE", "C13", 0.2)
//' df <- precursor_peak_calculator_DIY("PEPTIDE", "N15", 0.5)
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator_DIY(String AAstr, String Atom,
										double Prob)
{

	AerithParameters::validateAbundance(Prob);
	// Reset the compiled parameter object.
	AerithParameters::reset();
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
	IsotopeDistribution myIso;
	AerithParameters::current().isotopologue.computeIsotopicDistribution(AAstr, myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' Simple calculator of C H O N P S atom count of peptide
//' @param AAstrs a CharacterVector of peptides
//' @param pool Atom pool to count: `"total"` (default), `"sip"`, `"natural"`,
//' `"reagent"`, or `"solvent"`. Natural is the sum of reagent and solvent atoms.
//' @details Cysteine is IAA-blocked by default. An explicit `C/` annotation
//' represents the same fixed modification and is counted once.
//' @return a dataframe of C H O N P S atom count each row is for one peptide
//' @export
//' @examples
//' df <- calPepAtomCount(c("HKFL","ADCH"))
// [[Rcpp::export]]
DataFrame calPepAtomCount(StringVector AAstrs, String pool = "total")
{
    if (pool != "total" && pool != "sip" && pool != "natural" &&
        pool != "reagent" && pool != "solvent")
        stop("pool must be 'total', 'sip', 'natural', 'reagent', or 'solvent'");
	AerithParameters::reset();
	PeptideIsotopeCalculator peptideCalculator;
	vector<int> C(AAstrs.size(), 0);
	vector<int> H, O, N, P, S;
	H = O = N = P = S = C;
	for (int i = 0; i < AAstrs.size(); i++)
	{
		peptideCalculator.calPepAtomCounts(as<std::string>((AAstrs[i])));
        const auto counts = pool == "sip" ? peptideCalculator.pepComposition[IsotopeSource::Biosynthetic] :
            pool == "natural" ? peptideCalculator.pepComposition.naturalSourceTotal() :
            pool == "reagent" ? peptideCalculator.pepComposition[IsotopeSource::ReagentNatural] :
            pool == "solvent" ? peptideCalculator.pepComposition[IsotopeSource::DigestionSolvent] : peptideCalculator.pepAtomCounts;
		C[i] = counts[0];
		H[i] = counts[1];
		O[i] = counts[2];
		N[i] = counts[3];
		P[i] = counts[4];
		S[i] = counts[5];
	}
	DataFrame df = DataFrame::create(Named("C") = C, _("H") = H,
									 _("O") = O, _("N") = N,
									 _("P") = P, _("S") = S);
	return df;
}

//' Simple calculator of C H O N P S atom count and mass without isotope of B Y ions
//' @param AAstrs a CharacterVector of peptides
//' @return a list of data.frame of C H O N P S atom count and each data.frame is for one peptide
//' @export
//' @examples
//' peps <- calBYAtomCountAndBaseMass(c("HK~FL","AD!CH","~ILKMV"))
// [[Rcpp::export]]
List calBYAtomCountAndBaseMass(StringVector AAstrs)
{
	AerithParameters::reset();
	PeptideIsotopeCalculator peptideCalculator;
	List pepBYs(AAstrs.size());
	for (int i = 0; i < (int)AAstrs.size(); i++)
	{
		peptideCalculator.calBYionBaseMasses(as<std::string>((AAstrs[i])));
		int BYionsSize = peptideCalculator.BionsBaseMasses.size() + peptideCalculator.YionsBaseMasses.size();
		vector<int> C(BYionsSize, 0);
		vector<int> H, O, N, P, S;
		H = O = N = P = S = C;
		vector<string> BYkinds(BYionsSize);
		vector<double> BYbaseMasses(BYionsSize);
		for (size_t j = 0; j < peptideCalculator.BionsBaseMasses.size(); j++)
		{
			C[j] = peptideCalculator.BionsAtomCounts[j][0];
			H[j] = peptideCalculator.BionsAtomCounts[j][1];
			O[j] = peptideCalculator.BionsAtomCounts[j][2];
			N[j] = peptideCalculator.BionsAtomCounts[j][3];
			P[j] = peptideCalculator.BionsAtomCounts[j][4];
			S[j] = peptideCalculator.BionsAtomCounts[j][5];
			BYkinds[j] = "B" + to_string(j + 1);
			BYbaseMasses[j] = peptideCalculator.BionsBaseMasses[j];
		}
		int start = peptideCalculator.BionsBaseMasses.size();
		for (size_t j = 0; j < peptideCalculator.YionsBaseMasses.size(); j++)
		{
			C[j + start] = peptideCalculator.YionsAtomCounts[j][0];
			H[j + start] = peptideCalculator.YionsAtomCounts[j][1];
			O[j + start] = peptideCalculator.YionsAtomCounts[j][2];
			N[j + start] = peptideCalculator.YionsAtomCounts[j][3];
			P[j + start] = peptideCalculator.YionsAtomCounts[j][4];
			S[j + start] = peptideCalculator.YionsAtomCounts[j][5];
			BYkinds[j + start] = "Y" + to_string(j + 1);
			BYbaseMasses[j + start] = peptideCalculator.YionsBaseMasses[j];
		}
		DataFrame df = DataFrame::create(Named("C") = C, _("H") = H,
										 _("O") = O, _("N") = N,
										 _("P") = P, _("S") = S, _("Kind") = BYkinds,
										 _("BaseMass") = BYbaseMasses);
		pepBYs[i] = df;
	}
	pepBYs.names() = AAstrs;

	return pepBYs;
}

//' Estimate a representative precursor isotope mass
//' @details Uses the nominal shift of the isotope envelope's most abundant
//' peak and the mean isotope spacing across all atom sources.
//' This is an estimate of the modal peak mass.
//' @param AAstrs a CharacterVector of peptides
//' @param Atom a Character of "C13", "H2", "O18", "N15", or "S34"
//' @param Probs a NumericVector with the same length of AAstr for SIP abundances
//' @return a vector of peptide precursor masses
//' @export
//' @examples
//' masses <- calPepPrecursorMass(c("HKFL", "ADCH"), "C13", c(0.2, 0.3))
// [[Rcpp::export]]
NumericVector calPepPrecursorMass(StringVector AAstrs, String Atom, NumericVector Probs)
{
    if (Probs.size() != AAstrs.size())
        stop("AAstrs and Probs must have equal lengths.");
    const char atom = AerithParameters::isotopeElement(Atom.get_cstring());
    for (double probability : Probs) AerithParameters::validateAbundance(probability);
    AerithParameters::reset();
    PeptideIsotopeCalculator calculator;
    NumericVector v(AAstrs.size());
    for (int i = 0; i < AAstrs.size(); ++i)
    {
        calculator.changeAtomSIPabundance(atom, Probs[i]);
        v[i] = calculator.calPrecursorMass(as<std::string>(AAstrs[i]));
    }
	return v;
}

//' Simple calculator neutron mass by average delta mass of each isotope
//' @param AAstrs a CharacterVector of peptides
//' @param Atom a Character of "C13", "H2", "O18", "N15", or "S34"
//' @param Probs a NumericVector with the same length of AAstr for SIP abundances
//' @return a vector of peptide neutron masses
//' @export
//' @examples
//' masses <- calPepNeutronMass(c("HKFL", "ADCH"), "C13", c(0.2, 0.3))
// [[Rcpp::export]]
NumericVector calPepNeutronMass(StringVector AAstrs, String Atom, NumericVector Probs)
{
    if (Probs.size() != AAstrs.size())
        stop("AAstrs and Probs must have equal lengths.");
    const char atom = AerithParameters::isotopeElement(Atom.get_cstring());
    for (double probability : Probs) AerithParameters::validateAbundance(probability);
    AerithParameters::reset();
    PeptideIsotopeCalculator calculator;
    NumericVector v(AAstrs.size());
    for (int i = 0; i < AAstrs.size(); ++i)
    {
        calculator.changeAtomSIPabundance(atom, Probs[i]);
        v[i] = calculator.calNetronMass(as<std::string>(AAstrs[i]));
    }
	return v;
}

//' @title BY Ion Peak Calculator with User-Defined Isotopic Distribution
//' @description This function calculates the isotopic distribution of B and Y ions for a given amino acid string with a user-defined isotopic distribution and returns a DataFrame containing the mass, probability, and type of each ion.
//' @param AAstr A string representing the amino acid sequence.
//' @param Atom A string representing the isotope ("C13", "N15", "H2", "O18", "S34").
//' @param Prob A double representing the abundance of the specified isotope (0.0 to 1.0).
//' @return A DataFrame with three columns: "Mass" containing the mass of each ion, "Prob" containing the probability of each ion, and "Kind" indicating whether the ion is a B or Y ion.
//' @examples
//' # Example usage
//' df <- BYion_peak_calculator_DIY("PEPTIDE", "C13", 0.2)
//' df <- BYion_peak_calculator_DIY("PEPTIDE", "N15", 0.5)
//' @export
// [[Rcpp::export]]
DataFrame BYion_peak_calculator_DIY(String AAstr, String Atom,
									double Prob)
{
	AerithParameters::validateAbundance(Prob);
	// Reset the compiled parameter object.
	AerithParameters::reset();
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
    string AAstr_str = AAstr.get_cstring();
    // AA string format is [AAKRCI] for example
	AAstr_str = "[" + AAstr_str + "]";
	vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;
	AerithParameters::current().isotopologue.computeProductIon(AAstr_str, vvdYionMass,
														vvdYionProb, vvdBionMass, vvdBionProb);
	vector<double> masses, probs;
	vector<string> kinds;
	for (size_t i = 0; i < vvdYionMass.size(); i++)
	{
		for (size_t j = 0; j < vvdYionMass[i].size(); j++)
		{
			masses.push_back(vvdYionMass[i][j]);
			probs.push_back(vvdYionProb[i][j]);
			kinds.push_back("Y" + to_string(i + 1));
		}
	}
	for (size_t i = 0; i < vvdBionMass.size(); i++)
	{
		for (size_t j = 0; j < vvdBionMass[i].size(); j++)
		{
			masses.push_back(vvdBionMass[i][j]);
			probs.push_back(vvdBionProb[i][j]);
			kinds.push_back("B" + to_string(i + 1));
		}
	}
	DataFrame df =
		DataFrame::create(Named("Mass") = std::move(masses),
						  _["Prob"] = std::move(probs), _["Kind"] = std::move(kinds));
	return df;
}

//' Inspect Aerith's compiled parameters
//' @description Returns a snapshot of the compiled peptide chemistry and SIP
//' scoring defaults. Calculations initialize these
//' parameters directly and do not read configuration files.
//' @return A list containing the chemistry profile, fixed PTMs, source-specific
//' residue formulas, isotope distributions, and scoring defaults.
//' @examples
//' parameters <- getAerithParameters()
//' parameters$chemistryProfile
//' parameters$residues$reagent["C", ]
//' @export
// [[Rcpp::export]]
List getAerithParameters()
{
    const AerithParameters parameters;
    const auto &iso = parameters.isotopologue;
    List atoms(6), residues(3);
    CharacterVector elements = CharacterVector::create("C", "H", "O", "N", "P", "S");
    atoms.attr("names") = elements;
    for (int i = 0; i < 6; ++i)
        atoms[i] = DataFrame::create(_["Mass"] = iso.naturalAtomIsotopicDistribution[i].vMass,
            _["Prob"] = iso.naturalAtomIsotopicDistribution[i].vProb);
    for (int source = 0; source < 3; ++source)
    {
        IntegerMatrix counts(iso.mResidueCompositions.size(), 6);
        CharacterVector names(iso.mResidueCompositions.size());
        int row = 0;
        for (const auto &entry : iso.mResidueCompositions)
        {
            names[row] = entry.first;
            for (int element = 0; element < 6; ++element)
                counts(row, element) = entry.second.atoms[source][element];
            ++row;
        }
        counts.attr("dimnames") = List::create(names, elements);
        residues[source] = counts;
    }
    residues.attr("names") = CharacterVector::create("sip", "reagent", "solvent");
    return List::create(_["chemistryProfile"] = parameters.chemistryProfile,
        _["fixedPtms"] = parameters.fixedPtms, _["residues"] = residues,
        _["isotopes"] = atoms, _["searchType"] = parameters.searchType,
        _["fragmentToleranceDa"] = parameters.fragmentToleranceDa,
        _["parentToleranceDa"] = parameters.parentToleranceDa,
        _["deductionMinValue"] = parameters.deductionMinValue,
        _["deductionFold"] = parameters.deductionFold);
}
