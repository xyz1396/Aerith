#include "lib/initSIP.h"
#include "lib/averagine.h"
#include <Rcpp.h>

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
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	IsotopeDistribution myIso;
	ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr, myIso);
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
	if (Prob < 0 || Prob > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
	IsotopeDistribution myIso;
	auto residueIter = ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.find(residue);
	if (residueIter != ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution.end())
		myIso = residueIter->second;
	else
		Rcout << "Rediue not found!" << endl;
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
	
	if (Prob < 0 || Prob > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
	IsotopeDistribution myIso;
	ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr, myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' Simple calculator of C H O N P S atom count of peptide
//' @param AAstrs a CharacterVector of peptides
//' @return a dataframe of C H O N P S atom count each row is for one peptide
//' @export
//' @examples
//' df <- calPepAtomCount(c("HKFL","ADCH"))
// [[Rcpp::export]]
DataFrame calPepAtomCount(StringVector AAstrs)
{
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	averagine mAveragine(ProNovoConfig::getMinPeptideLength(),
						 ProNovoConfig::getMaxPeptideLength());
	vector<int> C(AAstrs.size(), 0);
	vector<int> H, O, N, P, S;
	H = O = N = P = S = C;
	for (int i = 0; i < AAstrs.size(); i++)
	{
		mAveragine.calPepAtomCounts(as<std::string>((AAstrs[i])));
		C[i] = mAveragine.pepAtomCounts[0];
		H[i] = mAveragine.pepAtomCounts[1];
		O[i] = mAveragine.pepAtomCounts[2];
		N[i] = mAveragine.pepAtomCounts[3];
		P[i] = mAveragine.pepAtomCounts[4];
		S[i] = mAveragine.pepAtomCounts[5];
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
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	averagine mAveragine(ProNovoConfig::getMinPeptideLength(),
						 ProNovoConfig::getMaxPeptideLength());
	List pepBYs(AAstrs.size());
	for (int i = 0; i < (int)AAstrs.size(); i++)
	{
		mAveragine.calBYionBaseMasses(as<std::string>((AAstrs[i])));
		int BYionsSize = mAveragine.BionsBaseMasses.size() + mAveragine.YionsBaseMasses.size();
		vector<int> C(BYionsSize, 0);
		vector<int> H, O, N, P, S;
		H = O = N = P = S = C;
		vector<string> BYkinds(BYionsSize);
		vector<double> BYbaseMasses(BYionsSize);
		for (size_t j = 0; j < mAveragine.BionsBaseMasses.size(); j++)
		{
			C[j] = mAveragine.BionsAtomCounts[j][0];
			H[j] = mAveragine.BionsAtomCounts[j][1];
			O[j] = mAveragine.BionsAtomCounts[j][2];
			N[j] = mAveragine.BionsAtomCounts[j][3];
			P[j] = mAveragine.BionsAtomCounts[j][4];
			S[j] = mAveragine.BionsAtomCounts[j][5];
			BYkinds[j] = "B" + to_string(j + 1);
			BYbaseMasses[j] = mAveragine.BionsBaseMasses[j];
		}
		int start = mAveragine.BionsBaseMasses.size();
		for (size_t j = 0; j < mAveragine.YionsBaseMasses.size(); j++)
		{
			C[j + start] = mAveragine.YionsAtomCounts[j][0];
			H[j + start] = mAveragine.YionsAtomCounts[j][1];
			O[j + start] = mAveragine.YionsAtomCounts[j][2];
			N[j + start] = mAveragine.YionsAtomCounts[j][3];
			P[j + start] = mAveragine.YionsAtomCounts[j][4];
			S[j + start] = mAveragine.YionsAtomCounts[j][5];
			BYkinds[j + start] = "Y" + to_string(j + 1);
			BYbaseMasses[j + start] = mAveragine.YionsBaseMasses[j];
		}
		DataFrame df = DataFrame::create(Named("C") = C, _("H") = H,
										 _("O") = O, _("N") = N,
										 _("P") = P, _("S") = S, _("Kind") = BYkinds,
										 _("BaseMass") = BYbaseMasses);
		pepBYs[i] = df;
	}
	pepBYs.names() = AAstrs;
	ProNovoConfig::unSetFilename();
	return pepBYs;
}

//' Simple calculator of peptide precursor mass by binomial NP
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
	bool goodInput = true;
	if (Probs.size() != AAstrs.size())
	{
		Rcout << "lengths of AAstr and Probs are not equal!" << endl;
		goodInput = false;
	}
	for (R_xlen_t i = 0; i < Probs.size(); i++)
	{
		if (Probs[i] < 0 || Probs[i] > 1)
		{
			Rcout << "Wrong isotopic percentage!" << endl;
			goodInput = false;
		}
	}
	char cAtom = 'C';
	if (Atom == "C13")
		cAtom = 'C';
	else if (Atom == "H2")
		cAtom = 'H';
	else if (Atom == "O18")
		cAtom = 'O';
	else if (Atom == "N15")
		cAtom = 'N';
    else if (Atom == "S34")
		cAtom = 'S';
	else
	{
		goodInput = false;
		Rcout << Atom.get_cstring() << " element not supported!" << endl;
	}
	NumericVector v(AAstrs.size(), 0);
	if (goodInput)
	{
		// read default config
		string config = get_extdata();
		ProNovoConfig::setFilename(config);
		averagine mAveragine(ProNovoConfig::getMinPeptideLength(),
							 ProNovoConfig::getMaxPeptideLength());
		for (int i = 0; i < AAstrs.size(); i++)
		{
			mAveragine.changeAtomSIPabundance(cAtom, Probs[i]);
			v[i] = mAveragine.calPrecursorMass(as<std::string>((AAstrs[i])));
		}
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
	bool goodInput = true;
	if (Probs.size() != AAstrs.size())
	{
		Rcout << "lengths of AAstr and Probs are not equal!" << endl;
		goodInput = false;
	}
	for (R_xlen_t i = 0; i < Probs.size(); i++)
	{
		if (Probs[i] < 0 || Probs[i] > 1)
		{
			Rcout << "Wrong isotopic percentage!" << endl;
			goodInput = false;
		}
	}
	char cAtom = 'C';
	if (Atom == "C13")
		cAtom = 'C';
	else if (Atom == "H2")
		cAtom = 'H';
	else if (Atom == "O18")
		cAtom = 'O';
	else if (Atom == "N15")
		cAtom = 'N';
    else if (Atom == "S34")
		cAtom = 'S';
	else
	{
		goodInput = false;
		Rcout << Atom.get_cstring() << " element not supported!" << endl;
	}
	NumericVector v(AAstrs.size(), 0);
	if (goodInput)
	{
		// read default config
		string config = get_extdata();
		ProNovoConfig::setFilename(config);
		averagine mAveragine(ProNovoConfig::getMinPeptideLength(),
							 ProNovoConfig::getMaxPeptideLength());
		for (int i = 0; i < AAstrs.size(); i++)
		{
			mAveragine.changeAtomSIPabundance(cAtom, Probs[i]);
			v[i] = mAveragine.calNetronMass(as<std::string>((AAstrs[i])));
		}
	}
	return v;
}

//' Simple peak calculator of user defined isotopic distribution of one peptide by averagine
//' @param AAstrs a CharacterVector of peptides
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @return a list of DataFrames of spectra
//' @export
// [[Rcpp::export]]
List precursor_peak_calculator_DIY_averagine(StringVector AAstrs, String Atom,
											 double Prob)
{
	bool goodInput = true;
	if (Prob < 0 || Prob > 1)
	{
		Rcout << "Wrong isotopic percentage!" << endl;
		goodInput = false;
	}
	char cAtom = 'C';
	if (Atom == "C13")
		cAtom = 'C';
	else if (Atom == "H2")
		cAtom = 'H';
	else if (Atom == "O18")
		cAtom = 'O';
	else if (Atom == "N15")
		cAtom = 'N';
    else if (Atom == "S34")
		cAtom = 'S';
	else
	{
		goodInput = false;
		Rcout << Atom.get_cstring() << " element not supported!" << endl;
	}
	List spectraList(AAstrs.size());
	if (goodInput)
	{
		// read default config
		string config = get_extdata();
		ProNovoConfig::setFilename(config);
		averagine mAveragine(ProNovoConfig::getMinPeptideLength(),
							 ProNovoConfig::getMaxPeptideLength());
		mAveragine.changeAtomSIPabundance(cAtom, Prob);
		mAveragine.calAveraginePepAtomCounts();
		mAveragine.calAveraginePepSIPdistributions();

		IsotopeDistribution mSIP;
		DataFrame df;
		for (int i = 0; i < AAstrs.size(); i++)
		{
			mAveragine.calPrecursorIsotopeDistribution(as<std::string>(AAstrs(i)), mSIP);
			df =
				DataFrame::create(Named("Mass") = move(mSIP.vMass), _["Prob"] = move(mSIP.vProb));
			spectraList[i] = move(df);
		}
	}
	return spectraList;
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
	if (Prob < 0 || Prob > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom, Prob);
    string AAstr_str = AAstr.get_cstring();
    // AA string format is [AAKRCI] for example
	AAstr_str = "[" + AAstr_str + "]";
	vector<vector<double>> vvdYionMass, vvdYionProb, vvdBionMass, vvdBionProb;
	ProNovoConfig::configIsotopologue.computeProductIon(AAstr_str, vvdYionMass,
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
		DataFrame::create(Named("Mass") = move(masses),
						  _["Prob"] = move(probs), _["Kind"] = move(kinds));
	return df;
}
