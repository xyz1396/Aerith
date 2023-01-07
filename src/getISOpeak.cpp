#include "lib/initSIP.h"
#include "lib/averagine.h"
#include <Rcpp.h>

using namespace Rcpp;

//' Simple peak calculator of natural isotopic distribution
//' @param AAstr a CharacterVector
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator(CharacterVector AAstr)
{
	if (AAstr.length() > 1)
		Rcerr << "only one string one time" << endl;
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	IsotopeDistribution myIso;
	string AAstr_str = as<std::string>(AAstr);
	ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr_str,
																  myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' Simple residue peak calculator of user defined isotopic distribution of one residue
//' @param residue residue name
//' @param Atom "C13" or "N15"
//' @param Prob its SIP abundance (0.0~1.0)
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

//' Simple peak calculator of user defined isotopic distribution of one peptide
//' @param AAstr a CharacterVector
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @export
// [[Rcpp::export]]
DataFrame precursor_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom,
										NumericVector Prob)
{
	// check and convert input
	if (AAstr.length() > 1 || Atom.length() > 1 || Prob.length() > 1)
		Rcerr << "only one string one time" << endl;
	string Atom_str = as<string>(Atom);
	double Prob_d = as<double>(Prob);
	if (Prob_d < 0 || Prob_d > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom_str, Prob_d);
	// test if prob change take effct
	// ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].print();
	// ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].print();
	IsotopeDistribution myIso;
	string AAstr_str = as<std::string>(AAstr);
	ProNovoConfig::configIsotopologue.computeIsotopicDistribution(AAstr_str,
																  myIso);
	DataFrame df =
		DataFrame::create(Named("Mass") = myIso.vMass, _["Prob"] = myIso.vProb);
	return df;
}

//' Simple peak calculator of user defined isotopic distribution of one peptide by averagine
//' @param AAstr a CharacterVector of peptides
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
	char cAtom;
	if (Atom == "C13")
		cAtom = 'C';
	else if (Atom == "N15")
		cAtom = 'N';
	else
	{
		goodInput = false;
		cout << Atom.get_cstring() << " element not support!" << endl;
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
			mAveragine.calPrecusorIsotopeDistribution(as<std::string>(AAstrs(i)), mSIP);
			df =
				DataFrame::create(Named("Mass") = move(mSIP.vMass), _["Prob"] = move(mSIP.vProb));
			spectraList[i] = move(df);
		}
	}
	return spectraList;
}

//' peak calculator of B Y ione from of one peptide using user defined isotopic distribution
//' @param AAstr a CharacterVector
//' @param Atom a CharacterVector C13 or N15
//' @param Prob a NumericVector for its abundance
//' @export
// [[Rcpp::export]]
DataFrame BYion_peak_calculator_DIY(CharacterVector AAstr, CharacterVector Atom,
									NumericVector Prob)
{
	// check and convert input
	if (AAstr.length() > 1 || Atom.length() > 1 || Prob.length() > 1)
		Rcerr << "only one string one time" << endl;
	string AAstr_str = as<std::string>(AAstr);
	string Atom_str = as<string>(Atom);
	double Prob_d = as<double>(Prob);
	if (Prob_d < 0 || Prob_d > 1)
		Rcout << "Wrong isotopic percentage" << endl;
	// read default config
	string config = get_extdata();
	ProNovoConfig::setFilename(config);
	// compute residue mass and prob again
	computeResidueMassIntensityAgain(Atom_str, Prob_d);
	// AA string format is [AAKRCI]
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
