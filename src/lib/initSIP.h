#pragma once
#include "proNovoConfig.h"
#include <Rcpp.h>

using namespace Rcpp;

// inline make it can be include by multiple .cpp files
inline string get_extdata()
{
	Environment base("package:base");
	Function sys_file = base["system.file"];
	StringVector resVector =
		sys_file("extdata", "SiprosConfig.cfg", _["package"] = "Aerith");
	string res = as<std::string>(resVector);
	return res;
}

inline void computeResidueMassIntensityAgain(const string Atom_str, double Prob_d)
{
	// change Prob
	if (Atom_str == "C13")
	{
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[0] =
			1.0 - Prob_d;
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[0].vProb[1] =
			Prob_d;
		ProNovoConfig::getSetSIPelement() = "C";
	}
	else if (Atom_str == "N15")
	{
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[0] =
			1.0 - Prob_d;
		ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[3].vProb[1] =
			Prob_d;
		ProNovoConfig::getSetSIPelement() = "N";
	}
	else
		Rcerr << "this element is not support" << endl;
	// compute residue mass and prob again
	map<string, vector<int>>::iterator ResidueIter;
	IsotopeDistribution tempIsotopeDistribution;
	for (ResidueIter =
			 ProNovoConfig::configIsotopologue.mResidueAtomicComposition.begin();
		 ResidueIter !=
		 ProNovoConfig::configIsotopologue.mResidueAtomicComposition.end();
		 ResidueIter++)
	{
		if (!ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
				ResidueIter->second, tempIsotopeDistribution))
		{
			Rcerr << "ERROR: cannot calculate the isotopic distribution for residue "
				  << ResidueIter->first << endl;
		}
		ProNovoConfig::configIsotopologue
			.vResidueIsotopicDistribution[ResidueIter->first] =
			tempIsotopeDistribution;
	}
}