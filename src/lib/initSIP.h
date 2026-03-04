#pragma once
#include "proNovoConfig.h"
#include "averagine.h"
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
	char sipAtom = '\0';
	if (Atom_str == "C13")
		sipAtom = 'C';
	else if (Atom_str == "H2")
		sipAtom = 'H';
	else if (Atom_str == "O18")
		sipAtom = 'O';
	else if (Atom_str == "N15")
		sipAtom = 'N';
	else if (Atom_str == "S34")
		sipAtom = 'S';
	else
	{
		Rcerr << "this element is not support" << endl;
		return;
	}

	averagine mAveragine;
	mAveragine.changeAtomSIPabundance(sipAtom, Prob_d);

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
