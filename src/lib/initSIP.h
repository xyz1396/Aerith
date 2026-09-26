#pragma once
#include "parameters.h"
#include "PeptideIsotopeCalculator.h"
#include <Rcpp.h>

using namespace Rcpp;

inline void computeResidueMassIntensityAgain(const string isotope, double abundance)
{
    const char atom = AerithParameters::isotopeElement(isotope);
    AerithParameters::validateAbundance(abundance);
    PeptideIsotopeCalculator calculator;
    calculator.changeAtomSIPabundance(atom, abundance);
    AerithParameters::current().isotopologue.refreshResidueDistributions();
}
