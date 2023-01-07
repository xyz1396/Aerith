#pragma once
#include "proNovoConfig.h"
#include <random>

using namespace std;

class sipGenerator
{
private:
    random_device rd;
    mt19937 gen;
    // <elements<elemen isotopic masses, element isotopic abundances converted to integer>>
    vector<pair<vector<double>, vector<int>>> elementMassAbundances;
    vector<discrete_distribution<>> elementDistributions;

public:
    sipGenerator(/* args */);
    ~sipGenerator(); 
    // append one element's isotopic mass and abundance to elementMassAbundanceTable
    void appendElementMassAbundanceTable(const vector<double> &masses, const vector<double> &probs);
    void initElementDistributions(); 
    double generateMass(const vector<int> &elementIXs, const vector<int> &elementCounts);
};
