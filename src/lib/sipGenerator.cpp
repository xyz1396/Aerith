#include "sipGenerator.h"

sipGenerator::sipGenerator(/* args */)
{
    gen = mt19937(rd());
}

sipGenerator::~sipGenerator()
{
}

void sipGenerator::initElementDistributions()
{
    for (auto elementMassAbundance : elementMassAbundances)
    {
        elementDistributions.push_back(discrete_distribution<>(elementMassAbundance.second.begin(),
                                                               elementMassAbundance.second.end()));
    }
}

double sipGenerator::generateMass(const vector<int> &elementIXs, const vector<int> &elementCounts)
{
    double mass = 0;
    for (size_t i = 0; i < elementIXs.size(); i++)
    {
        for (int j = 0; j < elementCounts[i]; j++)
        {
            mass += elementMassAbundances[elementIXs[i]].first[elementDistributions[elementIXs[i]](gen)];
        }
    }
    return mass;
}