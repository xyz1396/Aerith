#include "averagine.h"

averagine::averagine(const int minPepLen, const int maxPepLen)
    : minPepLen(minPepLen), maxPepLen(maxPepLen), pepLenRange(maxPepLen - minPepLen + 1)
{
    // get averagine's Atom content and mass distribution
    averagineAtomCounts = ProNovoConfig::configIsotopologue.mResidueAtomicComposition[averagineResidue];
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(averagineAtomCounts,
                                                                  averagineSIPdistribution);
    diffAtomCounts.resize(averagineAtomCounts.size());
}

averagine::~averagine()
{
}

void averagine::changeAtomSIPabundance(const char SIPatom, const double pct)
{
    bool SIPatomFound = false;
    for (size_t i = 0; i < SIPatoms.size(); i++)
    {
        if (SIPatoms[i] == SIPatom)
        {
            SIPatomIX = i;
            SIPatomFound = true;
            break;
        }
    }
    if (SIPatomFound)
    {
        ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[SIPatomIX].vProb[0] =
            1.0 - pct;
        ProNovoConfig::configIsotopologue.vAtomIsotopicDistribution[SIPatomIX].vProb[1] =
            pct;
        ProNovoConfig::configIsotopologue.computeIsotopicDistribution(
            ProNovoConfig::configIsotopologue.mResidueAtomicComposition[averagineResidue],
            averagineSIPdistribution);
    }
    else
    {
        cout << SIPatom << "element is not support!" << endl;
    }
}

void averagine::calAveraginePepAtomCounts()
{
    vector<int> averaginePepAtomCounts;
    averaginePepAtomCountss.reserve(pepLenRange);
    averaginePepAtomCounts.reserve(averagineAtomCounts.size());
    for (int pepLen = minPepLen; pepLen <= maxPepLen; pepLen++)
    {
        averaginePepAtomCounts.clear();
        for (size_t j = 0; j < averagineAtomCounts.size(); j++)
            averaginePepAtomCounts.push_back(averagineAtomCounts[j] * pepLen);
        averaginePepAtomCountss.push_back(averaginePepAtomCounts);
    }
}

vector<int> *averagine::getAveraginePepAtomCounts(const int pepLen)
{
    return &averaginePepAtomCountss[pepLen - minPepLen];
}

void averagine::calAveraginePepSIPdistributions()
{
    IsotopeDistribution NtermSIPdistribution =
        ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution["Nterm"];
    IsotopeDistribution CtermSIPdistribution =
        ProNovoConfig::configIsotopologue.vResidueIsotopicDistribution["Cterm"];
    averaginePepSIPdistributions.reserve(pepLenRange);
    IsotopeDistribution tempSIPdistribution = NtermSIPdistribution;
    ProNovoConfig::configIsotopologue.sum(CtermSIPdistribution, tempSIPdistribution);
    for (int i = 0; i < minPepLen; i++)
    {
        ProNovoConfig::configIsotopologue.sum(averagineSIPdistribution, tempSIPdistribution);
    }
    for (int i = 0; i <= pepLenRange; i++)
    {
        ProNovoConfig::configIsotopologue.sum(averagineSIPdistribution, tempSIPdistribution);
        averaginePepSIPdistributions.push_back(tempSIPdistribution);
    }
}

IsotopeDistribution *averagine::getAveraginePepSIPdistribution(const int pepLen)
{
    return &averaginePepSIPdistributions[pepLen - minPepLen];
}

void averagine::calDiffAtomCounts(const string &pepSeq)
{
    vector<int> residueAtomCounts;
    for (size_t i = 0; i < pepSeq.size(); i++)
    {
        residueAtomCounts = ProNovoConfig::configIsotopologue
                                .mResidueAtomicComposition[pepSeq.substr(i, 1)];
        for (size_t j = 0; j < residueAtomCounts.size(); j++)
            diffAtomCounts[j] += residueAtomCounts[j];
    }
    vector<int> *averaginePepAtomCountPtr = getAveraginePepAtomCounts(pepSeq.length());
    for (size_t k = 0; k < diffAtomCounts.size(); k++)
        diffAtomCounts[k] -= (*averaginePepAtomCountPtr)[k];
}

void averagine::calPrecusorIsotopeDistribution(const string &pepSeq, IsotopeDistribution &tempSIPdistribution)
{
    calDiffAtomCounts(pepSeq);
    ProNovoConfig::configIsotopologue.computeIsotopicDistribution(diffAtomCounts,
                                                                  tempSIPdistribution);
    ProNovoConfig::configIsotopologue.sum(*(getAveraginePepSIPdistribution(pepSeq.length())),
                                          tempSIPdistribution);
}