#include "peptidesFiltrator.h"

peptideInfo::peptideInfo(bool mIsDecoy, float mBestScore)
	: isDecoy(mIsDecoy), bestScore(mBestScore)
{
}

peptidesFiltrator::peptidesFiltrator(const vector<sipPSM> &sipPSMs, float mFDRthreshold)
	: FDRthreshold(mFDRthreshold)
{
	for (sipPSM psm : sipPSMs)
	{
		for (size_t i = 0; i < psm.identifiedPeptides.size(); i++)
		{
			if (psm.parentCharges[i] <= 2)
			{
				initFillPeptideMap(psm.identifiedPeptides[i], psm.scores[i], psm.proteinNames[i],
								   peptideMapCharge2);
			}
			else if (psm.parentCharges[i] == 3)
			{
				initFillPeptideMap(psm.identifiedPeptides[i], psm.scores[i], psm.proteinNames[i],
								   peptideMapCharge3);
			}
			else
			{
				initFillPeptideMap(psm.identifiedPeptides[i], psm.scores[i], psm.proteinNames[i],
								   peptideMapChargeLargerThan3);
			}
		}
	}
}

peptidesFiltrator::~peptidesFiltrator()
{
}

void peptidesFiltrator::splitString(const string mString)
{
	string sep = ",";
	size_t start = 0;
	size_t end = mString.find(sep);
	tokens.clear();
	while (end != std::string::npos)
	{
		tokens.push_back(mString.substr(start, end - start));
		start = end + sep.length();
		end = mString.find(sep, start);
	}
	tokens.push_back(mString.substr(start));
}

void peptidesFiltrator::initFillPeptideMap(const string &pepSeq,
									   const float score,
									   const string &proName,
									   unordered_map<string, peptideInfo> &mMap)
{
	auto pepIX = mMap.find(pepSeq);
	if (pepIX != mMap.end())
	{
		pepIX->second.scores.push_back(score);
		if (score > pepIX->second.bestScore)
			pepIX->second.bestScore = score;
		// if one protein is not decoy the peptide is not decoy
		// check it is decoy or not only when it is decoy before
		if (pepIX->second.isDecoy)
		{
			if (!detectDecoy(proName))
				pepIX->second.isDecoy = false;
		}
	}
	else
	{
		peptideInfo temp = peptideInfo(detectDecoy(proName), score);
		temp.scores.push_back(score);
		mMap.insert({pepSeq, temp});
	}
}

bool peptidesFiltrator::detectDecoy(const string &proteinName)
{
	// remove {} out of protein names then split it
	splitString(proteinName.substr(1, proteinName.size() - 2));
	// if one protein is not decoy the peptide is not decoy
	for (string token : tokens)
	{
		if (token.substr(0, 4) != "Rev_")
			return false;
	}
	return true;
}

void peptidesFiltrator::
	sortPeptideBestScore(const unordered_map<string, peptideInfo> &mPepMap, vector<pair<float, bool>> &bestScoreDecoyPairs)
{
	// convert unordered_map to vector of pairs
	bestScoreDecoyPairs.clear();
	bestScoreDecoyPairs.resize(mPepMap.size());
	size_t i = 0;
	for (auto &pepIX : mPepMap)
	{
		bestScoreDecoyPairs[i] = {pepIX.second.bestScore, pepIX.second.isDecoy};
		i++;
	}
	// descending sort
	sort(bestScoreDecoyPairs.begin(), bestScoreDecoyPairs.end(),
		 [](const pair<float, bool> &a, const pair<float, bool> &b) -> bool
		 {
			 return a.first > b.first;
		 });
}

void peptidesFiltrator::sortPeptideAllScore(const unordered_map<string, peptideInfo> &mPepMap, vector<pair<float, pair<string, bool>>> &scoreDecoyPairs)
{
	// convert unordered_map to vector of pairs
	scoreDecoyPairs.clear();
	for (auto &pepIX : mPepMap)
	{
		for (float score : pepIX.second.scores)
		{
			scoreDecoyPairs.push_back({score, {pepIX.first, pepIX.second.isDecoy}});
		}
	}
	// descending sort
	sort(scoreDecoyPairs.begin(), scoreDecoyPairs.end(),
		 [](const pair<float, pair<string, bool>> &a, const pair<float, pair<string, bool>> &b) -> bool
		 {
			 return a.first > b.first;
		 });
}

tuple<size_t, size_t, float> peptidesFiltrator::getDecoyCountScoreThreshold(vector<pair<float, bool>> &bestScoreDecoyPairs)
{
	size_t decoyCount = 0;
	vector<float> FDRs(bestScoreDecoyPairs.size());
	vector<int> decoyCounts(bestScoreDecoyPairs.size());
	for (size_t i = 0; i < bestScoreDecoyPairs.size(); i++)
	{
		if (bestScoreDecoyPairs[i].second)
			decoyCount++;
		// must convert to float or double
		FDRs[i] = decoyCount / (float)i;
		decoyCounts[i] = decoyCount;
	}
	// find ScoreThreshold from end of FDRs
	size_t i = FDRs.size();
	while (i > 0)
	{
		i--;
		if (FDRs[i] <= FDRthreshold)
			return {decoyCounts[i], i, bestScoreDecoyPairs[i].first};
	}
	// if Cannot find FDRthreshold
	cout << "Cannot find FDRthreshold" << endl;
	return {0, 0, 0};
}

tuple<size_t, size_t, float> peptidesFiltrator::getDecoyCountScoreThreshold(vector<pair<float, pair<string, bool>>> &scoreDecoyPairs)
{
	vector<int> pepCounts(scoreDecoyPairs.size());
	vector<int> decoyCounts(scoreDecoyPairs.size());
	unordered_set<string> pepSeqs;
	unordered_set<string> decoySeqs;
	vector<float> FDRs(scoreDecoyPairs.size());
	size_t i = 0;
	for (auto ix : scoreDecoyPairs)
	{
		pepSeqs.insert(ix.second.first);
		if (ix.second.second)
			decoySeqs.insert(ix.second.first);
		pepCounts[i] = pepSeqs.size();
		decoyCounts[i] = decoySeqs.size();
		// must convert to float or double
		FDRs[i] = decoySeqs.size() / (float)pepSeqs.size();
		i++;
	}
	// find ScoreThreshold from end of FDRs
	while (i > 0)
	{
		i--;
		if (FDRs[i] <= FDRthreshold)
			return {decoyCounts[i], pepCounts[i], scoreDecoyPairs[i].first};
	}
	// if Cannot find FDRthreshold
	cout << "Cannot find FDRthreshold" << endl;
	return {0, 0, 0};
}

void peptidesFiltrator::filterPeptideMap()
{
	// temp bestScoreDecoyPairs for sort
	// vector<pair<float, bool>> bestScoreDecoyPairs;
	// sortPeptideBestScore(peptideMapCharge2, bestScoreDecoyPairs);
	// tie(decoyCountCharge2, pepCountCharge2, scoreThresholdCharge2) = getDecoyCountScoreThreshold(bestScoreDecoyPairs);
	// sortPeptideBestScore(peptideMapCharge3, bestScoreDecoyPairs);
	// tie(decoyCountCharge3, pepCountCharge3, scoreThresholdCharge3) = getDecoyCountScoreThreshold(bestScoreDecoyPairs);
	// sortPeptideBestScore(peptideMapChargeLargerThan3, bestScoreDecoyPairs);
	// tie(decoyCountChargeLargerThan3, pepCountChargeLargerThan3, scoreThresholdChargeLargerThan3) = getDecoyCountScoreThreshold(bestScoreDecoyPairs);

	vector<pair<float, pair<string, bool>>> scoreDecoyPairs;
	sortPeptideAllScore(peptideMapCharge2, scoreDecoyPairs);
	tie(decoyCountCharge2, pepCountCharge2, scoreThresholdCharge2) = getDecoyCountScoreThreshold(scoreDecoyPairs);
	sortPeptideAllScore(peptideMapCharge3, scoreDecoyPairs);
	tie(decoyCountCharge3, pepCountCharge3, scoreThresholdCharge3) = getDecoyCountScoreThreshold(scoreDecoyPairs);
	sortPeptideAllScore(peptideMapChargeLargerThan3, scoreDecoyPairs);
	tie(decoyCountChargeLargerThan3, pepCountChargeLargerThan3, scoreThresholdChargeLargerThan3) = getDecoyCountScoreThreshold(scoreDecoyPairs);
}

void peptidesFiltrator::fillPeptideMapExtraInfo()
{
}