#include "PSMsFiltrator.h"

PSMsFiltrator::PSMsFiltrator(const string &msipPath, const string &mftPath)
	: sipPath(msipPath), ftPath(mftPath)
{
}

PSMsFiltrator::PSMsFiltrator(const string &workPath)
	: sipPath(workPath), ftPath(workPath)
{
}

PSMsFiltrator::~PSMsFiltrator()
{
}

void PSMsFiltrator::splitString(const string &mString)
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

tuple<int, bool, string> PSMsFiltrator::detectProDecoy(const string &proteinName)
{
	string targetName, decoyName;
	int targetCount = 0, decoyCount = 0;
	bool isDecoy = true;
	// remove {} out of protein names then split it
	splitString(proteinName.substr(1, proteinName.size() - 2));
	// if one protein is not decoy the peptide is not decoy
	for (string token : tokens)
	{
		if (token.substr(0, 4) != "Rev_")
		{
			isDecoy = false;
			targetName += token + "\t";
			targetCount++;
		}
		else
		{
			decoyName += token + "\t";
			decoyCount++;
		}
	}
	decoyName.pop_back();
	targetName.pop_back();
	if (isDecoy)
		return {decoyCount, true, decoyName};
	else
		return {targetCount, false, targetName};
}

float PSMsFiltrator::getRentionTime()
{
	return 0;
}

int PSMsFiltrator::getPct(const string &searchName)
{
	string pct;
	for (size_t i = 0; i < searchName.length(); i++)
	{
		if (searchName[i] == '_')
		{
			pct = searchName.substr(i + 1, 2);
			if (pct[1] < '0' || pct[1] > '9')
				pct.pop_back();
			return (stoi(pct));
		}
	}
	return (0);
}

tuple<int, string, string> PSMsFiltrator::getRealPep(const string &pepSeq)
{
	int length = pepSeq.length() - 2;
	int start = 1;
	string flankN = "n", flankC = "c";
	if (pepSeq[0] != '[')
	{
		length--;
		start++;
		flankN = pepSeq[0];
	}
	if (pepSeq.back() != ']')
	{
		length--;
		flankC = pepSeq.back();
	}
	string realPep = pepSeq.substr(start, length);
	return {length, realPep, flankN + "." + realPep + "." + flankC};
}

pair<vector<double>, vector<double>> PSMsFiltrator::getPrecursorIsotopicPeak()
{
	return pair<vector<double>, vector<double>>{};
}

sipPSMinfo PSMsFiltrator::convertFilesScansTopPSMs(
	const unordered_map<string, unordered_map<int, vector<topPSM>>> &filesScansTopPSMs)
{
	// converted from filesScansTopPSMs
	sipPSMinfo topPSMs;
	for (auto &fileIX : filesScansTopPSMs)
	{
		for (auto &scanIX : fileIX.second)
		{
			for (size_t i = 0; i < scanIX.second.size(); i++)
			{
				topPSMs.fileNames.push_back(fileIX.first);
				topPSMs.scanNumbers.push_back(scanIX.first);
				topPSMs.parentCharges.push_back(scanIX.second[i].parentCharge);
				topPSMs.measuredParentMasses.push_back(scanIX.second[i].measuredParentMass);
				topPSMs.calculatedParentMasses.push_back(scanIX.second[i].calculatedParentMass);
				topPSMs.searchNames.push_back(scanIX.second[i].searchName);
				topPSMs.scores.push_back(scanIX.second[i].score);
				topPSMs.identifiedPeptides.push_back(scanIX.second[i].identifiedPeptide);
				topPSMs.originalPeptides.push_back(scanIX.second[i].originalPeptide);
				topPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
				topPSMs.retentionTimes.push_back(getRentionTime());
				int proteinCount;
				bool isDecoy;
				string trimedProName;
				tie(proteinCount, isDecoy, trimedProName) = detectProDecoy(scanIX.second[i].proteinName);
				topPSMs.isDecoys.push_back(isDecoy);
				topPSMs.proCounts.push_back(proteinCount);
				topPSMs.trimedProteinNames.push_back(trimedProName);
				topPSMs.ranks.push_back(i + 1);
				int pct = getPct(scanIX.second[i].searchName);
				topPSMs.pcts.push_back(pct);
				int pepLen;
				string realPep, formatedPep;
				tie(pepLen, realPep, formatedPep) = getRealPep(scanIX.second[i].originalPeptide);
				topPSMs.pepLengths.push_back(pepLen);
				topPSMs.realPepSeqs.push_back(realPep);
				topPSMs.formatedPepSeqs.push_back(formatedPep);
				topPSMs.psmIDs.push_back(fileIX.first + "_" +
										 to_string(scanIX.first) + "_" +
										 scanIX.second[i].identifiedPeptide + "_" + to_string(pct));
				// cout << proteinCount << endl;
			}
		}
	}
	return topPSMs;
}

void PSMsFiltrator::writePercolatorTSV()
{
}

void PSMsFiltrator::readPercolatorTSV()
{
}
