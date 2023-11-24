#include "PSMsFiltrator.h"

PSMsFiltrator::PSMsFiltrator(const std::string &msipPath, const std::string &mftPath)
	: sipPath(msipPath), ftPath(mftPath)
{
}

PSMsFiltrator::PSMsFiltrator(const std::string &workPath)
	: sipPath(workPath), ftPath(workPath)
{
}

PSMsFiltrator::~PSMsFiltrator()
{
}

void PSMsFiltrator::splitString(const std::string &mString)
{
	std::string sep = ",";
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

std::tuple<int, bool, std::string> PSMsFiltrator::detectProDecoy(const std::string &proteinName)
{
	std::string targetName, decoyName;
	int targetCount = 0, decoyCount = 0;
	bool isDecoy = true;
	// remove {} out of protein names then split it
	splitString(proteinName.substr(1, proteinName.size() - 2));
	// if one protein is not decoy the peptide is not decoy
	for (std::string token : tokens)
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

int PSMsFiltrator::getPct(const std::string &searchName)
{
	std::string pct;
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

std::tuple<int, std::string, std::string> PSMsFiltrator::getRealPep(const std::string &pepSeq)
{
	int length = pepSeq.length() - 2;
	int start = 1;
	std::string flankN = "n", flankC = "c";
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
	std::string realPep = pepSeq.substr(start, length);
	return {length, realPep, flankN + "." + realPep + "." + flankC};
}

std::pair<std::vector<double>, std::vector<double>> PSMsFiltrator::getPrecursorIsotopicPeak()
{
	return std::pair<std::vector<double>, std::vector<double>>{};
}

sipPSMinfo PSMsFiltrator::convertFilesScansTopPSMs(
	const std::unordered_map<std::string, std::unordered_map<int, std::vector<scanTopPSM>>> &filesScansTopPSMs)
{
	// converted from filesScansTopPSMs
	sipPSMinfo scanTopPSMs;
	for (auto &fileIX : filesScansTopPSMs)
	{
		for (auto &scanIX : fileIX.second)
		{
			for (size_t i = 0; i < scanIX.second.size(); i++)
			{
				scanTopPSMs.fileNames.push_back(fileIX.first);
				scanTopPSMs.scanNumbers.push_back(scanIX.first);
				scanTopPSMs.parentCharges.push_back(scanIX.second[i].parentCharge);
				scanTopPSMs.measuredParentMasses.push_back(scanIX.second[i].measuredParentMass);
				scanTopPSMs.calculatedParentMasses.push_back(scanIX.second[i].calculatedParentMass);
				scanTopPSMs.searchNames.push_back(scanIX.second[i].searchName);
				scanTopPSMs.scores.push_back(scanIX.second[i].score);
				scanTopPSMs.identifiedPeptides.push_back(scanIX.second[i].identifiedPeptide);
				scanTopPSMs.originalPeptides.push_back(scanIX.second[i].originalPeptide);
				scanTopPSMs.proteinNames.push_back(scanIX.second[i].proteinName);
				scanTopPSMs.retentionTimes.push_back(getRentionTime());
				int proteinCount;
				bool isDecoy;
				std::string trimedProName;
				tie(proteinCount, isDecoy, trimedProName) = detectProDecoy(scanIX.second[i].proteinName);
				scanTopPSMs.isDecoys.push_back(isDecoy);
				scanTopPSMs.proCounts.push_back(proteinCount);
				scanTopPSMs.trimedProteinNames.push_back(trimedProName);
				scanTopPSMs.ranks.push_back(i + 1);
				int pct = getPct(scanIX.second[i].searchName);
				scanTopPSMs.pcts.push_back(pct);
				int pepLen;
				std::string realPep, formatedPep;
				tie(pepLen, realPep, formatedPep) = getRealPep(scanIX.second[i].originalPeptide);
				scanTopPSMs.pepLengths.push_back(pepLen);
				scanTopPSMs.realPepSeqs.push_back(realPep);
				scanTopPSMs.formatedPepSeqs.push_back(formatedPep);
				scanTopPSMs.psmIDs.push_back(fileIX.first + "_" +
											 std::to_string(scanIX.first) + "_" +
											 scanIX.second[i].identifiedPeptide + "_" + std::to_string(pct));
				// cout << proteinCount << endl;
			}
		}
	}
	return scanTopPSMs;
}

void PSMsFiltrator::writePercolatorTSV()
{
}

void PSMsFiltrator::readPercolatorTSV()
{
}
