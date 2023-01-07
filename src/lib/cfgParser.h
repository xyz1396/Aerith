#pragma once
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <filesystem>
#include <unordered_map>
#include <algorithm>
namespace fs = std::filesystem;
using namespace std;

class cfgParser
{
private:
public:
    ifstream cfgFileStream;
    vector<string> lines;
    vector<string> tokens;
    // unordered_map<parameter name, splited parameter content>
    unordered_map<string, vector<string>> parametersMap;
    // unordered_map<parameter name, parameter line IX>
    unordered_map<string, size_t> parametersIXMap;
    size_t Search_NameIX, Parent_Mass_WindowsIX, Element_PercentIX;
    string Element_PercentName, newFileName;
    cfgParser(const string &cfgFileName);
    ~cfgParser();
    void splitString(const string &mString);
    void parseParameters();
    size_t findParameter(const string &parameter);
    void setSearch_NameIX();
    void setParent_Mass_WindowsIX();
    void setElement_PercentIX(const string &element);
    void changeSearchName(const string &nameSuffix);
    void changeMassWindowsCenter(const int center, const int windowsSize);
    void changeSIPabundance(const int sipAbundance);
    void writeFile(const string &folderPath);
};
