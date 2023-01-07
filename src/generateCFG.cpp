#include "lib/cfgParser.h"
#include <Rcpp.h>
using namespace Rcpp;

//' generateOneCFG
//' @param cfgPath a full path of .cfg file
//' @param outPath a full path for .cfg file output
//' @param element a string of element name, "N" for example
//' @param pct a integer of element SIP abundance
//' @param center a integer of mass window center
//' @param width a integer of mass half window width
//' @return a bool value if generate succeed or not
//' @export
// [[Rcpp::export]]
bool generateOneCFG(String cfgPath, String outPath, String element,
                    int pct, int center, int width)
{
    cfgParser parser(cfgPath);
    // make sure read cfg file succeed
    if (parser.lines.size() == 0)
        return false;
    if (element != "C" && element != "N")
        return false;
    parser.setSearch_NameIX();
    parser.changeSearchName(to_string(pct));
    parser.setParent_Mass_WindowsIX();
    parser.changeMassWindowsCenter(center, width);
    parser.setElement_PercentIX(element);
    parser.changeSIPabundance(pct);
    parser.writeFile(outPath);
    return true;
}

//' generateCFGs
//' @param cfgPath a full path of .cfg file
//' @param outPath a full path for .cfg file output
//' @param element a string of element name, "N" for example
//' @param pct a integer of element SIP abundance
//' @param center a integer of mass window center
//' @param width a integer of mass half window width
//' @return a bool value if generate succeed or not
//' @export
// [[Rcpp::export]]
bool generateCFGs(String cfgPath, String outPath, String element)
{
    cfgParser parser(cfgPath);
    // make sure read cfg file succeed
    if (parser.lines.size() == 0)
        return false;
    parser.setSearch_NameIX();
    parser.setParent_Mass_WindowsIX();
    vector<int> centers, widths;
    if (element == "C")
    {
        centers = {0, 0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 6, 6,
                   7, 7, 8, 9, -2, -1, -1, 0, 0, 1, 1, 2, 3, 3, 4,
                   4, 5, 5, 6, 7, 7, -3, -3, -2, -2, -1, -1, 0, 1,
                   1, 2, 2, 3, 3, 4, 4, 5, -5, -5, -4, -4, -3, -3,
                   -2, -1, -1, 0, 0, 1, 1, 2, 2, 3, 4, -7, -6, -6,
                   -5, -5, -4, -4, -3, -2, -2, -1, -1, 0, 0, 1, 2,
                   2, -8, -8, -7, -7, -6, -6, -5, -4, -4, -3, -3,
                   -2, -2, -1, -1, 0, 0};
        widths = {2, 2, 2, 3, 4, 4, 4, 5, 4, 5, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 6, 7, 7, 7, 8, 7, 8, 7, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
                  8, 8, 8, 8, 8, 7, 8, 7, 7, 7, 7, 6, 7, 6, 6, 6, 6, 6,
                  6, 6, 6, 6, 5, 5, 4, 4, 4, 4, 4, 3, 3, 2};
    }
    else if (element == "N")
    {
        centers = {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2,
                   2, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
                   5, 5, 6, 6, 6, -5, -5, -5, -5, -5, -4, -4, -4, -4, -4, -4, -4, -4,
                   -4, -3, -3, -3, -3, -3, -3, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2,
                   -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        widths = {2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                  4, 4, 4, 4, 4, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                  4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                  4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 2};
    }
    else
    {
        Rcout << element.get_cstring() << " element is not supported!" << endl;
        return false;
    }
    parser.setElement_PercentIX(element);
    for (size_t pct = 0; pct < 101; pct++)
    {
        parser.changeSearchName(to_string(pct) + "Pct");
        parser.changeMassWindowsCenter(centers[pct], widths[pct]);
        // set center 0
        // parser.changeMassWindowsCenter(0, widths[pct]);
        parser.changeSIPabundance(pct);
        parser.writeFile(outPath);
    }
    return true;
}