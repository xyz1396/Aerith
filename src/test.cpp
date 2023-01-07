#include "lib/ftFileReader.h"
#include <Rcpp.h>

using namespace Rcpp;

//' Simple sum
//' @param x a numeric vector
//' @export
// [[Rcpp::export]]
double calc_sum(NumericVector x)
{
	double sum = 0;
	for (int i = 0; i < x.size(); ++i)
	{
		std::cout << i << std::endl;
		sum += x[i];
	}
	return sum;
}

//' test ftFileReader
//' @param ftFile a ft file's full path
//' @export
// [[Rcpp::export]]
void test_ftFileReader(CharacterVector ftFile)
{
	ftFileReader reader(as<std::string>(ftFile));
	if (reader.detectPrecursorAndCharge())
		reader.printFileInfo();
	else
		cout << "test failed" << endl;
}
