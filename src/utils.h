#ifndef UTILS_H
#define UTILS_H
#include <map>
#include <string>
using namespace std;

map<string, double> FillMap(const string & inFile, const string & outDir);
void Normalization(string name_file, string name_file_1, int num, int np, int nc);
void PrintAcceptance(double accrate, int iblck, const string & outDir );
#endif // UTILS_H
