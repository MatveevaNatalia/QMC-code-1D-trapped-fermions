#ifndef STATISTICS_H
#define STATISTICS_H

#include <string>
using namespace std;
void Normalization(string name_file, string name_file_1, int num, int np, int nc);
void Stat_Energy(string name_file, string name_file_1, string name_file_2, int ncol, int nel, int n_notused,  int nfull);
void statistics_cor(string name_file, string name_file_1, int npt, int ncomp, int nfora, int nel);

#endif
