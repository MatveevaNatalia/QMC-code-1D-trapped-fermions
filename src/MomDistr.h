#ifndef MOMDISTR_H
#define MOMDISTR_H

#include "qmc.h"
#include <string>
#include "ConfigDistr.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

class MomentDistr{
private:
    ConfigDistr average;
    vector<ConfigDistr> oldPage, newPage;
    int num_points;
    double step;
public:
    ConfigDistr auxil;
    MomentDistr(const CorFunParam& param);
    void SetZero();
    void SetZeroAx();
    void NotAccept(int ipop);
    void WalkerMatch();
    void WalkerCollect(int nsons);
    void Normalization(const ParamModel& param, int norm);
    void PrintDistr(const string& name_file);
    void PageSwap();
};

#endif // MOMDISTR_H
