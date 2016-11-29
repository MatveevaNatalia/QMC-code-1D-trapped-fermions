#ifndef OBDM_H
#define OBDM_H

#include "qmc.h"
#include "Locals.h"
#include "MomDistr.h"
#include "Wave_fun.h"
#include "ConfigDistr.h"
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

class OBDM{
private:
    int num_points;
    double step;
    ConfigDistr auxil, average;
    ConfigDistr auxil_norm,  average_norm;
    vector<ConfigDistr> oldPage, newPage;
    vector<ConfigDistr> oldPage_norm, newPage_norm;
public:
    OBDM(const CorFunParam& );
    void SetZero();
    void SetZeroAx();
    void NotAccept(int ipop);
    void WalkerMatch();
    void WalkerCollect(int nsons);
    void Normalization();
    void PrintDistr(const string& name_file);
    void OBDM_Calc( ParamModel& param_model, const Configuration& xaux, WaveFunction & wave_func, MomentDistr& moment_distr, const CorFunParam&  mom_distr_param);
    void PageSwap();

};

#endif
