#ifndef OBDM_H
#define OBDM_H

#include "qmc.h"
#include "Locals.h"
#include "MomDistr.h"
#include "Wave_fun.h"
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
    void setzero(double * x){
        for(int i = 0; i < num_points+1; i++) x[i] = 0;
    }
public:
    double *fr, *fra, *nfr, *nfra, **frlocal, **nfrlocal;
    OBDM(const CorFunParam& obdm);

    void SetZero();
    void SetZeroAx();
    void NotAccept(int ipop);
    void WalkerMatch(int jpop);
    void WalkerCollect(int nsons);
    void Normalization();
    void PrintDistr(const string& name_file);
    void OBDM_Calc( ParamModel& param_model, const Configuration& xaux, WaveFunction & wave_func, MomentDistr& moment_distr, const CorFunParam&  mom_distr_param);

    ~OBDM();
};


#endif
