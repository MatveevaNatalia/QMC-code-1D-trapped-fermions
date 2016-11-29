#ifndef WAVE_FUN_H
#define WAVE_FUN_H

#include "qmc.h"
#include "Locals.h"
#include <vector>

using namespace std;

class WaveFunction{
private:
    double metrop, auxil;
    vector<double> oldPage, newPage;
public:
    double GetOld(int ipop){
        return oldPage[ipop];
    }
    double GetMetrop(){
        return metrop;
    }
    void Calc(const ParamModel& param_model, Locals& coordinates);
    void Accept(){
        auxil = metrop;
    }
    void NotAccept(int ipop){
        auxil = oldPage[ipop];
    }

    void WalkerMatch(){
        newPage.push_back(auxil);
    }
    double WaveFunction_MC(const ParamModel& param_model, const Configuration& xaux, int ipmac, double xm, int ncomp_MC);

    void PageSwap();

};


#endif
