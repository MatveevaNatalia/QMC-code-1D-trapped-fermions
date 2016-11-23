#ifndef WAVE_FUN_H
#define WAVE_FUN_H

#include "qmc.h"
#include "Locals.h"

using namespace std;

class WaveFunction{
private:
    double fmetrop, fnew, fold;
public:
    double **flocal;
    WaveFunction();
    double GetOld(){
        return fold;
    }
    double GetMetrop(){
        return fmetrop;
    }
    void Calc(const ParamModel& param_model, Locals& coordinates);
    void Accept(){
        fnew = fmetrop;
    }
    void NotAccept(){
        fnew = fold;
    }
    void WalkerMatch(int jpop, int io){
        flocal[jpop][io] = fnew;
    }
    void SetOldZero(){
        fold = 0.;
    }
    void SetOld(int ipop, int in){
        fold = flocal[ipop][in];
    }
    double WaveFunction_MC(const ParamModel& param_model, const Configuration& xaux, int ipmac, double xm, int ncomp_MC);
    ~WaveFunction();
};

#endif
