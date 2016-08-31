#ifndef MOMDISTR_H
#define MOMDISTR_H

#include "qmc.h"
#include <string>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

class MomentDistr{
private:
    int num_points, step;
    void setzero(double * x){
        for(int i = 0; i < num_points; i++) x[i] = 0;
    }
public:
    double *dnkup, *dnkupa, **dnkuplocal;

    MomentDistr(const CorFunParam& mom_distr);

    void SetZero();
    void SetZeroAx();

    void NotAccept(int ipop);

    void WalkerMatch(int jpop);

    void WalkerCollect(int nsons);

    void Normalization(int np, int nkuppt);

    void PrintDistr(const string& name_file);

    ~MomentDistr();
};



#endif // MOMDISTR_H
