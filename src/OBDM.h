#ifndef OBDM_H
#define OBDM_H

#include "qmc.h"

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

    ~OBDM();
};

void OBDM1D_11(double, int, int, double, double, double, double *, double **, double *, double *, double, double *, int, double, long *, int, double);
void OBDM1D_22(double, int, int, double, double, double, double *, double **, double *, double *, double, double *, int, double, long *, int, double);

#endif
