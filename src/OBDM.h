#ifndef OBDM_H
#define OBDM_H

#include "qmc.h"

#include <string>

using namespace std;

class OBDM{
private:
    int stored_num_points;
    void setzero(double * x){
        for(int i = 0; i < stored_num_points+1; i++) x[i] = 0;
    }
public:
    double *fr, *fra, *nfr, *nfra, **frlocal, **nfrlocal;
    OBDM(int num_points);

    void SetZero();

    void SetZeroAx();

    void NotAccept(int ipop);

    void WalkerMatch(int jpop);

    void WalkerCollect(int nsons);

    void Normalization();

    ~OBDM();
};

void OBDM1D_11(double, int, int, double, double, double, double *, double **, double *, double *, double, double *, int, double, long *, int, double);
void OBDM1D_22(double, int, int, double, double, double, double *, double **, double *, double *, double, double *, int, double, long *, int, double);

#endif
