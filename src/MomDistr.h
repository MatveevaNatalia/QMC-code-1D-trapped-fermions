#ifndef MOMDISTR_H
#define MOMDISTR_H

#include "qmc.h"
#include <string>

using namespace std;

class MomentDistr{
private:
    int stored_num_points;
    void setzero(double * x){
        for(int i = 0; i < stored_num_points; i++) x[i] = 0;
    }
public:
    double *dnkup, *dnkupa, **dnkuplocal;

    MomentDistr(int num_points);

    void SetZero();
    void SetZeroAx();

    void NotAccept(int ipop);

    void WalkerMatch(int jpop);

    void WalkerCollect(int nsons);

    void Normalization(int np, int nkuppt);

    ~MomentDistr();
};



#endif // MOMDISTR_H
