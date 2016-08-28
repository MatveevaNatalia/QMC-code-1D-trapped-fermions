#ifndef DISTRR_H
#define DISTRR_H

#include "qmc.h"
using namespace std;


class DistributionR{
    int stored_num_points;
    void setzero(double * x){
        for(int i = 0; i < stored_num_points+1; i++) x[i] = 0;
    }

public:
    double *dr, *dra, *draMT;
    double ***drlocal;

    DistributionR(int num_points);

    void SetZero();

    void SetZeroAx();

    void Accept();

   void NotAccept(int ipop, int in);

   void WalkerMatch(int jpop, int io);

   void WalkerCollect(int nsons);

   void NormalizationGR(int ngr, int ncomp, int np, double step);
   void NormalizationNR(int ngr, float step);

    ~DistributionR();
};

#endif // DISTRR_H
