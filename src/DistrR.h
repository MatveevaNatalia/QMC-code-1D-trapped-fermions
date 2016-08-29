#ifndef DISTRR_H
#define DISTRR_H

#include <cmath>
#include "qmc.h"
using namespace std;


class DistributionR{
    int num_points;
    int step;
    void setzero(double * x){
        for(int i = 0; i < num_points+1; i++) x[i] = 0;
    }

public:
    double *dr, *dra, *draMT;
    double ***drlocal;

    DistributionR(const CorFunParam& pair_distr);

    void SetZero();

    void SetZeroAx();

    void Accept();

    void NotAccept(int ipop, int in);

    void WalkerMatch(int jpop, int io);

    void WalkerCollect(int nsons);

    void NormalizationGR(int ngr, int ncomp, int np, double step);
    void NormalizationNR(int ngr, float step);

    void PairDistrFirst(double **xMT, const ParamModel& param_model);
    void PairDistrSecond(double **xMT, const ParamModel& param_model);
    void PairDistrCross(double **xMT, const ParamModel& param_model );

    ~DistributionR();
};

#endif // DISTRR_H
