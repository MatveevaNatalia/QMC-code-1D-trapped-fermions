#ifndef DISTRR_H
#define DISTRR_H

#include <cmath>
#include "qmc.h"
#include "Locals.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

class DistributionR{
    int num_points;
    double step;
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

    void PairDistrFirst(const Configuration& xMT, const ParamModel& param_model);    
    void PairDistrCross(const Configuration& xMT, const ParamModel& param_model);
    void DensityFirst(const Configuration& xMT, const ParamModel& param_model);
    void PrintDistr(const string& name_file);

    ~DistributionR();
};

#endif // DISTRR_H
