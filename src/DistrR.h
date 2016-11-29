#ifndef DISTRR_H
#define DISTRR_H

#include <cmath>
#include <vector>
#include "qmc.h"
#include "Locals.h"
#include "ConfigDistr.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

class DistributionR{
private:
    ConfigDistr average, auxil;
    vector<ConfigDistr> oldPage, newPage;
    int num_points;
    double step;
public:
    DistributionR(const CorFunParam& param);
    void SetZero();
    void SetZeroAx();
    void NotAccept(int ipop);
    void WalkerMatch();
    void WalkerCollect(int nsons);
    void NormalizationGR(int nsons_total, const ParamModel& param);
    void NormalizationNR(int nsons_total);

    void PairDistrFirst(const Configuration& xMT, const ParamModel& param_model);
    void PairDistrCross(const Configuration& xMT, const ParamModel& param_model);
    void DensityFirst(const Configuration& xMT, const ParamModel& param_model);
    void PrintDistr(const string& name_file);
    void PageSwap();
};

#endif // DISTRR_H
