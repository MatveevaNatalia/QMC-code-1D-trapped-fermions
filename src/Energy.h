#ifndef ENERGY_H
#define ENERGY_H
#include "qmc.h"
#include "Locals.h"
using namespace std;

class EnergyState
{
public:
    double tot, pot, kin;
    EnergyState(){
        SetZero();
    }
    void SetZero();
};


class Energy{
private:
    EnergyState average, auxil, metrop, mean, average_part;
    double average_walk;
    vector<EnergyState> oldPage, newPage;
public:
    Energy();
    void SetZeroAverage(){
        average.SetZero();
    }
    void SetZeroMetrop(){
        metrop.SetZero();
    }
    void SetZeroMean(){
        mean.SetZero();
    }
    double GetAverageWalker(){
        return average_walk;
    }

    double GetNewEnergy(){
        return auxil.tot;
    }
    double GetOldEnergy(int i){
        return oldPage[i].tot;
    }
    void Calc(Locals& coord, Configuration& FMT, const ParamModel& param_model);
    void Accept();
    void NotAccept(int ipop);
    void WalkerMatch();
    void WalkerCollect(int nsons);
    void Normalization(const ParamModel& param_model, int jpop);
    void Average();
    void Print(int nwrite, int npopmean, int ntemps, const string & outDir);
    void PageSwap();

};







/*class EnergyVector
{
public:
    double tot, pot, kin;
    void SetZero();
};

class Energy{
private:
    EnergyVector eav, epar, emean, enew, eold, emetrop;
    double ewalk;
    EnergyVector **elocal;
public:
    Energy();
    void SetZeroAverage(){
        eav.SetZero();
    }
    void SetZeroMetrop(){
        emetrop.SetZero();
    }
    void SetZeroMean(){
        emean.SetZero();
    }
    double GetWalkerEnergy(){
        return ewalk;
    }
    double GetNewEnergy(){
        return enew.tot;
    }
    double GetOldEnergy(){
        return eold.tot;
    }
    double GetTotalEnergy(int ipop, int in){
        return elocal[ipop][in].tot;
    }

    void SetOldConf(int ipop, int in);
    void Calc(Locals& coord, Configuration& FMT, const ParamModel& param_model);
    void Accept();
    void NotAccept();
    void WalkerMatch(int jpop, int io);
    void WalkerCollect(int nsons){
        eav.tot = eav.tot + enew.tot * nsons;
    }
    void Normalization(const ParamModel& param_model, int jpop);
    void Average();
    void Print(int nwrite, int npopmean, int ntemps, const string & outDir);

};*/

double EnergyPartial(double xk, double xi, double scat_length);
double ForcePartial(double xk, double xi, double scat_length);

#endif
