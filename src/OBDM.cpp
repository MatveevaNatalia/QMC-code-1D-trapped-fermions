#include "OBDM.h"
#include "Wave_fun.h"
#include "Algorithm.h"

using namespace std;

OBDM::OBDM(const CorFunParam& param):
    average(param.num_points),
    auxil(param.num_points),
    average_norm(param.num_points),
    auxil_norm(param.num_points),
    oldPage(dmnpop, ConfigDistr(param.num_points)),
    oldPage_norm(dmnpop, ConfigDistr(param.num_points))
{
    num_points = param.num_points;
    step = param.step;
}

void OBDM::SetZero(){
    average.setzero();
    average_norm.setzero();
}

void OBDM::SetZeroAx(){
    auxil.setzero();
    auxil_norm.setzero();
}

void OBDM::NotAccept(int ipop){
    auxil = oldPage[ipop];
    auxil_norm = oldPage_norm[ipop];
}

void OBDM::WalkerMatch(){
    newPage.push_back(auxil);
    newPage_norm.push_back(auxil_norm);
}

void OBDM::WalkerCollect(int nsons){
    for(int i = 1; i < (num_points+1); i++){
        average[i] += nsons * auxil[i];
        average_norm[i] += nsons * auxil_norm[i];
    }
}

void OBDM::Normalization(){
    for (int i = 1; i < (num_points+1); i++){
        if(average_norm[i] > 0)
            average[i] = average[i]/float(average_norm[i]);
    }
}

void OBDM::PageSwap(){
     oldPage = newPage;
     newPage.clear();
     oldPage_norm = newPage_norm;
     newPage_norm.clear();

}

void OBDM::PrintDistr(const string& name_file)
{
    double coord_bin;
    fstream outfile(name_file, fstream::out| ios::app );
    for(int i = 1; i < (num_points+1);i++)
    {
        coord_bin = float(i) * step + step/2.0;
        outfile<<coord_bin<<" "<<average[i]<<"\n";
   }
    outfile.close();
}

void OBDM::OBDM_Calc( ParamModel& param_model, const Configuration& xaux, WaveFunction & wave_func, MomentDistr& moment_distr, const CorFunParam&  mom_distr_param)
{
    double x1m, x1ax, r1ax;
    double Psi_MC = 0.0;
    double quoc, prod;
    int ipmac, idr, num_MC = 0;
    double xi, xj, Lmax, Lmax_ro;

    int sign_old = 1, sign_MC;
    double dx;
    int nc, np;

    Lmax_ro = step * num_points;
    np = param_model.num_part;
    nc = param_model.num_comp;
    Lmax = param_model.Lmax;

    for(int alpha = 0; alpha < nc; alpha++)
    {
        for(int i = 0; i < np; i++)
        {
            for(int j = i+1; j < np; j++)
            {
                xi = xaux.GetParticleComp(alpha,i);
                xj = xaux.GetParticleComp(alpha,j);

                dx = xi - xj;
                sign_old = sign_old * int(dx/fabs(dx));
            }
        }
    }

    for(int nmac = 1; nmac <= 100; nmac++)
    {
        sign_MC = 1;

        //The coefficients here are the optimized ones
        //for the algorithm

        x1m = -5.0*Lmax + 10.0*ran2(&(param_model.seed)) * Lmax;
        ipmac = int(ran2(&(param_model.seed)) * np);

        x1ax = x1m - xaux.GetParticleComp(num_MC, ipmac);

        r1ax = fabs(x1ax);

        Psi_MC = wave_func.WaveFunction_MC(param_model, xaux, ipmac, x1m, num_MC);

        for(int alpha = 0; alpha < nc; alpha++)
        {
            for(int i = 0; i < np; i++)
            {
                for(int j = i+1; j < np; j++)
                {
                    if(alpha == 0 && i == ipmac) {xi = x1m;}
                    else{xi = xaux.GetParticleComp(alpha, i);}

                    if(alpha == 0 && j == ipmac) {xj = x1m;}
                    else{xj = xaux.GetParticleComp(alpha, j);}

                    dx = xi - xj;
                    sign_MC = sign_MC * int(dx/fabs(dx));
                }
            }
        }

        quoc = exp(Psi_MC - wave_func.GetMetrop())*sign_MC/sign_old;

        if(r1ax < Lmax_ro)
        {
            idr = int(fabs(r1ax)/step);            
            auxil[idr] += quoc;
            auxil_norm[idr] ++;
        }

        for(int ik = 0; ik < mom_distr_param.num_points; ik++)
        {
            prod = x1ax * float(ik)*mom_distr_param.step;
            moment_distr.auxil[ik] += cos(prod) * quoc;
        }
    }
}




