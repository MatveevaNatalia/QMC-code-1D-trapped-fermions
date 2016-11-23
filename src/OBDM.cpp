
#include "OBDM.h"
#include "Wave_fun.h"
#include "Algorithm.h"

using namespace std;

OBDM::OBDM(const CorFunParam& obdm){
    num_points = obdm.num_points;
    step = obdm.step;
    fr = new double[num_points+1];
    fra = new double[num_points+1];
    nfr = new double[num_points+1];
    nfra = new double[num_points+1];

    frlocal = new double*[num_points+1];
    for(int i = 0; i < (num_points+1); i++) frlocal[i] = new double[dmnpop];

    nfrlocal = new double*[num_points+1];
    for(int i = 0; i < (num_points+1); i++) nfrlocal[i] = new double[dmnpop];
}

void OBDM::SetZero(){
    setzero(fr);
    setzero(nfr);
}

void OBDM::SetZeroAx(){
    setzero(fra);
    setzero(nfra);
}

void OBDM::NotAccept(int ipop){
    for(int i = 1; i < (num_points+1); i++)
    {
        fra[i] = frlocal[i][ipop];
        nfra[i] = nfrlocal[i][ipop];
    }
}

void OBDM::WalkerMatch(int jpop){
    for(int i = 1; i < (num_points+1); i++)
    {
        frlocal[i][jpop] = fra[i];
        nfrlocal[i][jpop] = nfra[i];
    }
}

void OBDM::WalkerCollect(int nsons){
    for(int i = 1; i < (num_points+1); i++)
    {

        fr[i] = fr[i] + nsons * fra[i];
        nfr[i] = nfr[i] + nsons * nfra[i];
    }
}

void OBDM::Normalization(){
    for (int ir = 1; ir < (num_points+1); ir++)
        if(nfr[ir] > 0) fr[ir] = fr[ir]/float(nfr[ir]);
}

void OBDM::PrintDistr(const string& name_file)
{
    double coord_bin;
    fstream outfile(name_file, fstream::out| ios::app );
    for(int i = 1; i < (num_points+1);i++)
    {
        coord_bin = float(i) * step + step/2.0;
        outfile<<coord_bin<<" "<<fr[i]<<"\n";
   }
    outfile.close();
}


OBDM::~OBDM(){
    delete [] fr;
    delete [] fra;
    delete [] nfr;
    delete [] nfra;

    for(int i = 0; i < (num_points+1); i++)
        delete [] frlocal[i];
    delete [] frlocal;

    for(int i = 0; i < (num_points+1); i++)
        delete [] nfrlocal[i];
    delete [] nfrlocal;
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
            fra[idr] = fra[idr] + quoc;
            nfra[idr] = nfra[idr] + 1;
        }

        for(int ik = 0; ik < mom_distr_param.num_points; ik++)
        {
            prod = x1ax * float(ik)*mom_distr_param.step;
            moment_distr.dnkupa[ik] += cos(prod) * quoc;
        }
    }

}



