#include "Algorithm.h"
using namespace std;

void MetropolisDif(int ipop, ParamModel& param_model, WaveFunction& wave_func, Locals& coordinates, Locals& force, int& accepta, int& nprova, int ntemps, int in, int i_VMC)
{
    double  fdif, QQ, DDF1, DDF2, DDS, dte = param_model.alfa/4;
    double force_old, force_new, x_old, x_new;
    if(ntemps == 1)
    {
        wave_func.SetOldZero();
        QQ = 0.0;
    }
    else
    {        
        wave_func.SetOld(ipop,in);
        QQ = 0.0;
        for(int inc = 0; inc < param_model.num_comp; inc++)
        {
            for(int ip = 0; ip < param_model.num_part; ip ++)
            {

                force_old = force.oldPage[ipop].GetParticleComp(inc, ip);
                force_new = force.metrop.GetParticleComp(inc, ip);
                x_old = coordinates.metrop.GetParticleComp(inc, ip);
                x_new = coordinates.oldPage[ipop].GetParticleComp(inc, ip);

                DDF1 = force_old + force_new;
                DDF2 = force_old - force_new;
                DDS = x_new - x_old;
                QQ = QQ + 0.5 * DDF1 * (0.5 * dte * DDF2 - DDS);
            }
        }
    }

    if(i_VMC == 1)
        fdif = 2.0 * (wave_func.GetMetrop() - wave_func.GetOld()); //For_VMC_calculations
    if(i_VMC == 0)
        fdif = 2.0 * (wave_func.GetMetrop() - wave_func.GetOld()) + QQ; //SVMC and DMC, QQ because of drift jump

    accepta = 1;
    if(ntemps > 1)
    {
        nprova = nprova + 1;
        if (fdif >= 0)
            accepta = 1;
        else
        {
            accepta = 0;
            if(ran2(&(param_model.seed))< exp(fdif))
                accepta = 1;
        }
    }
}


//****************************************************************************//


//BranchingCalc(&nsons, accepta, ntemps, nacc, nprova, param_model.alfa/4.0, nwalk_mean, ewalk, enew.tot, eold.tot, &param_model.seed, param_model.num_walk);

void BranchingCalc(ParamModel& param_model, Energy& energy, int *nsons, int accepta, int ntemps,int nacc, int nprova)
//void BranchingCalc(int *nsons, int accepta, int ntemps,int nacc, int nprova, double dte, int nwa, double E, double enew, double eold, long *kkk, int npop)
{
    int nwalk_min, nwalk_max, nwalk_mean, nsaux, npop;
    double redu, ampi, accrate, rsons, dte, dteff;
    double energy_walk, enew, eold;

    npop = param_model.num_walk;
    nwalk_mean = param_model.nwalk_mean;
    nwalk_min = nwalk_mean - nwalk_mean/10.;
    nwalk_max = nwalk_mean + nwalk_mean/10.;
    redu = 0.5 * float ( nwalk_min + nwalk_max ) / float ( nwalk_max );
    ampi = 0.5 * float ( nwalk_min + nwalk_max ) / float ( nwalk_min );

    energy_walk = energy.GetWalkerEnergy();
    enew = energy.GetNewEnergy();
    eold = energy.GetOldEnergy();

    dte = param_model.alfa/4.0;
    *nsons = 1;

    if(accepta == 1)
    {
        if(ntemps > 1)
        {
            accrate = float(nacc)/float(nprova);

            dteff = dte * accrate;
            if(nwalk_mean > 1)
            {
                rsons = exp(2.0*dteff * (energy_walk - 0.5 * (enew + eold)));
                *nsons = int(rsons + ran2(&(param_model.seed)));
            }
            else {*nsons = 1;}

            if(nsons > 0 && nwalk_mean != 1)
            {
                if(npop > nwalk_max)
                {
                    nsaux = *nsons * redu + ran2(&(param_model.seed));
                    *nsons = nsaux;

                }
                if(npop < nwalk_min)
                {
                    nsaux = *nsons * ampi + ran2(&(param_model.seed));
                    *nsons = nsaux;
                }
            }
        }
    }
}

//***************************************************************************//

void Gauss1D(double * x, double alfa, long *kkk)
{
    double fi, zzz, r;
    fi = 6.283185307 * ran2(kkk);
    zzz = ran2(kkk);
    while( zzz == 0.0) zzz = ran2(kkk);
    r = sqrt(-alfa * log(zzz));
    *x = r * cos(fi);   
}

//**********************************************//
// Taken from Numerical receipts
//**********************************************//
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
float ran2(long *idum) {
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

