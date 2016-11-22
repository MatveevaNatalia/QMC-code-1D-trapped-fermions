#include "Algorithm.h"

using namespace std;

void MetropolisDif(int ipop, ParamModel& param_model, double PsiTotal, double **flocal, Locals& coordinates, Locals& force, int& accepta, int& nprova, double& fvella, int ntemps, int in, int i_VMC)
{
    double fdif, QQ, DDF1, DDF2, DDS, dte = param_model.alfa/4;
    double force_old, force_new, x_old, x_new;
    if(ntemps == 1)
    {
        fvella = 0.0;
        QQ = 0.0;
    }
    else
    {
        fvella = flocal[ipop][in];
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
        fdif = 2.0 * (PsiTotal - fvella); //For_VMC_calculations
    if(i_VMC == 0)
        fdif = 2.0 * (PsiTotal - fvella) + QQ; //SVMC and DMC, QQ because of drift jump

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


void BranchingCalc(int *nsons, int accepta, int ntemps,int nacc, int nprova, double dte, int icrit, double E, double enew, double eold, long *kkk, int npop)
{
    int icrmin, icrmax, nsaux;
    double redu, ampi, accrate, rsons, dteff;
    icrmin = icrit - icrit/10.;
    icrmax = icrit + icrit/10.;
    redu = 0.5 * float ( icrmin + icrmax ) / float ( icrmax );
    ampi = 0.5 * float ( icrmin + icrmax ) / float ( icrmin );

    *nsons = 1;

    if(accepta == 1)
    {
        if(ntemps > 1)
        {
            accrate = float(nacc)/float(nprova);

            dteff = dte * accrate;
            if(icrit > 1)
            {
                rsons = exp(2.0*dteff * (E - 0.5 * (enew + eold)));
                *nsons = int(rsons + ran2(kkk));
            }
            else {*nsons = 1;}

            if(nsons > 0 && icrit != 1)
            {
                if(npop > icrmax)
                {
                    nsaux = *nsons * redu + ran2(kkk);
                    *nsons = nsaux;

                }
                if(npop < icrmin)
                {
                    nsaux = *nsons * ampi + ran2(kkk);
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

