#include "Wave_fun.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>


using namespace std;

double WaveFunction(const ParamModel& param_model, Locals& coordinates)
{

    double PsiG = 0, PsiJ1 = 0, PsiJ2 = 0;
    double xi, xj;

    for(int alpha = 0; alpha < param_model.num_comp; alpha++)
    {
        for(int i = 0; i < param_model.num_part; i++)
        {
            xi = coordinates.metrop[alpha][i];

            PsiG = PsiG -param_model.width*xi*xi;

            for(int j = i+1; j < param_model.num_part; j++)
            {
                xj = coordinates.metrop[alpha][j];
                PsiJ1 = PsiJ1 + log(fabs(fabs(xi - xj)-param_model.scat_lenght_bos));
            }
        }
    }

    for(int alpha = 0; alpha < param_model.num_comp; alpha++)
    {
        for(int beta = alpha+1; beta < param_model.num_comp; beta++)
        {
            for(int i = 0; i < param_model.num_part; i++)
            {
                for(int j = 0; j < param_model.num_part; j++)
                {
                    xi = coordinates.metrop[alpha][i];
                    xj = coordinates.metrop[beta][j];
                    PsiJ2 = PsiJ2 + log(fabs(fabs(xi-xj)-param_model.scat_lenght));
                }
            }
        }
    }

    double PsiTotal = PsiG + PsiJ1 + PsiJ2;

    return PsiTotal;

    //cout<<" LogPsiG= "<<PsiG<<" LogPsiJ1= "<<PsiJ1<<" LogPsiJ2= "<<PsiJ2<<"\n";

}

//****************************************************************


double WaveFunction_MC(int ncomp, int np, double aB, double a, double width, double **xaux, double *Psi_MC, int ipmac, double xm, int ncomp_MC) 
{

    double xi, xj;
    double PsiG = 0, PsiJ1 = 0, PsiJ2 = 0;

    //cout<<" From MC wf: "<<" ipmac= "<<ipmac<<"xm= "<<xm<<"\n";

    for(int alpha = 0; alpha < ncomp; alpha++)
    {
        for(int i = 0; i < np; i++)
        {
            if(alpha == ncomp_MC && i == ipmac) {xi = xm;}
            else {xi = xaux[alpha][i];}

            PsiG = PsiG -width*xi*xi;

            for(int j = i+1; j < np; j++)
            {
                if(alpha == ncomp_MC && j == ipmac) {xj = xm;}
                else {xj = xaux[alpha][j];}

                //			cout<<" From MC wf: "<< "xi= "<<xi<<" xj= "<<xj<<"\n";

                PsiJ1 = PsiJ1 + log(fabs(fabs(xi - xj)-aB));
            }
        }
    }



    for(int alpha = 0; alpha < ncomp; alpha++)
    {
        for(int beta = alpha+1; beta < ncomp; beta++)
        {
            for(int i = 0; i < np; i++)
            {
                for(int j = 0; j < np; j++)
                {
                    if(alpha == ncomp_MC && i == ipmac) {xi = xm;}
                    else {xi = xaux[alpha][i];}

                    if(beta == ncomp_MC && j == ipmac) {xj = xm;}
                    else {xj = xaux[beta][j];}

                    PsiJ2 = PsiJ2 + log(fabs(fabs(xi-xj)-a));
                }
            }


        }
    }

    *Psi_MC = PsiG + PsiJ1 + PsiJ2;

    //cout<<" LogPsiG= "<<PsiG<<" LogPsiJ1= "<<PsiJ1<<" LogPsiJ2= "<<PsiJ2<<"\n";

}


