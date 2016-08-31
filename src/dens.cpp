
#include "dens.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void DensityDistribution_calc(double **xMT, double *nraMT_11, double *nraMT_22, int mgr, double dr, int ncomp, int np)
{
    // In this form only for two-component case, I do not do the xact counting how many times
    // the distance is less than Lmax
    double Lmax, x1, x2;
    int igr1, igr2;
    Lmax = dr*mgr;


    /*if(ncomp == 2)
{
    for(int i = 0; i < np; i++)
    {
        x1 = fabs(xMT[0][i]);
        x2 = fabs(xMT[1][i]);

       // cout<<"x1= "<<x1<<" "<<" Lmax= "<<Lmax<<endl;
        if(x1 < Lmax)
        {
            igr1 = x1/dr;
            nraMT_11[igr1] = nraMT_11[igr1] + 1;
        }

        if(x2 < Lmax)
        {
            igr2 = x2/dr;
            nraMT_22[igr2] = nraMT_22[igr2] + 1;
        }

    }
}*/

    //if(ncomp == 1)
    //{
    for(int i = 0; i < np; i++)// PairDistribution::CalculateDensity
    {
        x1 = fabs(xMT[0][i]);

        if(x1 < Lmax)
        {
            igr1 = x1/dr;
            nraMT_11[igr1] = nraMT_11[igr1] + 1;
        }

    }
    //}
}


/*void DensitySecond(double **xMT, double *nraMT, int mgr, double dr, int np)
{
    double Lmax, deltax;
    int bin_num;
    Lmax = dr*mgr;

    for(int i = 0; i < np; i++)
    {
        deltax = fabs(xMT[0][i]);

        if(deltax < Lmax)
        {
            bin_num = deltax/dr;
            nraMT[bin_num] += 1;
        }
    }
}*/

