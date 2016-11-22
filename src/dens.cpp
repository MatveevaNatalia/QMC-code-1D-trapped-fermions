
#include "dens.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void DensityDistribution_calc(double **xMT, double *nraMT_11,  int mgr, double dr, int np)
{

    double Lmax, x1;
    int igr1;
    Lmax = dr*mgr;

    for(int i = 0; i < np; i++)// PairDistribution::CalculateDensity
    {
        x1 = fabs(xMT[0][i]);

        if(x1 < Lmax)
        {
            igr1 = x1/dr;
            nraMT_11[igr1] = nraMT_11[igr1] + 1;
        }

    }

}




