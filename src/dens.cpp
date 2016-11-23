
#include "dens.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void DensityDistribution_calc(double **xMT, double *nraMT,  int mgr, double dr, int np)
{
    double Lmax, x1;
    int igr1;
    Lmax = dr*mgr;
    for(int i = 0; i < np; i++)
    {
        x1 = fabs(xMT[0][i]);

        if(x1 < Lmax)
        {
            igr1 = x1/dr;
            nraMT[igr1] ++;
        }
    }

}




