
#include "g2D.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void PairDistribution_calc(double **xMT, double *graMT_11, double *graMT_22, double *graMT_12, int mgr, double dr, int ncomp, int np)
{

    double Lmax, x1, x2, x12;
    int igr1, igr12;
    Lmax = dr*mgr;

    for(int i = 0; i < np; i++)
    {
        for(int j = i+1; j < np; j++)
        {
            x1 = fabs(xMT[0][i] - xMT[0][j]);
            x2 = fabs(xMT[1][i] - xMT[1][j]);

            if(x1 < Lmax)
            {
                igr1 = x1/dr;
                graMT_11[igr1] = graMT_11[igr1] + 1;
            }

        }
    }

    for(int icomp = 1; icomp < ncomp; icomp++)
    {
        for(int i = 0; i < np; i++)
        {
            for(int j = 0; j < np; j++)
            {
                x12 = fabs(xMT[0][i] - xMT[icomp][j]);
                if(x12 < Lmax)
                {
                    igr12 = x12/dr;
                    graMT_12[igr12] = graMT_12[igr12] + 1;
                }
            }
        }
    }
}










