
#include "g2D.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void PairDistribution_calc(double **xMT, double *graMT_11, double *graMT_22, double *graMT_12, int mgr, double dr, int ncomp, int np)
{
// I do not do the xact counting how many times
// the distance is less than Lmax
double Lmax, x1, x2, x12;
int igr1, igr2, igr12;
Lmax = dr*mgr;

for(int i = 0; i < np; i++) // PairDistribution::CalculateSelf
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

//		if(x2 < Lmax)
//		{
//			igr2 = x2/dr;
//			graMT_22[igr2] = graMT_22[igr2] + 1; 
//		}
	}
}

for(int icomp = 1; icomp < ncomp; icomp++) // PairDistribution::CalculateCross
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


void CalculateFirst(double **xMT, double *graMT_11, int mgr, double dr, int ncomp, int np)
{
    double Lmax, deltax;
    int bin_number;
    Lmax = dr*mgr;

    for(int i = 0; i < np; i++)
    {
        for(int j = i+1; j < np; j++)
        {
            deltax = fabs(xMT[0][i] - xMT[0][j]);

            if(deltax < Lmax)
            {
                bin_number = deltax/dr;
                graMT_11[bin_number] = graMT_11[bin_number] + 1;
            }

        }
    }
}

void CalculateSecond(double **xMT, double *graMT_22, int mgr, double dr, int ncomp, int np)
{
    double Lmax, deltax;
    int bin_number;
    Lmax = dr*mgr;

    for(int i = 0; i < np; i++)
    {
        for(int j = i+1; j < np; j++)
        {
            deltax = fabs(xMT[1][i] - xMT[1][j]);

            if(deltax < Lmax)
            {
                bin_number = deltax/dr;
                graMT_22[bin_number] = graMT_22[bin_number] + 1;
            }

        }
    }
}





void CalculateCross(double **xMT, double *graMT_12, int mgr, double dr, int ncomp, int np)
{
    double Lmax, x12;
    int  bin_number;
    Lmax = dr*mgr;
    for(int icomp = 1; icomp < ncomp; icomp++)
    {
        for(int i = 0; i < np; i++)
        {
            for(int j = 0; j < np; j++)
            {
                x12 = fabs(xMT[0][i] - xMT[icomp][j]);
                if(x12 < Lmax)
                {
                    bin_number = x12/dr;
                    graMT_12[bin_number] = graMT_12[bin_number] + 1;
                }
            }
        }
    }
}
