
#include "Wave_fun.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

double WaveFunction(int ncomp, int np, double aB, double a, double width, double **xaux, double *PsiTotal) // Returns the logarithm of the total trial wave function
{

double PsiG = 0, PsiJ1 = 0, PsiJ2 = 0;


for(int alpha = 0; alpha < ncomp; alpha++)
{
	for(int i = 0; i < np; i++)
	{
		PsiG = PsiG -width*xaux[alpha][i]*xaux[alpha][i];

		for(int j = i+1; j < np; j++)
		{
			PsiJ1 = PsiJ1 + log(fabs(fabs(xaux[alpha][i] - xaux[alpha][j])-aB));
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
				PsiJ2 = PsiJ2 + log(fabs(fabs(xaux[alpha][i]-xaux[beta][j])-a));
			}
		}


	}
}

*PsiTotal = PsiG + PsiJ1 + PsiJ2;

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
