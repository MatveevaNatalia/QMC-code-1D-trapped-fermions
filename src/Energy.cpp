
#include "Energy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

//void Energy_calc(double **xaux, double **FT, double *Etot, double *Epot_trap, double *Ekin,  int ncomp, int np, double aB, double a, double width)
void Energy_calc(double **xaux, double **FT, Energy& energy,  int ncomp, int np, double aB, double a, double width)
{

double Ekin_J1 = 0, Ekin_J2 = 0, Ekin_G = 0, Ekin_GJ1 = 0, Ekin_GJ2 = 0, Ekin_J1J2 = 0, Epot = 0;
double xi, xk;
double **FG, **FJ1, **FJ2;
double sum = 0;

FG = new double*[ncomp];
for(int i = 0; i < ncomp; i++) FG[i] = new double[np];

FJ1 = new double*[ncomp];
for(int i = 0; i < ncomp; i++) FJ1[i] = new double[np];

FJ2 = new double*[ncomp];
for(int i = 0; i < ncomp; i++) FJ2[i] = new double[np];

for(int alpha = 0; alpha < ncomp; alpha++)
{
	for(int i = 0; i < np; i++)
	{
		FG[alpha][i] = 0.0;
		FJ1[alpha][i] = 0.0;
		FJ2[alpha][i] = 0.0;
	}
}

//Forces (WITHOUT factor of 2)

for(int gamma = 0; gamma < ncomp; gamma++)
{
	for(int k = 0; k < np; k++)
	{
		xk = xaux[gamma][k];

		FG[gamma][k] = -2.0*width*xk;

//               -----------------------

		for(int i = 0; i < np; i++)
		{
			if(i != k) {xi = xaux[gamma][i]; FJ1[gamma][k] = FJ1[gamma][k] + (fabs(xk-xi)-aB)*(xk-xi)/(fabs(fabs(xk-xi)-aB)*fabs(fabs(xk-xi)-aB)*fabs(xk-xi));}
		}

//               -----------------------

		for(int alpha = 0; alpha < ncomp; alpha++)
		{
			if(alpha != gamma)
			{
				for(int i  = 0; i < np; i++)
				{
					xi = xaux[alpha][i];
					FJ2[gamma][k] = FJ2[gamma][k] + (fabs(xk-xi)-a)*(xk-xi)/(fabs(fabs(xk-xi)-a)*fabs(fabs(xk-xi)-a)*fabs(xk-xi));
				}
			}
		}		

	FT[gamma][k] = 2.0*(FG[gamma][k] + FJ1[gamma][k] + FJ2[gamma][k]);
	
//	cout<<"FG: "<<"\n";
//	cout<<"gamma= "<<gamma<<" k= "<<k<<" "<<FG[gamma][k]<<" "<< FJ1[gamma][k]<<" "<<FJ2[gamma][k]<<" "<<FT[gamma][k]<<"\n";	

	}	
}


/////////////////////////////////////////////////////////////////////////////////////////////

for(int gamma = 0; gamma < ncomp; gamma++)
{
	for(int k = 0; k < np; k++)
	{

		xk = xaux[gamma][k];

		Ekin_J1 = Ekin_J1 - 0.5*FJ1[gamma][k]*FJ1[gamma][k];

		for(int i = 0; i < np; i++)
		{
			if( i!= k){xi = xaux[gamma][i]; Ekin_J1 = Ekin_J1 + 0.5/(fabs(fabs(xk-xi)-aB)*fabs(fabs(xk-xi)-aB));}
		}
	

		Ekin_J2 = Ekin_J2 - 0.5 * FJ2[gamma][k]*FJ2[gamma][k];

		for(int alpha = 0; alpha < ncomp; alpha++)
		{
			if(alpha != gamma)
			{
				for(int i = 0; i < np; i++)
				{
					xi = xaux[alpha][i];	

					Ekin_J2 = Ekin_J2 + 0.5 / (fabs(fabs(xk-xi)-a)*fabs(fabs(xk-xi)-a));			
				}

			}
		} 


		Ekin_G = Ekin_G + 0.5*(2.0*width - FG[gamma][k]*FG[gamma][k]);

		Ekin_GJ1 = Ekin_GJ1 - FG[gamma][k]*FJ1[gamma][k];

		Ekin_GJ2 = Ekin_GJ2 - FG[gamma][k]*FJ2[gamma][k];	  

		Ekin_J1J2 = Ekin_J1J2 - FJ1[gamma][k]*FJ2[gamma][k];

		Epot = Epot + 0.5*xk*xk;
	}
}

energy.tot = Ekin_J1 + Ekin_J2 + Ekin_G + Ekin_GJ1 + Ekin_GJ2 + Ekin_J1J2 + Epot;


for(int gamma = 0; gamma < ncomp; gamma++)
{
	for(int k = 0; k < np; k++)
	{
		sum = sum + 0.25*FT[gamma][k]*FT[gamma][k];  
	}
}

energy.kin = 0.5*sum;
//*Epot_trap = Epot;

double x1, x2;
x1 = xaux[0][0];
x2 = xaux[1][0];
//It is not Epot_trap, it is supposed to be the interaction energy
//*Epot_trap = -fabs(x1-x2)/((fabs(x1-x2)-a)*(fabs(x1-x2)-a));
energy.pot = Epot;

//cout<<"Ekin_J1 = "<<Ekin_J1<<" Ekin_J2= "<<Ekin_J2<<" Ekin_G= "<<Ekin_G<<" Ekin_GJ1= "<< Ekin_GJ1<<" Ekin_GJ2= "<<Ekin_GJ2 <<" Ekin_J1J2= "<<Ekin_J1J2<<" Epot_trap = "<<Epot<<"\n";

//cout<< "Etotal= "<<*Eest1<<"\n";


for(int ic = 0; ic < ncomp; ic++)
delete [] FG[ic];
delete [] FG;

for(int ic = 0; ic < ncomp; ic++)
delete [] FJ1[ic];
delete [] FJ1;

for(int ic = 0; ic < ncomp; ic++)
delete [] FJ2[ic];
delete [] FJ2;
	
}
