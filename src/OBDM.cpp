
#include "OBDM.h"
#include "Wave_fun.h"
#include "Algorithm.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

OBDM::OBDM(int num_points){
    stored_num_points = num_points;
    fr = new double[num_points+1];
    fra = new double[num_points+1];
    nfr = new double[num_points+1];
    nfra = new double[num_points+1];

    frlocal = new double*[num_points+1];
    for(int i = 0; i < (num_points+1); i++) frlocal[i] = new double[dmnpop];

    nfrlocal = new double*[num_points+1];
    for(int i = 0; i < (num_points+1); i++) nfrlocal[i] = new double[dmnpop];
}

void OBDM::SetZero(){
    setzero(fr);
    setzero(nfr);
}

void OBDM::SetZeroAx(){
    setzero(fra);
    setzero(nfra);
}

void OBDM::NotAccept(int ipop){
    for(int i = 1; i < (stored_num_points+1); i++)
    {
        fra[i] = frlocal[i][ipop];
        nfra[i] = nfrlocal[i][ipop];
    }
}

void OBDM::WalkerMatch(int jpop){
    for(int i = 1; i < (stored_num_points+1); i++)
    {
        frlocal[i][jpop] = fra[i];
        nfrlocal[i][jpop] = nfra[i];
    }
}

void OBDM::WalkerCollect(int nsons){
    for(int i = 1; i < (stored_num_points+1); i++)
    {

        fr[i] = fr[i] + nsons * fra[i];
        nfr[i] = nfr[i] + nsons * nfra[i];
    }
}

void OBDM::Normalization(){
    for (int ir = 1; ir < (stored_num_points+1); ir++)
        if(nfr[ir] > 0) fr[ir] = fr[ir]/float(nfr[ir]);
}

OBDM::~OBDM(){
    delete [] fr;
    delete [] fra;
    delete [] nfr;
    delete [] nfra;

    for(int i = 0; i < (stored_num_points+1); i++)
        delete [] frlocal[i];
    delete [] frlocal;

    for(int i = 0; i < (stored_num_points+1); i++)
        delete [] nfrlocal[i];
    delete [] nfrlocal;
}

void OBDM1D_11(double Lmax, int nc, int np, double width, double aB, double a, double *kr, double **xaux, double *fra, double *nfra, double Psi_old, double *dnkupa, int numks, double dk, long *kkk, int mgr, double dr)
{
double x1m, x1ax, r1ax;
double x1, r1, a0f, a1f, a2f, ulog;
double Psi_MC = 0.0;
double quocsign, quoc, prod;
int ipmac, idr;
double xi, xj, Lmax_ro;

int sign_old = 1, sign_MC;
double dx;
Lmax_ro = dr*mgr;

// To get the sign for the old wf (contributions for all components)

for(int alpha = 0; alpha < nc; alpha++)
{
	for(int i = 0; i < np; i++)
	{
		for(int j = i+1; j < np; j++)
		{
			dx = xaux[alpha][i] - xaux[alpha][j];
			sign_old = sign_old * int(dx/fabs(dx));
		}
	}
}

//cout<<" log(Psi_old)= "<<Psi_old<<" sign_old="<<sign_old<<"\n";


//For the first component
for(int nmac = 1; nmac <= 100; nmac++)     //
{
	sign_MC = 1;

	x1m = -5.0*Lmax + 10.0*ran2(kkk) * Lmax; //Should be in the range from [-Lmax:Lmax] - main thing that it should be the s
	ipmac = int(ran2(kkk) * np); 

	x1ax = x1m - xaux[0][ipmac];

	r1ax = fabs(x1ax);


	WaveFunction_MC(nc, np, aB, a, width, xaux, &Psi_MC, ipmac, x1m, 0);

	for(int alpha = 0; alpha < nc; alpha++)
	{
		for(int i = 0; i < np; i++)
		{	
			for(int j = i+1; j < np; j++)
			{
				if(alpha == 0 && i == ipmac) {xi = x1m;}
				else{xi = xaux[alpha][i];} 	

				if(alpha == 0 && j == ipmac) {xj = x1m;}
				else{xj = xaux[alpha][j];}	

				dx = xi - xj; 
				sign_MC = sign_MC * int(dx/fabs(dx));
			}
		}
	}
	

	quoc = exp(Psi_MC - Psi_old)*sign_MC/sign_old; //may be int type can cause troubles???

	
	if(r1ax < Lmax_ro)  // HERE MUST BE Lmax_ro= dr*mgr
	{
		idr = int(fabs(r1ax)/dr);
		fra[idr] = fra[idr] + quoc;
		nfra[idr] = nfra[idr] + 1;
	}

	for(int ik = 0; ik < numks; ik++)
	{
	
		prod = x1ax * float(ik)*dk;
		dnkupa[ik] = dnkupa[ik] + cos(prod) * quoc;
	}

//		*nkppta = *nkppta + 1;
}

//cout<<" From OBDM_11 "<<nkppta<<" ";
}

///////////////////////////////////////////////////////////////////////////////////////////////////

/*void OBDM1D_22(double Lmax, int nc, int np, double width, double aB, double a, double *kr, double **xaux, double *fra, double *nfra, double Psi_old, double *dnkupa, int numks, double dk, long *kkk, int mgr, double dr)
{
double x1m, x1ax, r1ax;
double x1, r1, a0f, a1f, a2f, ulog;
double Psi_MC = 0.0;
double quocsign, quoc, prod;
int ipmac, idr;
double xi, xj, Lmax_ro;

int sign_old = 1, sign_MC;
double dx;
Lmax_ro = dr*mgr;

// To get the sign for the old wf (contributions for both components)


for(int alpha = 0; alpha < nc; alpha++)
{
	for(int i = 0; i < np; i++)
	{
		for(int j = i+1; j < np; j++)
		{
			dx = xaux[alpha][i] - xaux[alpha][j];
			sign_old = sign_old * int(dx/fabs(dx));
		}
	}
}

//cout<<" log(Psi_old)= "<<Psi_old<<" sign_old="<<sign_old<<"\n";


//For the first component
for(int nmac = 1; nmac <= 100; nmac++)     //
{
	sign_MC = 1;

	x1m = -5.0*Lmax + 10.0*ran2(kkk) * Lmax; //Should be in the range from [-Lmax:Lmax]???
	ipmac = int(ran2(kkk) * np); 

	x1ax = x1m - xaux[1][ipmac]; //ncomp_mac

	r1ax = fabs(x1ax);

//	cout<<" r1ax from OBDM_22 "<<r1ax<<"\n";

 
	WaveFunction_MC(2, np, aB, a, width, xaux, &Psi_MC, ipmac, x1m, 1); //ncomp_mac

	for(int alpha = 0; alpha < nc; alpha++)
	{
		for(int i = 0; i < np; i++)
		{	
			for(int j = i+1; j < np; j++)
			{
				if(alpha == 1 && i == ipmac) {xi = x1m;} //ncomp_mac
				else{xi = xaux[alpha][i];} 	

				if(alpha == 1 && j == ipmac) {xj = x1m;} //ncomp_mac
				else{xj = xaux[alpha][j];}	

				dx = xi - xj; 
				sign_MC = sign_MC * int(dx/fabs(dx));
			}
		}
	}
	

	quoc = exp(Psi_MC - Psi_old)*sign_MC/sign_old; //may be int type can cause troubles???


	if(r1ax < Lmax_ro)
	{
		idr = int(fabs(r1ax)/dr);
		fra[idr] = fra[idr] + quoc;
		nfra[idr] = nfra[idr] + 1;
	
			
	
//	cout<<"nmac= "<<nmac<<" ipmac= "<<ipmac<<" x1m= "<<x1m<<" Psi_MC= "<<Psi_MC<<"\n"; 
//	cout<<" sign_new= "<<sign_MC<<" quoc= "<<quoc<<" idr "<<idr<<"\n";	
	}

	for(int ik = 0; ik < numks; ik++)
	{
	
		prod = x1ax * float(ik)*dk;
		dnkupa[ik] = dnkupa[ik] + cos(prod) * quoc;
	}

}
//	*nkppta = *nkppta + 1;
//	cout<<" from OBDM_22 "<<nkppta<<"\n";
}*/

