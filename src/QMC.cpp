// Both for bosons and fermions,
// h^2/(m a_osc^2) = 1, D = 1/2


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Algorithm.h"
#include "Energy.h"
#include "g2D.h"
#include "dens.h"
#include "OBDM.h"
#include "Statistics.h"
#include "Wave_fun.h"
#include "qmc.h"

using namespace std;
#define BIG  1e30 
const int dmnpop = 300;
const double pi = 3.141592653589793; 

void Energy::SetZero(){
    tot = 0; pot = 0; kin = 0;
};


class DistributionR{
    int stored_num_points;
    void setzero(double * x){
        for(int i = 0; i < stored_num_points+1; i++) x[i] = 0;
    }

public:
    double *dr, *dra, *draMT;
    double ***drlocal;

    DistributionR(int num_points){
        stored_num_points = num_points;
        dr = new double[num_points+1];
        dra = new double[num_points+1];
        draMT = new double[num_points+1];

        drlocal = new double**[num_points+1];
        for(int i = 0; i < (num_points+1); i++){
            drlocal[i] = new double*[dmnpop];
            for(int j = 0; j < dmnpop; j++ ) drlocal[i][j] =  new double[2];
        }
    }

    void SetZero(){
        setzero(dr);
    }

    void SetZeroAx(){
        setzero(dra);
        setzero(draMT);
    }

    void Accept(){
        for(int i = 1; i < (stored_num_points+1); i++)
            dra[i] = draMT[i];
    }

   void NotAccept(int ipop, int in){
       for(int i = 1; i < (stored_num_points+1); i++)
           dra[i] = drlocal[i][ipop][in];
   }

   void WalkerMatch(int jpop, int io){
       for(int i = 1; i < (stored_num_points+1); i++)
           drlocal[i][jpop][io] = dra[i];
   }

   void WalkerCollect(int nsons){
       for(int i = 1; i < (stored_num_points+1); i++)
           dr[i] = dr[i] + nsons * dra[i];
   }

   void NormalizationGR(int ngr, int ncomp, int np, double step){
       double r1;
       for (int ir = 1; ir < (stored_num_points+1); ir++)
       {
           r1 = float(ir) * step;
           if (ngr > 0)
               dr[ir] = dr[ir]/float(ngr);
           else
               dr[ir] = 0.0;
           dr[ir] = dr[ir]/(ncomp*np*step);
       }
   }

   void NormalizationNR(int ngr, float step){
       double r1;
       for (int ir = 1; ir < (stored_num_points+1); ir++)
       {
           r1 = float(ir) * step;
           if (ngr > 0)
               dr[ir] = dr[ir]/float(ngr);
           else
               dr[ir] = 0.0;

           dr[ir] = dr[ir]/(2*step);
       }
   }

    ~DistributionR(){
        delete[] dra;
        delete[] dr;
        delete[] draMT;
        for(int i = 0; i < (stored_num_points+1); i++)
        {
            for (int j = 0; j < dmnpop; j++) delete [] drlocal[i][j];
            delete [] drlocal[i];
        }
        delete [] drlocal;

     }
};

class OBDM{
private:
    int stored_num_points;
    void setzero(double * x){
        for(int i = 0; i < stored_num_points+1; i++) x[i] = 0;
    }
public:
    double *fr, *fra, *nfr, *nfra, **frlocal, **nfrlocal;
    OBDM(int num_points){
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

    void SetZero(){
        setzero(fr);
        setzero(nfr);
    }

    void SetZeroAx(){
        setzero(fra);
        setzero(nfra);
    }

    void NotAccept(int ipop){
        for(int i = 1; i < (stored_num_points+1); i++)
        {
            fra[i] = frlocal[i][ipop];
            nfra[i] = nfrlocal[i][ipop];
        }
    }

    void WalkerMatch(int jpop){
        for(int i = 1; i < (stored_num_points+1); i++)
        {
            frlocal[i][jpop] = fra[i];
            nfrlocal[i][jpop] = nfra[i];
        }
    }

    void WalkerCollect(int nsons){
        for(int i = 1; i < (stored_num_points+1); i++)
        {

            fr[i] = fr[i] + nsons * fra[i];
            nfr[i] = nfr[i] + nsons * nfra[i];
        }
    }

    void Normalization(){
        for (int ir = 1; ir < (stored_num_points+1); ir++)
            if(nfr[ir] > 0) fr[ir] = fr[ir]/float(nfr[ir]);
    }



    ~OBDM(){
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
};


class MomentDistr{
private:
    int stored_num_points;
    void setzero(double * x){
            for(int i = 0; i < stored_num_points; i++) x[i] = 0;
    }
public:
    double *dnkup, *dnkupa, **dnkuplocal;

    MomentDistr(int num_points){
        stored_num_points = num_points;
        dnkup = new double[num_points];
        dnkupa = new double[num_points];

        dnkuplocal = new double*[num_points];
        for(int i = 0; i < num_points; i++) dnkuplocal[i] = new double[dmnpop];
    }

    void SetZero(){
        setzero(dnkup);
    }
    void SetZeroAx(){
        setzero(dnkupa);
    }

    void NotAccept(int ipop){
        for(int i = 0; i < stored_num_points; i++)
            dnkupa[i] = dnkuplocal[i][ipop];
    }

    void WalkerMatch(int jpop){
        for(int i = 0; i < stored_num_points; i++)
        {
            dnkuplocal[i][jpop] = dnkupa[i];
        }
    }

    void WalkerCollect(int nsons){
        for(int i = 0; i < stored_num_points; i++)
            dnkup[i] = dnkup[i] + dnkupa[i]*nsons;
    }

    void Normalization(int np, int nkuppt){
        for(int ik = 0; ik < stored_num_points;ik++)
            dnkup[ik] = dnkup[ik]*np/(float(nkuppt)); //why np is not clear
    }


    ~MomentDistr(){
        delete [] dnkup;
        delete [] dnkupa;

        for(int i = 0; i < stored_num_points; i++)
            delete [] dnkuplocal[i];
        delete [] dnkuplocal;
    }
};

class CorFun{
public:
    double step, max_value, num_points;
};

struct ParamModel{
    int num_comp, num_part, num_walk;
    long seed;
    double scat_lenght, scat_lenght_bos, width, alfa;
};


class Coordinates{
private:
    int num_comp_saved, num_part_saved, num_walk_saved;
    long seed_saved;
    double alfa_saved, step_jump;
public:
    double ****x, **xaux, **xMT;
    Coordinates(ParamModel param_model){

        num_comp_saved = param_model.num_comp;
        num_part_saved = param_model.num_part;
        num_walk_saved = param_model.num_walk;
        seed_saved = param_model.seed;
        alfa_saved = param_model.alfa;
        step_jump = alfa_saved/4.0;

        xaux = new double*[num_comp_saved];
        for(int i = 0; i < num_comp_saved; i++) xaux[i] = new double[num_part_saved];

        x = new double***[num_comp_saved];
        for(int i = 0; i < num_comp_saved; i++)
        {
            x[i] = new double**[num_part_saved];
            for(int j = 0; j < num_part_saved; j++)
            {
                 x[i][j] = new double*[dmnpop];
                 for(int m = 0; m < dmnpop; m++)
                    x[i][j][m] = new double[2];
            }
        }

        xMT = new double*[num_comp_saved];
        for(int i = 0; i < num_comp_saved; i++) xMT[i] = new double[num_part_saved];

//        xauxT = new double*[num_comp_saved];
//        for(int i = 0; i < num_comp_saved; i++) xauxT[i] = new double[num_part_saved];

   }
   void ReadInitial(const string & startingConfig){
       double etemp;
       fstream filestr(startingConfig, fstream::in | fstream::out);
       filestr>>setprecision(18);
       filestr >> num_walk_saved;
       //cout << num_walk_saved;
       for (int ipop = 0; ipop < num_walk_saved; ipop++ )
       {
           filestr >> etemp;
           for(int ic = 0; ic < num_comp_saved;ic++)
           {
               for(int ip = 0; ip < num_part_saved; ip++)
               {
                   filestr >> x[ic][ip][ipop][0];
                   //if(ipop == 0)
                     //  cout << x[ic][ip][ipop][0] << "\n";
               }
           }

       }
       filestr.close();
   }

   void GenerateInitial()
   {
       long seed = -13156;
       for (int ipop = 0; ipop < num_walk_saved; ipop++ )
       {
           for(int ic = 0; ic < num_comp_saved; ic++)
           {
               for(int ip = 0; ip < num_part_saved; ip++) {x[ic][ip][ipop][0] = ran2(&seed);}

           }
       }

       fstream outfile("../1D_Qt/Data/inprev.dat", fstream::out );
       outfile << setprecision(17);
       outfile<<num_walk_saved<<"\n";
       for(int ipop = 0; ipop < num_walk_saved;ipop++)
       {
           outfile<<ipop<<"\n";
           for(int ic = 0; ic < num_comp_saved; ic++)
           {
               for(int ip = 0; ip < num_part_saved; ip++)
                   outfile<<x[ic][ip][ipop][0]<<"\n";
           }
       }
       outfile.close();
   }

   void GaussianJump(int ntemps, int in, int i_VMC, int ipop, double ****FF){
       double xgaus;

       for (int ic = 0; ic < num_comp_saved; ic++)
       {
           for(int ip = 0; ip < num_part_saved; ip++)
           {
               if (ntemps > 1)
                   Gauss1D(&xgaus, alfa_saved, &seed_saved);
               else
                   xgaus = 0.0;

               if (i_VMC == 1){
                   xMT[ic][ip] = x[ic][ip][ipop][in] + xgaus;

               }
               else
                   xMT[ic][ip] = x[ic][ip][ipop][in] + xgaus + FF[ic][ip][ipop][in] * step_jump;
           }
        }
   }

   void Accept(){
       for(int ic = 0; ic < num_comp_saved; ic++)
       {
           for(int ip = 0; ip < num_part_saved; ip++)
           {
               xaux[ic][ip] = xMT[ic][ip];
           }
       }
   }

    void NotAccept(int ipop, int in){
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++)
                xaux[ic][ip] = x[ic][ip][ipop][in];
        }
    }

    void WalkerMatch(int jpop, int io){
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++)
                x[ic][ip][jpop][io] = xaux[ic][ip];
         }
    }


   ~Coordinates(){
       for(int i = 0; i < num_comp_saved; i++)
       {
           for(int j = 0; j < num_part_saved; j++)
           {
               for(int m = 0; m < dmnpop; m++)
                   delete [] x[i][j][m];

               delete[] x[i][j];
            }
              delete [] x[i];
        }
        delete [] x;

       for(int i = 0; i < num_comp_saved; i++)
         delete [] xaux[i];
       delete [] xaux;


       for(int i = 0; i < num_comp_saved; i++)
         delete [] xMT[i];
       delete [] xMT;

   }

};




void run (const string & inFile, const string & startingConfig, const string & outDir )
{
    Energy eav, epar, emean, enew, eold, emtnew;
    double ewalk, etemp; // previously it was E
	long  kkk = -2175;
    double alfa;
    double PsiTotal, fnew, fvella;
    int ntemps, accepta, nacc, nprova, nwrite, nmean;
    double xgaus, accrate;

    CorFun pair_distr, obdm, dens_distr, mom_distr;

    double Lmax;

    int in, io, ii, jpop, nsons, npopmean;

    ParamModel param_model;
    param_model.seed = kkk;



//	double a; // multi-component scattering length in oscillatory units
//	double aB; //in-component scattering length in oscillatory units, for fermions aB = 0

//	int np;
    int nblck, niter;

//	double dt, dte;

    //int npop;

	int init;

//	int ncomp;

	int icrit;

//	double width;

	string str;

	int i_VMC, i_Smart = 0, i_Drift = 0, i_FNDMC, i_OBDM_1, i_OBDM_2=0, i_stat_cor;

    fstream filestr2(inFile, fstream::in | fstream::out);
    filestr2>>str>>i_VMC>>str>>i_FNDMC>>str>>param_model.num_comp;
    filestr2>>str>>param_model.num_part>>str>>param_model.width;
    filestr2>>str>>param_model.scat_lenght_bos>>str>>param_model.scat_lenght;
    filestr2>>str>>param_model.alfa>>str>>niter>>str>>nblck>>str;
    filestr2>>param_model.num_walk>>str>>init>>str>>nwrite>>str>>icrit;
    filestr2>>str>>i_OBDM_1>>str>>pair_distr.num_points>>str>>obdm.num_points;
    filestr2>>str>>dens_distr.num_points>>str>>pair_distr.max_value;
    filestr2>>str>>obdm.max_value>>str>>dens_distr.max_value;
    filestr2>>str>>mom_distr.num_points>>str>>mom_distr.max_value>>str>>Lmax;
    filestr2.close();


    Coordinates coordinates(param_model);

    pair_distr.step = pair_distr.max_value/pair_distr.num_points;
    obdm.step = obdm.max_value/obdm.num_points;
    dens_distr.step = dens_distr.max_value/dens_distr.num_points;

    mom_distr.step = mom_distr.max_value/mom_distr.num_points;

    fstream outfile2(outDir+"/Output.dat", fstream::out| ios::app );
	outfile2<<setprecision(18);	
    outfile2<<" i_VMC= "<<i_VMC<<" i_FNDMC= "<<i_FNDMC<<" ncomp= "<<param_model.num_comp;
    outfile2<<" np= "<<param_model.num_part<<" width= "<<param_model.width;
    outfile2<<" aB= "<<param_model.scat_lenght_bos<<" a= "<<param_model.scat_lenght<<"\n";
    outfile2<<" alfa= "<<param_model.alfa<<" niter= "<<niter<<" nblck= ";
    outfile2<<nblck<<" npop= "<<param_model.num_walk;
    outfile2<<" init= "<<init<<" nwrite= "<<nwrite<<" icrit= "<<icrit<<"\n";
	outfile2<<" i_OBDM= "<<i_OBDM_1<<"\n";
    outfile2<<" mgr_g(r)= "<<pair_distr.num_points<<" mgr_OBDM= "<<obdm.num_points;
    outfile2<<" mgr_dens= "<<dens_distr.num_points<<"\n";
    outfile2<<" Lmax_g(r)= "<<pair_distr.max_value<<" Lmax_OBDM= "<<obdm.max_value;
    outfile2<<" Lmax_dens= "<<dens_distr.max_value<<"\n";
    outfile2<<" numks= "<<mom_distr.num_points<<" kmax= "<<mom_distr.max_value;
    outfile2<<" Lmax_McM= "<<Lmax<<"\n";

	outfile2.close();

//	alfa = dt;
//	dte = dt * 0.25;



//************************************************
//	double ****x;
//	double ** xaux;
	double ** F;
	
//	double **elocal, **elocal1, **elocal2;
    Energy **elocal;


//************************************************
//	double **xauxT;
	double ****FF;
	double **flocal;
//	double **xMT;
    double **FMT;

    DistributionR pair_distr_1(pair_distr.num_points);
    DistributionR pair_distr_2(pair_distr.num_points);
    DistributionR pair_distr_12(pair_distr.num_points);

    DistributionR dens_distr_1(dens_distr.num_points);
    DistributionR dens_distr_2(dens_distr.num_points);

    OBDM obdm_1(obdm.num_points), obdm_2(obdm.num_points);

    double nkuppt;

    MomentDistr moment_distr_1(mom_distr.num_points), moment_distr_2(mom_distr.num_points);

	double *faiJ;

	double *kr; 
	
    int numks_temp = mom_distr.num_points;
    kr = new double[numks_temp];

    faiJ = new double[param_model.num_part];

    int ngr;

/*    xaux = new double*[ncomp];
	for(int i = 0; i < ncomp; i++) xaux[i] = new double[np];

	x = new double***[ncomp];
	for(int i = 0; i < ncomp; i++)
	{
		x[i] = new double**[np];
		for(int j = 0; j < np; j++)
		{
			 x[i][j] = new double*[dmnpop];
			 for(int m = 0; m < dmnpop; m++)
				x[i][j][m] = new double[2];
		}
    }*/

    FF = new double***[param_model.num_comp];
    for(int i = 0; i < param_model.num_comp; i++)
	{
        FF[i] = new double**[param_model.num_part];
        for(int j = 0; j < param_model.num_part; j++)
		{
			 FF[i][j] = new double*[dmnpop];
			 for(int m = 0; m < dmnpop; m++)
				FF[i][j][m] = new double[2];
		}
	}

	flocal = new double*[dmnpop];
	for(int i = 0; i < dmnpop; i++) flocal[i] = new double[2];

/*	xMT = new double*[ncomp];
    for(int i = 0; i < ncomp; i++) xMT[i] = new double[np];*/

    FMT = new double*[param_model.num_comp];
    for(int i = 0; i < param_model.num_comp; i++) FMT[i] = new double[param_model.num_part];

    F = new double*[param_model.num_comp];
    for(int i = 0; i < param_model.num_comp; i++) F[i] = new double[param_model.num_part];

/*	xauxT = new double*[ncomp];
    for(int i = 0; i < ncomp; i++) xauxT[i] = new double[np];*/

    elocal = new Energy*[dmnpop];
    for(int i = 0; i < dmnpop; i++) elocal[i] = new Energy[2];


//    fold = new double[dmnpop];

	

//********************************************************
    if(init == 1)
    {
/*        fstream filestr(startingConfig, fstream::in | fstream::out);
		filestr>>setprecision(18);
		filestr >> npop;
		for (int ipop = 0; ipop < npop; ipop++ )
		{
			filestr >> etemp;
			for(int ic = 0; ic < ncomp;ic++)
			{
   				for(int ip = 0; ip < np; ip++)
				{
			  		filestr >> x[ic][ip][ipop][0];
					//if(ipop==96){cout<<x[ic][ip][ipop][0]<<"\n";}
				}
			}
		}
        filestr.close(); */
    coordinates.ReadInitial(startingConfig);
    }
    //else {Initial(x, npop, ncomp, np);}
    else {coordinates.GenerateInitial();}


//*********************************************************


	ntemps = 0;
	
	in = 0;
	io = 1;	

	nacc = 0;
	nprova = 0;

    for(int j = 0; j < param_model.num_walk; j++ ) // It will be in constructor of Force
	{
        for(int ic = 0; ic < param_model.num_comp; ic++)
		{
            for(int ip = 0; ip < param_model.num_part; ip++)
			{
				FF[ic][ip][j][in] = 0.0;
			}
		}
	}

//    cout << coordinates.x[0][0][0][0] <<"\n";

for(int iblck = 0; iblck < nblck; iblck++)
{


	ngr = 0;

    pair_distr_1.SetZero();
    pair_distr_2.SetZero();
    pair_distr_12.SetZero();

    obdm_1.SetZero();
    obdm_2.SetZero();

    dens_distr_1.SetZero();
    dens_distr_2.SetZero();

    moment_distr_1.SetZero();
    moment_distr_2.SetZero();

	nkuppt = 0;

    for(int iter = 0; iter < niter; iter++)
   	{
 //           cout << iter <<"\n";
     		ntemps = ntemps + 1;

            eav.SetZero();

            jpop = 0;
 
            for(int ipop = 0; ipop < param_model.num_walk; ipop++)
		{

            for(int ic = 0; ic < param_model.num_comp; ic++) // Force::SetZero
			{
                for(int ip = 0; ip < param_model.num_part; ip++)
				{
                    FMT[ic][ip] = 0.0;
				}
			}

            emtnew.SetZero();

            eold.tot = elocal[ipop][in].tot; // Energy::SetOldConf
            eold.kin = elocal[ipop][in].kin;
            eold.pot = elocal[ipop][in].pot;

            pair_distr_1.SetZeroAx();
            pair_distr_2.SetZeroAx();
            pair_distr_12.SetZeroAx();

            dens_distr_1.SetZeroAx();
            dens_distr_2.SetZeroAx();

            obdm_1.SetZeroAx();
            obdm_2.SetZeroAx();

            for (int i = 0; i < param_model.num_part; i++) faiJ[i] = 0.0; // To make a method

            moment_distr_1.SetZeroAx();
            moment_distr_2.SetZeroAx();



            coordinates.GaussianJump(ntemps, in, i_VMC, ipop, FF);

//            cout << coordinates.xMT[0][0] <<"\n";

//              cout << coordinates.x[0][0][0][0] <<"\n";
            // Coordinate::GaussJump

/*			for (int ic = 0; ic < ncomp; ic++)
			{
				for(int ip = 0; ip < np; ip++)
				{
	  				if (ntemps > 1)	Gauss1D(&xgaus, alfa, &kkk);
	  				else {xgaus = 0.0;}	
					if (i_VMC == 1) {xMT[ic][ip] = x[ic][ip][ipop][in] + xgaus;}
					else		
					 {xMT[ic][ip] = x[ic][ip][ipop][in] + xgaus + FF[ic][ip][ipop][in] * dte;
					//cout<<setprecision(18);	
					 //cout<<ic<<" "<<ip<<" "<<xgaus<<" "<<xMT[ic][ip]<<"\n";
					}
				}
            } */

            // ...

            //          WaveFunction(ncomp, np, aB, a, width, xMT, &PsiTotal); // To leave as function
            WaveFunction(param_model.num_comp, param_model.num_part, param_model.scat_lenght_bos, param_model.scat_lenght, param_model.width, coordinates.xMT, &PsiTotal); // To leave as function

//            WaveFunction(param_model, coordinates, &PsiTotal); // To leave as function



//            Energy_calc(xMT, FMT, emtnew, ncomp, np, aB, a, width); // To leave as function
            Energy_calc(coordinates.xMT, FMT, emtnew, param_model.num_comp, param_model.num_part, param_model.scat_lenght_bos, param_model.scat_lenght, param_model.width); // To leave as function

//            Energy_calc(coordinates, FMT, emtnew, param_model); // To leave as function


            // To make as methods
            PairDistribution_calc(coordinates.xMT, pair_distr_1.draMT, pair_distr_2.draMT, pair_distr_12.draMT, pair_distr.num_points, pair_distr.step, param_model.num_comp, param_model.num_part);
            DensityDistribution_calc( coordinates.xMT, dens_distr_1.draMT, dens_distr_2.draMT, dens_distr.num_points, dens_distr.step, param_model.num_comp, param_model.num_part);

			if(i_Drift == 0)
			{
            //MetropolisDif(ipop, ncomp, np, PsiTotal, flocal, xMT, x, FF, FMT, &accepta, &nprova, &fvella, ntemps, &kkk, in, dte, i_VMC);
                MetropolisDif(ipop, param_model.num_comp, param_model.num_part, PsiTotal, flocal, coordinates.xMT, coordinates.x, FF, FMT, &accepta, &nprova, &fvella, ntemps, &param_model.seed, in, param_model.alfa/4.0, i_VMC);

            }
			if(i_Drift == 1)
			{
			accepta = 1;
			}

//            cout << "accepta= " << accepta <<"\n";
            if(accepta == 1)
			{
				if(ntemps > 1) nacc = nacc + 1;

                fnew = PsiTotal; // Accept method of WaveFunction

                enew.tot = emtnew.tot; //Accept method of Energy
                enew.kin = emtnew.kin;
                enew.pot = emtnew.pot;

                coordinates.Accept();


                for(int ic = 0; ic < param_model.num_comp; ic++)
				{
                    for(int ip = 0; ip < param_model.num_part; ip++)
					{
                        //xaux[ic][ip] = xMT[ic][ip]; // Make Accept method
                        F[ic][ip] = FMT[ic][ip]; // Make Accept method, F should be replaced to Faux
					}
				}

                pair_distr_1.Accept();
                pair_distr_2.Accept();
                pair_distr_12.Accept();

                dens_distr_1.Accept();
                dens_distr_2.Accept();

				if(i_OBDM_1 == 1)
				{
					
                //OBDM1D_11(Lmax,ncomp,np,width, aB, a, kr, xaux,fra,nfra,PsiTotal, dnkupa, numks, dk, &kkk, mgr2, dr2);
                    OBDM1D_11(Lmax,param_model.num_comp,param_model.num_part,param_model.width, param_model.scat_lenght_bos, param_model.scat_lenght, kr, coordinates.xaux, obdm_1.fra, obdm_1.nfra, PsiTotal, moment_distr_1.dnkupa, mom_distr.num_points, mom_distr.step, &param_model.seed, obdm.num_points, obdm.step);

                }
			//	if(i_OBDM_2 == 1)
			//	{
            //	OBDM1D_22(Lmax,ncomp,np,width, aB, a, kr, xaux,fra_22,nfra_22,PsiTotal, dnkupa_22, numks, dk, &kkk, mgr2, dr2);
			//	}

			} 
			else 
			{
                fnew = fvella; // NotAccept method of WaveFunction

                enew.tot = eold.tot; // NotAccept method of energy
                enew.tot = eold.tot;
                enew.tot = eold.tot;

                coordinates.NotAccept(ipop, in);

                for(int ic = 0; ic < param_model.num_comp; ic++)
				{
                    for(int ip = 0; ip < param_model.num_part; ip++)
					{
                        //xaux[ic][ip] = x[ic][ip][ipop][in];
						F[ic][ip] = FF[ic][ip][ipop][in];
					}
				}


/*                for(int i = 1; i < (mgr1+1); i++)
				{
					gra_11[i] = grlocal_11[i][ipop][in];
					gra_22[i] = grlocal_22[i][ipop][in];
					gra_12[i] = grlocal_12[i][ipop][in];
                }*/

                // 3 loops instead of one?

                pair_distr_1.NotAccept(ipop, in);
                pair_distr_2.NotAccept(ipop, in);
                pair_distr_12.NotAccept(ipop, in);


                dens_distr_1.NotAccept(ipop, in);
                dens_distr_2.NotAccept(ipop, in);

                obdm_1.NotAccept(ipop);
                obdm_2.NotAccept(ipop);

                moment_distr_1.NotAccept(ipop);
                moment_distr_2.NotAccept(ipop);

			}

			if (i_FNDMC == 1)
			{
            BranchingCalc(&nsons, accepta, ntemps, nacc, nprova, param_model.alfa/4.0, icrit, ewalk, enew.tot, eold.tot, &kkk, param_model.num_walk);
			}

			else {nsons = 1;}
//			cout<<nsons;

//			nsons = 1;

			if(nsons > 0)
			{
				for(int js = 0; js < nsons; js++)
				{
					
					if(jpop >= dmnpop)
					{
		cout <<"The population has grown too much!"<<"nsons="<<nsons<<"jpop="<<jpop<<"\n";
					}

                    // Begining of WalkerMatch for Energy, Coordinate and WaveFunction

                    elocal[jpop][io].tot = enew.tot;
                    elocal[jpop][io].kin = enew.kin;
                    elocal[jpop][io].pot = enew.pot;


                    flocal[jpop][io] = fnew;

                    coordinates.WalkerMatch(jpop, io);

                    for(int ic = 0; ic < param_model.num_comp; ic++)
					{
                        for(int ip = 0; ip < param_model.num_part; ip++)
						{
                        //	x[ic][ip][jpop][io] = xaux[ic][ip];
							FF[ic][ip][jpop][io] = F[ic][ip];
						}
					}

                    // End of WalkerMatch for Energy, Coordinate and WaveFunction

                    pair_distr_1.WalkerMatch(jpop, io);
                    pair_distr_2.WalkerMatch(jpop, io);
                    pair_distr_12.WalkerMatch(jpop, io);

                    dens_distr_1.WalkerMatch(jpop, io);
                    dens_distr_2.WalkerMatch(jpop, io);

                    obdm_1.WalkerMatch(jpop);
                    obdm_2.WalkerMatch(jpop);

                    moment_distr_1.WalkerMatch(jpop);
                    moment_distr_2.WalkerMatch(jpop);

					jpop = jpop + 1;				
				}
			}

            eav.tot = eav.tot + enew.tot * nsons; // WalkerCollect for Energy

            ngr = ngr + nsons;

            pair_distr_1.WalkerCollect(nsons);
            pair_distr_2.WalkerCollect(nsons);
            pair_distr_12.WalkerCollect(nsons);

            dens_distr_1.WalkerCollect(nsons);
            pair_distr_2.WalkerCollect(nsons);

            obdm_1.WalkerCollect(nsons);
            obdm_2.WalkerCollect(nsons);

            moment_distr_1.WalkerCollect(nsons);
            moment_distr_2.WalkerCollect(nsons);

            nkuppt = nkuppt + nsons * 100; //Number of McMillan points


		}
	
        param_model.num_walk = jpop;
		
		ii = in;
		in = io;
		io = ii;

        // Energy::Normalization

        ewalk = eav.tot/jpop;
	   
        epar.tot = eav.tot /(jpop*param_model.num_comp*param_model.num_part);
        epar.kin = eav.kin/(jpop*param_model.num_comp*param_model.num_part);
        epar.pot = eav.pot/(jpop*param_model.num_comp*param_model.num_part);

        // ...

		if(ntemps == 1)
		{
            emean.SetZero();

			npopmean = 0;
			nmean = 0;
		}

		nmean = nmean + 1;

        // Energy::Average

        emean.tot = emean.tot + epar.tot;
        emean.kin = emean.kin + epar.kin;
        emean.pot = emean.pot + epar.pot;

        // ...

        npopmean = npopmean + param_model.num_walk;

		if(nmean==nwrite)
		{
            //Energy::Print

            fstream outfile(outDir + "/Energy.dat", fstream::out| ios::app );

			outfile<<setprecision(18);

            outfile<< ntemps/nwrite <<" "<<emean.tot/float(nwrite)<<" "<<emean.kin/float(nwrite)<<" "<<emean.pot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";
            cout<< ntemps/nwrite <<" "<<emean.tot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";

            outfile.close();

            emean.SetZero();

            nmean = 0;
			npopmean = 0;	
	 	}
	}
//		 cout << ntemps <<" "<<Epar/(2.0*dd)<<" "<<Epar1/(2.0*dd)<<"\n";		

//********************************************************************************************

    pair_distr_1.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr.step);
    pair_distr_2.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr.step);
    pair_distr_12.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr.step);
	
    double r1;


    // It should be CorFun::Print

    std::fstream outfile9(outDir + "/gr.dat", std::fstream::out| std::ios::app );
    for(int i = 1; i < (pair_distr.num_points+1);i++)
	{
        r1 = float(i) * pair_distr.step + pair_distr.step/2.0; //(float(i)+1) is in order to coincide with fortran code
            outfile9<<r1<<" "<<pair_distr_1.dr[i]<<"\n";
	       // cout<<r1<<" "<<gr_11[i]<<" "<<gr_12[i]<<"\n"; 	
	}
	outfile9.close();

    std::fstream outfile19(outDir + "/gr_12.dat", std::fstream::out| std::ios::app );
    for(int i = 1; i < (pair_distr.num_points+1);i++)
	{
        r1 = float(i) * pair_distr.step + pair_distr.step/2.0; //(float(i)+1) is in order to coincide with fortran code
            outfile19<<r1<<" "<<pair_distr_12.dr[i]<<"\n";
	       // cout<<r1<<" "<<gr_11[i]<<" "<<gr_12[i]<<"\n"; 	
	}
	outfile19.close();

    dens_distr_1.NormalizationNR(ngr, dens_distr.step);
    dens_distr_2.NormalizationNR(ngr, dens_distr.step);




    std::fstream outfile22(outDir + "/nr.dat", std::fstream::out| std::ios::app );
    for(int i = 1; i < (dens_distr.num_points+1);i++)
	{
        r1 = float(i) * dens_distr.step + dens_distr.step/2.0; //(float(i)+1) is in order to coincide with fortran code
//	        outfile22<<r1<<" "<<nr_11[i]<<" "<<nr_22[i]<<"\n"; 	
            outfile22<<r1<<" "<<dens_distr_1.dr[i]<<"\n";
	}
	outfile22.close();

    obdm_1.Normalization();
    obdm_2.Normalization();

    fstream outfile13(outDir + "/fr.dat", fstream::out| ios::app );
    for(int i = 1; i < (obdm.num_points+1);i++)
	{
        r1 = float(i) * obdm.step + obdm.step/2.0;
//		outfile13<<r1<<" "<<fr[i]<<" "<<fr_22[i]<<"\n";		
        outfile13<<r1<<" "<<obdm_1.fr[i]<<"\n";
	}
	outfile13.close();

	
//	double n0 =  dnkup[0]*np/float(nkuppt);
//	double n0_22 =  dnkup_22[0]*np/float(nkuppt);

	double k1;

    moment_distr_1.Normalization(param_model.num_part, nkuppt);
    moment_distr_2.Normalization(param_model.num_part, nkuppt);

    fstream outfile14(outDir + "/nk.dat", fstream::out| ios::app );
    for(int ik = 0; ik < mom_distr.num_points;ik++)
	{

        k1 = float(ik) * mom_distr.step;
        outfile14<<k1<<" "<<moment_distr_1.dnkup[ik]<<"\n";
		//cout<<k1<<" "<<dnkup[ik]<<"\n";
	
	}
	outfile14.close();



//***************************************************************************

    fstream outfile2(startingConfig, fstream::out ); //Coordinate::Print
    outfile2<<param_model.num_walk<<"\n";
	outfile2<<setprecision(18);
    for(int i = 0; i < param_model.num_walk;i++)
	{
        outfile2<<elocal[i][in].tot<<"\n";
        for(int ic = 0; ic < param_model.num_comp; ic++)
			{
                for(int ip = 0; ip < param_model.num_part; ip++ ) {outfile2<<coordinates.x[ic][ip][i][in]<<"\n"; }
			}
	}
	outfile2.close();

	accrate = 100.0 * float(nacc)/float(nprova);
    fstream outfile1(outDir + "/Accept.dat", fstream::out| ios::app );
	outfile1<< iblck <<" "<<accrate<<"\n";
	outfile1.close();

//*************************************************************** //
//  Here we take the average of all observables after each block *
//**************************************************************//

	int n_notused,  nfull, nfull_rew, ncol_rew, ncol, nel;

    // Parameters for the average of the energy//

    ncol = 4;
    n_notused = 0;
    //cout<<"n_notused= "<<n_notused<<"\n";
    nel = 1; // number of outputs for one block of statistics
    nfull = 0;

    Stat_Energy(outDir + "/Energy.dat", outDir + "/Energy_av.dat", outDir + "/Energy_full.dat", ncol, nel, n_notused,  nfull);

//  Parameters for correlation functions
//  (files "nr.dat", "mom_distr.dat", "fr.dat", "g2D.dat" )
    int nel_cor, ncol_cor,  n_notused_cor;
    
    nel_cor = 1; // number of elements per block for statistics
    n_notused_cor = 1; 
    
    ncol_cor = 1; // number of columns for 'mom_dist.dat'
    statistics_cor(outDir + "/nk.dat", outDir + "/nk_av.dat", mom_distr.num_points, ncol_cor, n_notused_cor, nel_cor);

    ncol_cor = 1; //number of columns for 'gr.dat'
    statistics_cor(outDir + "/gr.dat", outDir + "/gr_av.dat", pair_distr.num_points, ncol_cor, n_notused_cor, nel_cor);

    ncol_cor = 1; //number of columns for 'gr.dat'
    statistics_cor(outDir + "/gr_12.dat", outDir + "/gr_12_av.dat", pair_distr.num_points, ncol_cor, n_notused_cor, nel_cor);

    ncol_cor = 1; //number of columns for 'fr.dat'
    statistics_cor(outDir + "/fr.dat", outDir + "/fr_av.dat", obdm.num_points, ncol_cor, n_notused_cor, nel_cor);

    ncol_cor = 1; //number of columns for 'nr.dat'
    statistics_cor(outDir + "/nr.dat", outDir + "/nr_av.dat", dens_distr.num_points, ncol_cor, n_notused_cor, nel_cor);


    Normalization(outDir + "/nk_av.dat", outDir + "/nk_av_norm.dat", mom_distr.num_points, param_model.num_part, param_model.num_comp);
    Normalization(outDir + "/nr_av.dat", outDir + "/nr_av_norm.dat", dens_distr.num_points, param_model.num_part, param_model.num_comp);

/*    if(i_FNDMC == 1)
    {
	Extrapolation("../VMC/Measur/nk_av_norm.dat","Measur/nk_av_norm.dat", "../Extr/nk_extr.dat", "../Extr/nk_extr_lin.dat", numks);
    Extrapolation("../VMC/Measur/gr_av.dat","Measur/gr_av.dat", "../Extr/gr_extr.dat", "../Extr/gr_extr_lin.dat", mgr1);
    Extrapolation("../VMC/Measur/gr_12_av.dat","Measur/gr_12_av.dat", "../Extr/gr_12_extr.dat", "../Extr/gr_12_extr_lin.dat", mgr1);
    Extrapolation("../VMC/Measur/fr_av.dat","Measur/fr_av.dat", "../Extr/fr_extr.dat", "../Extr/fr_extr_lin.dat", mgr2);
    Extrapolation("../VMC/Measur/nr_av_norm.dat","Measur/nr_av_norm.dat", "../Extr/nr_extr.dat", "../Extr/nr_extr_lin.dat", mgr3);
    }*/
}

/*	for(int i = 0; i < ncomp; i++)
	{
	  for(int j = 0; j < np; j++)
		{
		    	for(int m = 0; m < dmnpop; m++)
			delete [] x[i][j][m];

			delete[] x[i][j];					
		}
	   delete [] x[i];
	}
    delete [] x;*/

    for(int i = 0; i < param_model.num_comp; i++)
	{
      for(int j = 0; j < param_model.num_part; j++)
		{
		    	for(int m = 0; m < dmnpop; m++)
			delete [] FF[i][j][m];

			delete[] FF[i][j];					
		}
	   delete [] FF[i];
	}
	delete [] FF;

    /*for(int i = 0; i < ncomp; i++)
	  delete [] xaux[i];
	delete [] xaux;	
	
	
	for(int i = 0; i < ncomp; i++)
	  delete [] xMT[i];
    delete [] xMT; */

    for(int i = 0; i < param_model.num_comp; i++)
	  delete [] F[i];
	delete [] F;

	for(int i = 0; i < dmnpop; i++)
	  delete [] flocal[i];
	delete [] flocal;		

	for(int i = 0; i < dmnpop; i++)
	  delete [] elocal[i];
	delete [] elocal;	

//	delete [] fold;

	delete [] faiJ;	

	delete [] kr;

}
