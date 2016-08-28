// Both for bosons and fermions,
// h^2/(m a_osc^2) = 1, D = 1/2


#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "Locals.h"
#include "Algorithm.h"
#include "Energy.h"
#include "g2D.h"
#include "dens.h"
#include "OBDM.h"
#include "MomDistr.h"
#include "DistrR.h"
#include "Statistics.h"
#include "Wave_fun.h"
#include "qmc.h"

using namespace std;

void Energy::SetZero(){
    tot = 0; pot = 0; kin = 0;
};

void run (const string & inFile, const string & startingConfig, const string & outDir )
{
    Energy eav, epar, emean, enew, eold, emtnew;
    double ewalk;
	long  kkk = -2175;
    double PsiTotal, fnew, fvella;
    int ntemps, accepta, nacc, nprova, nwrite, nmean;
    double accrate;

    CorFunParam pair_distr, obdm, dens_distr, mom_distr;

    double Lmax;

    int in, io, ii, jpop, nsons, npopmean;

    ParamModel param_model;
    param_model.seed = kkk;

    int nblck, niter;

	int init;
	int icrit;
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


    Locals coordinates(param_model);
    Locals force(param_model);

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
    Energy **elocal;
    double **flocal;

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

	flocal = new double*[dmnpop];
	for(int i = 0; i < dmnpop; i++) flocal[i] = new double[2];

    elocal = new Energy*[dmnpop];
    for(int i = 0; i < dmnpop; i++) elocal[i] = new Energy[2];

    if(init == 1)
        coordinates.ReadInitial(startingConfig);
    else
        coordinates.GenerateInitial();


	ntemps = 0;
	
	in = 0;
	io = 1;	

	nacc = 0;
	nprova = 0;

    force.SetZeroForceTotal(in);

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

            force.SetZeroForce();

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

            for (int i = 0; i < param_model.num_part; i++)
                faiJ[i] = 0.0; // To make a method

            moment_distr_1.SetZeroAx();
            moment_distr_2.SetZeroAx();

            coordinates.GaussianJump(ntemps, in, i_VMC, ipop, force.total);

            PsiTotal = WaveFunction(param_model, coordinates); // To leave as function

//            Energy_calc(xMT, FMT, emtnew, ncomp, np, aB, a, width); // To leave as function
            Energy_calc(coordinates.metrop, force.metrop, emtnew, param_model.num_comp, param_model.num_part, param_model.scat_lenght_bos, param_model.scat_lenght, param_model.width); // To leave as function

//            Energy_calc(coordinates, FMT, emtnew, param_model); // To leave as function


            // To make as methods
            PairDistribution_calc(coordinates.metrop, pair_distr_1.draMT, pair_distr_2.draMT, pair_distr_12.draMT, pair_distr.num_points, pair_distr.step, param_model.num_comp, param_model.num_part);
            DensityDistribution_calc( coordinates.metrop, dens_distr_1.draMT, dens_distr_2.draMT, dens_distr.num_points, dens_distr.step, param_model.num_comp, param_model.num_part);

			if(i_Drift == 0)
			{
            //MetropolisDif(ipop, ncomp, np, PsiTotal, flocal, xMT, x, FF, FMT, &accepta, &nprova, &fvella, ntemps, &kkk, in, dte, i_VMC);
                MetropolisDif(ipop, param_model.num_comp, param_model.num_part, PsiTotal, flocal, coordinates.metrop, coordinates.total, force.total, force.metrop, &accepta, &nprova, &fvella, ntemps, &param_model.seed, in, param_model.alfa/4.0, i_VMC);

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

                force.Accept();

                pair_distr_1.Accept();
                pair_distr_2.Accept();
                pair_distr_12.Accept();

                dens_distr_1.Accept();
                dens_distr_2.Accept();

				if(i_OBDM_1 == 1)
				{
					
                //OBDM1D_11(Lmax,ncomp,np,width, aB, a, kr, xaux,fra,nfra,PsiTotal, dnkupa, numks, dk, &kkk, mgr2, dr2);
                    OBDM1D_11(Lmax,param_model.num_comp,param_model.num_part,param_model.width, param_model.scat_lenght_bos, param_model.scat_lenght, kr, coordinates.auxil, obdm_1.fra, obdm_1.nfra, PsiTotal, moment_distr_1.dnkupa, mom_distr.num_points, mom_distr.step, &param_model.seed, obdm.num_points, obdm.step);

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
                force.NotAccept(ipop, in);


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
                    force.WalkerMatch(jpop, io);

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
                for(int ip = 0; ip < param_model.num_part; ip++ ) {outfile2<<coordinates.total[ic][ip][i][in]<<"\n"; }
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

	for(int i = 0; i < dmnpop; i++)
	  delete [] flocal[i];
	delete [] flocal;		

	for(int i = 0; i < dmnpop; i++)
	  delete [] elocal[i];
	delete [] elocal;	

    delete [] faiJ;

	delete [] kr;
}
