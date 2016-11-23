// Both for bosons and fermions,
// h^2/(m a_osc^2) = 1, D = 1/2


#include <iostream>
#include <map>
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
#include "utils.h"

using namespace std;

void Energy::SetZero(){
    tot = 0; pot = 0; kin = 0;
};

void run (const string & inFile, const string & startingConfig, const string & outDir )
{
    Energy eav, epar, emean, enew, eold, emtnew;
    double ewalk;
    int ntemps, accepta, nacc, nprova, nwrite, nmean;
    double accrate;

    int in, io, ii, jpop, nsons, npopmean;
    int nblck, niter, init, nwalk_mean;
    int i_VMC, i_Drift, i_FNDMC, i_OBDM;

    map<string, double> paramMap;

    paramMap = FillMap(inFile, outDir);

    i_VMC   = paramMap["i_VMC"];
    i_Drift = paramMap["i_Drift"];
    i_FNDMC = paramMap["i_FNDMC"];
    niter = paramMap["niter"];
    nblck = paramMap["nblck"];
    init = paramMap["init"];
    nwrite = paramMap["nwrite"];
    nwalk_mean = paramMap["icrit"];
    i_OBDM = paramMap["i_OBDM"];

    ParamModel param_model(paramMap);

    CorFunParam pair_distr_param(paramMap["Lmax_gr"], paramMap["mgr_g(r)"]);
    CorFunParam obdm_param(paramMap["Lmax_OBDM"], paramMap["mgr_OBDM"]);
    CorFunParam dens_distr_param(paramMap["Lmax_OBDM"], paramMap["mgr_dens"]);
    CorFunParam mom_distr_param(paramMap["kmakx"], paramMap["numks"]);

    Locals coordinates(param_model);
    Locals force(param_model);

    Energy **elocal;

    WaveFunction wave_func;

    DistributionR pair_distr(pair_distr_param);
    DistributionR pair_distr_cross(pair_distr_param);
    DistributionR dens_distr(dens_distr_param);
    OBDM obdm(obdm_param);

    double nkuppt;

    MomentDistr moment_distr(mom_distr_param);

    double *kr;

    int numks_temp = mom_distr_param.num_points;
    kr = new double[numks_temp];

    int ngr;

    elocal = new Energy*[dmnpop];
    for(int i = 0; i < dmnpop; i++) elocal[i] = new Energy[2];

    if(init == 1)
        coordinates.ReadInitial(startingConfig);
    else
        coordinates.GenerateInitial(startingConfig);


    ntemps = 0;

    in = 0;
    io = 1;

    nacc = 0;
    nprova = 0;

    force.SetZero();

    for(int iblck = 0; iblck < nblck; iblck++)
    {
        ngr = 0;

        pair_distr.SetZero();
        pair_distr_cross.SetZero();
        obdm.SetZero();
        dens_distr.SetZero();
        moment_distr.SetZero();

        nkuppt = 0;

        for(int iter = 0; iter < niter; iter++)
        {

            ntemps = ntemps + 1;

            eav.SetZero();

            jpop = 0;

            for(int ipop = 0; ipop < param_model.num_walk; ipop++)
            {

                force.metrop.SetZero();

                emtnew.SetZero();

                eold.tot = elocal[ipop][in].tot; // Energy::SetOldConf
                eold.kin = elocal[ipop][in].kin;
                eold.pot = elocal[ipop][in].pot;

                pair_distr.SetZeroAx();
                pair_distr_cross.SetZeroAx();
                dens_distr.SetZeroAx();
                obdm.SetZeroAx();
                moment_distr.SetZeroAx();

                coordinates.GaussianJump(ntemps, i_VMC, ipop, force);

                wave_func.Calc(param_model, coordinates);

                Energy_calc(coordinates.metrop, force.metrop, emtnew, param_model);


                if(i_Drift == 0)
                {
                     // This function is ugly, too many parameters ...
                     MetropolisDif(ipop, param_model, wave_func, coordinates, force, accepta, nprova, ntemps, in, i_VMC);
                }
                if(i_Drift == 1)
                {
                    accepta = 1;
                }

                if(accepta == 1)
                {
                    if(ntemps > 1) nacc = nacc + 1;

                    wave_func.Accept();

                    enew.tot = emtnew.tot; //Accept method of Energy
                    enew.kin = emtnew.kin;
                    enew.pot = emtnew.pot;

                    coordinates.Accept();

                    force.Accept();

                    // this will be converted to virtual call to "ObserveAuxilaryState" in
                    // coordinates.ObserveAuxilaryState(pair_distr, param_model);
                    // coordinates.ObserveAuxilaryState(pair_distr_cross, param_model);
                    pair_distr.PairDistrFirst(coordinates.auxil, param_model);
                    pair_distr_cross.PairDistrCross(coordinates.auxil, param_model);
                    dens_distr.DensityFirst(coordinates.auxil, param_model);

                    if(i_OBDM == 1)
                        obdm.OBDM_Calc(param_model, coordinates.auxil, wave_func, moment_distr, mom_distr_param);

                }
                else
                {

                    wave_func.NotAccept();

                    enew.tot = eold.tot; // NotAccept method of energy
                    enew.tot = eold.tot;
                    enew.tot = eold.tot;

                    coordinates.NotAccept(ipop);
                    force.NotAccept(ipop);

                    pair_distr.NotAccept(ipop, io);
                    pair_distr_cross.NotAccept(ipop, io);
                    dens_distr.NotAccept(ipop, io);
                    obdm.NotAccept(ipop);
                    moment_distr.NotAccept(ipop);

                }

                if (i_FNDMC == 1)               
                    BranchingCalc(&nsons, accepta, ntemps, nacc, nprova, param_model.alfa/4.0, nwalk_mean, ewalk, enew.tot, eold.tot, &param_model.seed, param_model.num_walk);


                else {nsons = 1;}

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

                        wave_func.WalkerMatch(jpop, io);

                        coordinates.WalkerMatch();
                        force.WalkerMatch();

                        pair_distr.WalkerMatch(jpop, io);
                        pair_distr_cross.WalkerMatch(jpop, io);
                        dens_distr.WalkerMatch(jpop, io);
                        obdm.WalkerMatch(jpop);
                        moment_distr.WalkerMatch(jpop);

                        jpop = jpop + 1;
                    }
                }

                eav.tot = eav.tot + enew.tot * nsons; // WalkerCollect for Energy

                ngr = ngr + nsons;

                pair_distr.WalkerCollect(nsons);
                pair_distr_cross.WalkerCollect(nsons);
                dens_distr.WalkerCollect(nsons);
                obdm.WalkerCollect(nsons);
                moment_distr.WalkerCollect(nsons);

                nkuppt = nkuppt + nsons * 100; //Number of McMillan points

            }

            param_model.num_walk = jpop;

            ii = in;
            in = io;
            io = ii;

            coordinates.PageSwap();
            force.PageSwap();

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
        //********************************************************************************************



        pair_distr.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);
        pair_distr_cross.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);

        // It should be CorFun::Print

        pair_distr.PrintDistr( outDir + "/gr.dat");
        pair_distr_cross.PrintDistr(outDir + "/gr_12.dat");

        dens_distr.NormalizationNR(ngr, dens_distr_param.step);
        dens_distr.PrintDistr(outDir + "/nr.dat");

        obdm.Normalization();
        obdm.PrintDistr(outDir + "/fr.dat");

        moment_distr.Normalization(param_model.num_part, nkuppt);
        moment_distr.PrintDistr(outDir + "/nk.dat");

        //***************************************************************************

        fstream outfile2(startingConfig, fstream::out ); //Coordinate::Print
        outfile2<<param_model.num_walk<<"\n";
        outfile2<<setprecision(18);
        for(int i = 0; i < param_model.num_walk;i++)
        {
            outfile2<<elocal[i][in].tot<<"\n";
            for(int ic = 0; ic < param_model.num_comp; ic++)
            {
                for(int ip = 0; ip < param_model.num_part; ip++ )

                    {outfile2<<coordinates.oldPage[i].GetParticleComp(ic,ip)<<"\n"; }
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

        int n_notused,  nfull, ncol, nel;

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
        statistics_cor(outDir + "/nk.dat", outDir + "/nk_av.dat", mom_distr_param.num_points, ncol_cor, n_notused_cor, nel_cor);

        ncol_cor = 1; //number of columns for 'gr.dat'
        statistics_cor(outDir + "/gr.dat", outDir + "/gr_av.dat", pair_distr_param.num_points, ncol_cor, n_notused_cor, nel_cor);

        ncol_cor = 1; //number of columns for 'gr.dat'
        statistics_cor(outDir + "/gr_12.dat", outDir + "/gr_12_av.dat", pair_distr_param.num_points, ncol_cor, n_notused_cor, nel_cor);

        ncol_cor = 1; //number of columns for 'fr.dat'
        statistics_cor(outDir + "/fr.dat", outDir + "/fr_av.dat", obdm_param.num_points, ncol_cor, n_notused_cor, nel_cor);

        ncol_cor = 1; //number of columns for 'nr.dat'
        statistics_cor(outDir + "/nr.dat", outDir + "/nr_av.dat", dens_distr_param.num_points, ncol_cor, n_notused_cor, nel_cor);


        Normalization(outDir + "/nk_av.dat", outDir + "/nk_av_norm.dat", mom_distr_param.num_points, param_model.num_part, param_model.num_comp);
        Normalization(outDir + "/nr_av.dat", outDir + "/nr_av_norm.dat", dens_distr_param.num_points, param_model.num_part, param_model.num_comp);

    }

    for(int i = 0; i < dmnpop; i++)
        delete [] elocal[i];
    delete [] elocal;

    delete [] kr;
}
