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

void run (const string & inFile, const string & startingConfig, const string & outDir )
{
    int ntemps, nacc, nprova, nwrite;
    double acept_rate;
    int jpop, nsons, npopmean;
    int nblck, niter;
    bool accept, init, i_VMC, i_Drift, i_FNDMC, i_OBDM;

    map<string, double> paramMap;
    paramMap = FillMap(inFile, outDir);

    i_VMC   = paramMap["i_VMC"];
    i_Drift = paramMap["i_Drift"];
    i_FNDMC = paramMap["i_FNDMC"];
    niter = paramMap["niter"];
    nblck = paramMap["nblck"];
    init = paramMap["init"];
    nwrite = paramMap["nwrite"];
    i_OBDM = paramMap["i_OBDM"];

    ParamModel param_model(paramMap);
    CorFunParam pair_distr_param(paramMap["Lmax_gr"], paramMap["mgr_g(r)"]);
    CorFunParam obdm_param(paramMap["Lmax_OBDM"], paramMap["mgr_OBDM"]);
    CorFunParam dens_distr_param(paramMap["Lmax_OBDM"], paramMap["mgr_dens"]);
    CorFunParam mom_distr_param(paramMap["kmakx"], paramMap["numks"]);

    Locals coordinates(param_model);
    Locals force(param_model);

    Energy energy;
    WaveFunction wave_func;

    DistributionR pair_distr(pair_distr_param);
    DistributionR pair_distr_cross(pair_distr_param);
    DistributionR dens_distr(dens_distr_param);
    OBDM obdm(obdm_param);
    MomentDistr moment_distr(mom_distr_param);

    double n_rand_points;
    int nsons_total;

    if(init)
        coordinates.ReadInitial(startingConfig);
    else
        coordinates.GenerateInitial(startingConfig);

    ntemps = 0;
    nacc = 0;
    nprova = 0;
    npopmean = 0;

    nsons_total = 0;
    n_rand_points = 0;

    for(int iter = 0; iter < niter*nblck; iter++)
    {
        ntemps++;
        energy.SetZeroAverage();
        jpop = 0;

        for(int ipop = 0; ipop < param_model.num_walk; ipop++)
        {

            // Same method
            force.metrop.SetZero();
            energy.SetZeroMetrop();

            pair_distr.SetZeroAx();
            pair_distr_cross.SetZeroAx();
            dens_distr.SetZeroAx();
            obdm.SetZeroAx();
            moment_distr.SetZeroAx();
            // ----

            coordinates.GaussianJump(ntemps, i_VMC, ipop, force);
            wave_func.Calc(param_model, coordinates);
            energy.Calc(coordinates, force.metrop, param_model);

            if(i_Drift)
                accept = true;
            else
                accept = MetropolisDif(ipop, param_model, wave_func, coordinates, force, nprova, ntemps, i_VMC);

            if(accept)
            {
                if(ntemps > 1)
                    nacc ++;

                // Same method
                wave_func.Accept();
                energy.Accept();
                coordinates.Accept();
                force.Accept();
                // ---

                // Same method
                pair_distr.PairDistrFirst(coordinates.auxil, param_model);
                pair_distr_cross.PairDistrCross(coordinates.auxil, param_model);
                dens_distr.DensityFirst(coordinates.auxil, param_model);
                // ----

                if(i_OBDM)
                    obdm.OBDM_Calc(param_model, coordinates.auxil, wave_func, moment_distr, mom_distr_param);
            }
            else
            {

                // Same method
                wave_func.NotAccept(ipop);
                energy.NotAccept(ipop);
                coordinates.NotAccept(ipop);
                force.NotAccept(ipop);
                pair_distr.NotAccept(ipop);
                pair_distr_cross.NotAccept(ipop);
                dens_distr.NotAccept(ipop);
                obdm.NotAccept(ipop);
                moment_distr.NotAccept(ipop);
                // ----

            }

            if (i_FNDMC)
                nsons = BranchingCalc(param_model, energy, accept, ntemps, nacc, nprova, ipop);
            else
                nsons = 1;

            if(nsons > 0)
            {
                for(int js = 0; js < nsons; js++)
                {

                    if(jpop >= dmnpop)
                        cout <<"The population has grown too much!"<<"nsons="<<nsons<<"jpop="<<jpop<<"\n";

                    energy.WalkerMatch();
                    wave_func.WalkerMatch();
                    coordinates.WalkerMatch();
                    force.WalkerMatch();
                    pair_distr.WalkerMatch();
                    pair_distr_cross.WalkerMatch();
                    dens_distr.WalkerMatch();
                    obdm.WalkerMatch();
                    moment_distr.WalkerMatch();

                    jpop++;
                }
            }

            energy.WalkerCollect(nsons);

            nsons_total += nsons;

            pair_distr.WalkerCollect(nsons);
            pair_distr_cross.WalkerCollect(nsons);
            dens_distr.WalkerCollect(nsons);
            obdm.WalkerCollect(nsons);
            moment_distr.WalkerCollect(nsons);

            n_rand_points += nsons * 100; //Number of McMillan points
        }
        param_model.num_walk = jpop;

        coordinates.PageSwap();
        force.PageSwap();
        wave_func.PageSwap();
        pair_distr.PageSwap();
        pair_distr_cross.PageSwap();
        dens_distr.PageSwap();
        obdm.PageSwap();
        moment_distr.PageSwap();
        energy.PageSwap();

        energy.Normalization(param_model, jpop);

        energy.Average();
        npopmean = npopmean + param_model.num_walk;

        if(iter % nwrite == 0)
        {
            energy.Print(nwrite, npopmean, ntemps, outDir);
            energy.SetZeroMean();
            npopmean = 0;
        }

        if(iter % niter == 0)
        {

            pair_distr.NormalizationGR(nsons_total, param_model);
            pair_distr_cross.NormalizationGR(nsons_total, param_model);

            pair_distr.PrintDistr( outDir + "/gr.dat");
            pair_distr_cross.PrintDistr(outDir + "/gr_12.dat");

            dens_distr.NormalizationNR(nsons_total);
            dens_distr.PrintDistr(outDir + "/nr.dat");

            obdm.Normalization();
            obdm.PrintDistr(outDir + "/fr.dat");

            moment_distr.Normalization(param_model, n_rand_points);
            moment_distr.PrintDistr(outDir + "/nk.dat");

            coordinates.PrintAll(startingConfig);

            acept_rate = 100.0 * float(nacc)/float(nprova);
            PrintAcceptance(acept_rate, iter/niter, outDir);

            pair_distr.SetZero();
            pair_distr_cross.SetZero();
            obdm.SetZero();
            dens_distr.SetZero();
            moment_distr.SetZero();
            nsons_total = 0;
            n_rand_points = 0;
        }
    }
}
