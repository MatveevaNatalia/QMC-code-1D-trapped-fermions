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
    int ntemps, accepta, nacc, nprova, nwrite, nmean;
    double accrate;
    int in, io, ii, jpop, nsons, npopmean;
    int nblck, niter, init;
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

    double nkuppt;
    double *kr;

    int numks_temp = mom_distr_param.num_points;
    kr = new double[numks_temp];
    int ngr;

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

            energy.SetZeroAverage();

            jpop = 0;

            for(int ipop = 0; ipop < param_model.num_walk; ipop++)
            {

                force.metrop.SetZero();
                energy.SetZeroMetrop();
                energy.SetOldConf(ipop,in);

                pair_distr.SetZeroAx();
                pair_distr_cross.SetZeroAx();
                dens_distr.SetZeroAx();
                obdm.SetZeroAx();
                moment_distr.SetZeroAx();

                coordinates.GaussianJump(ntemps, i_VMC, ipop, force);
                wave_func.Calc(param_model, coordinates);
                energy.Calc(coordinates, force.metrop, param_model);

                if(i_Drift == 0)                                    
                     MetropolisDif(ipop, param_model, wave_func, coordinates, force, accepta, nprova, ntemps, in, i_VMC);

                if(i_Drift == 1)                
                    accepta = 1;

                if(accepta == 1)
                {
                    if(ntemps > 1) nacc = nacc + 1;

                    wave_func.Accept();
                    energy.Accept();
                    coordinates.Accept();
                    force.Accept();

                    pair_distr.PairDistrFirst(coordinates.auxil, param_model);
                    pair_distr_cross.PairDistrCross(coordinates.auxil, param_model);
                    dens_distr.DensityFirst(coordinates.auxil, param_model);

                    if(i_OBDM == 1)
                        obdm.OBDM_Calc(param_model, coordinates.auxil, wave_func, moment_distr, mom_distr_param);
                }
                else
                {

                    wave_func.NotAccept();
                    energy.NotAccept();
                    coordinates.NotAccept(ipop);
                    force.NotAccept(ipop);
                    pair_distr.NotAccept(ipop, io);
                    pair_distr_cross.NotAccept(ipop, io);
                    dens_distr.NotAccept(ipop, io);
                    obdm.NotAccept(ipop);
                    moment_distr.NotAccept(ipop);

                }

                if (i_FNDMC == 1)               
                    nsons = BranchingCalc(param_model, energy, accepta, ntemps, nacc, nprova);
                else
                    nsons = 1;

                if(nsons > 0)
                {
                    for(int js = 0; js < nsons; js++)
                    {

                        if(jpop >= dmnpop)                       
                            cout <<"The population has grown too much!"<<"nsons="<<nsons<<"jpop="<<jpop<<"\n";

                        energy.WalkerMatch(jpop, io);
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

                energy.WalkerCollect(nsons);

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

            energy.Normalization(param_model, jpop);

            if(ntemps == 1)
            {            
                energy.SetZeroMean();
                npopmean = 0;
                nmean = 0;
            }

            nmean = nmean + 1;
            energy.Average();

            npopmean = npopmean + param_model.num_walk;

            if(nmean==nwrite)
            {        
                energy.Print(nwrite, npopmean, ntemps, outDir);
                energy.SetZeroMean();
                nmean = 0;
                npopmean = 0;
            }
        }

        pair_distr.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);
        pair_distr_cross.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);

        pair_distr.PrintDistr( outDir + "/gr.dat");
        pair_distr_cross.PrintDistr(outDir + "/gr_12.dat");

        dens_distr.NormalizationNR(ngr, dens_distr_param.step);
        dens_distr.PrintDistr(outDir + "/nr.dat");

        obdm.Normalization();
        obdm.PrintDistr(outDir + "/fr.dat");

        moment_distr.Normalization(param_model.num_part, nkuppt);
        moment_distr.PrintDistr(outDir + "/nk.dat");

        coordinates.PrintAll(startingConfig, in);

        accrate = 100.0 * float(nacc)/float(nprova);
        PrintAcceptance(accrate, iblck, outDir);
    }
    delete [] kr;
}
