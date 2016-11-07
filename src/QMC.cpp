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
//#include "metrop.h"
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
//
//                               .-----.
//                              /7  .  (
//                             /   .-.  \
//                            /   /   \  \
//                           / `  )   (   )
//                          / `   )   ).  \
//                        .'  _.   \_/  . |
//       .--.           .' _.' )`.        |
//      (    `---...._.'   `---.'_)    ..  \
//       \            `----....___    `. \  |
//        `.           _ ----- _   `._  )/  |
//          `.       /"  \   /"  \`.  `._   |
//            `.    ((O)` ) ((O)` ) `.   `._\
//              `-- '`---'   `---' )  `.    `-.
//                 /                  ` \      `-.
//               .'                      `.       `.
//              /                     `  ` `.       `-.
//       .--.   \ ===._____.======. `    `   `. .___.--`     .''''.
//      ' .` `-. `.                )`. `   ` ` \          .' . '  8)
//     (8  .  ` `-.`.               ( .  ` `  .`\      .'  '    ' /
//      \  `. `    `-.               ) ` .   ` ` \  .'   ' .  '  /
//       \ ` `.  ` . \`.    .--.     |  ` ) `   .``/   '  // .  /
//        `.  ``. .   \ \   .-- `.  (  ` /_   ` . / ' .  '/   .'
//          `. ` \  `  \ \  '-.   `-'  .'  `-.  `   .  .'/  .'
//            \ `.`.  ` \ \    ) /`._.`       `.  ` .  .'  /
//             |  `.`. . \ \  (.'               `.   .'  .'
//          __/  .. \ \ ` ) \                     \.' .. \__
//   .-._.-'     '"  ) .-'   `.                   (  '"     `-._.--.
//  (_________.-====' / .' /\_)`--..__________..-- `====-. _________)
//                   (.'(.'
//



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

    CorFunParam pair_distr_param, obdm_param, dens_distr_param, mom_distr_param;

    double Lmax;

    int in, io, ii, jpop, nsons, npopmean;

    ParamModel param_model;
    param_model.seed = kkk;

    int nblck, niter, init, icrit;

    string str;

    int i_VMC, i_Smart = 0, i_Drift = 0, i_FNDMC, i_OBDM_1, i_OBDM_2=0, i_stat_cor;

    fstream filestr2(inFile, fstream::in | fstream::out);
    filestr2>>str>>i_VMC>>str>>i_FNDMC>>str>>param_model.num_comp;
    filestr2>>str>>param_model.num_part>>str>>param_model.width;
    filestr2>>str>>param_model.scat_lenght_bos>>str>>param_model.scat_lenght;
    filestr2>>str>>param_model.alfa>>str>>niter>>str>>nblck>>str;
    filestr2>>param_model.num_walk>>str>>init>>str>>nwrite>>str>>icrit;
    filestr2>>str>>i_OBDM_1>>str>>pair_distr_param.num_points>>str>>obdm_param.num_points;
    filestr2>>str>>dens_distr_param.num_points>>str>>pair_distr_param.max_value;
    filestr2>>str>>obdm_param.max_value>>str>>dens_distr_param.max_value;
    filestr2>>str>>mom_distr_param.num_points>>str>>mom_distr_param.max_value>>str>>Lmax;
    filestr2.close();


    Locals coordinates(param_model);
    Locals force(param_model);

    pair_distr_param.step = pair_distr_param.max_value/pair_distr_param.num_points;

    obdm_param.step = obdm_param.max_value/obdm_param.num_points;
    dens_distr_param.step = dens_distr_param.max_value/dens_distr_param.num_points;

    mom_distr_param.step = mom_distr_param.max_value/mom_distr_param.num_points;

    fstream outfile2(outDir+"/Output.dat", fstream::out| ios::app );
    outfile2<<setprecision(18);
    outfile2<<" i_VMC= "<<i_VMC<<" i_FNDMC= "<<i_FNDMC<<" ncomp= "<<param_model.num_comp;
    outfile2<<" np= "<<param_model.num_part<<" width= "<<param_model.width;
    outfile2<<" aB= "<<param_model.scat_lenght_bos<<" a= "<<param_model.scat_lenght<<"\n";
    outfile2<<" alfa= "<<param_model.alfa<<" niter= "<<niter<<" nblck= ";
    outfile2<<nblck<<" npop= "<<param_model.num_walk;
    outfile2<<" init= "<<init<<" nwrite= "<<nwrite<<" icrit= "<<icrit<<"\n";
    outfile2<<" i_OBDM= "<<i_OBDM_1<<"\n";
    outfile2<<" mgr_g(r)= "<<pair_distr_param.num_points<<" mgr_OBDM= "<<obdm_param.num_points;
    outfile2<<" mgr_dens= "<<dens_distr_param.num_points<<"\n";
    outfile2<<" Lmax_g(r)= "<<pair_distr_param.max_value<<" Lmax_OBDM= "<<obdm_param.max_value;
    outfile2<<" Lmax_dens= "<<dens_distr_param.max_value<<"\n";
    outfile2<<" numks= "<<mom_distr_param.num_points<<" kmax= "<<mom_distr_param.max_value;
    outfile2<<" Lmax_McM= "<<Lmax<<"\n";

    outfile2.close();
    Energy **elocal;
    double **flocal;



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

        pair_distr.SetZero();
        pair_distr_cross.SetZero();
        obdm.SetZero();
        dens_distr.SetZero();
        moment_distr.SetZero();

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

                pair_distr.SetZeroAx();
                pair_distr_cross.SetZeroAx();
                dens_distr.SetZeroAx();
                obdm.SetZeroAx();
                moment_distr.SetZeroAx();

                coordinates.GaussianJump(ntemps, in, i_VMC, ipop, force.total);

                PsiTotal = WaveFunction(param_model, coordinates);

                Energy_calc(coordinates.metrop, force.metrop, emtnew, param_model);

                pair_distr.PairDistrFirst(coordinates.metrop, param_model);
                pair_distr_cross.PairDistrCross(coordinates.metrop, param_model);
                dens_distr.DensityFirst( coordinates.metrop, param_model);


                if(i_Drift == 0)
                {
                     MetropolisDif(ipop, param_model, PsiTotal, flocal, coordinates, force, accepta, nprova, fvella, ntemps, in, i_VMC);
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

                    pair_distr.Accept();
                    pair_distr_cross.Accept();

                   // cout <<ntemps <<" "<<pair_distr_1.dra[5] << endl;

                    dens_distr.Accept();

                    if(i_OBDM_1 == 1)
                    {

                        //OBDM1D_11(Lmax,ncomp,np,width, aB, a, kr, xaux,fra,nfra,PsiTotal, dnkupa, numks, dk, &kkk, mgr2, dr2);
                        OBDM1D_11(Lmax,param_model.num_comp,param_model.num_part,param_model.width, param_model.scat_lenght_bos, param_model.scat_lenght, kr, coordinates.auxil, obdm.fra, obdm.nfra, PsiTotal, moment_distr.dnkupa, mom_distr_param.num_points, mom_distr_param.step, &param_model.seed, obdm_param.num_points, obdm_param.step);

                    }

                }
                else
                {
                    fnew = fvella; // NotAccept method of WaveFunction

                    enew.tot = eold.tot; // NotAccept method of energy
                    enew.tot = eold.tot;
                    enew.tot = eold.tot;

                    coordinates.NotAccept(ipop, in);
                    force.NotAccept(ipop, in);

                    pair_distr.NotAccept(ipop, in);
                    pair_distr_cross.NotAccept(ipop, in);
                    dens_distr.NotAccept(ipop, in);
                    obdm.NotAccept(ipop);
                    moment_distr.NotAccept(ipop);

                }

                if (i_FNDMC == 1)
                {
                    BranchingCalc(&nsons, accepta, ntemps, nacc, nprova, param_model.alfa/4.0, icrit, ewalk, enew.tot, eold.tot, &param_model.seed, param_model.num_walk);
                }

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

                        flocal[jpop][io] = fnew;

                        coordinates.WalkerMatch(jpop, io);
                        force.WalkerMatch(jpop, io);

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

        double r1;

        pair_distr.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);
        pair_distr_cross.NormalizationGR(ngr, param_model.num_comp, param_model.num_part, pair_distr_param.step);

        // It should be CorFun::Print

        pair_distr.PrintDistr( outDir + "/gr.dat");
        pair_distr_cross.PrintDistr(outDir + "/gr_12.dat");

        dens_distr.NormalizationNR(ngr, dens_distr_param.step);
        dens_distr.PrintDistr(outDir + "/nr.dat");

        obdm.Normalization();
        obdm.PrintDistr(outDir + "/fr.dat");


        //	double n0 =  dnkup[0]*np/float(nkuppt);
        //	double n0_22 =  dnkup_22[0]*np/float(nkuppt);

        double k1;

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
                    {outfile2<<coordinates.total[ic][ip][i][in]<<"\n"; }
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

    delete [] kr;
}
