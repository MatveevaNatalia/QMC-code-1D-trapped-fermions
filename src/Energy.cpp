#include "Energy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void EnergyVector::SetZero(){
    tot = 0; pot = 0; kin = 0;
};

Energy::Energy(){
    elocal = new EnergyVector*[dmnpop];
    for(int i = 0; i < dmnpop; i++) elocal[i] = new EnergyVector[2];
    ewalk = 0.;
}

void Energy::SetOldConf(int ipop, int in){
    eold.tot = elocal[ipop][in].tot;
    eold.kin = elocal[ipop][in].kin;
    eold.pot = elocal[ipop][in].pot;
}

void Energy::Accept(){
    enew.tot = emetrop.tot;
    enew.kin = emetrop.kin;
    enew.pot = emetrop.pot;
}

void Energy::NotAccept(){
    enew.tot = eold.tot;
    enew.tot = eold.tot;
    enew.tot = eold.tot;
}

void Energy::WalkerMatch(int jpop, int io){
    elocal[jpop][io].tot = enew.tot;
    elocal[jpop][io].kin = enew.kin;
    elocal[jpop][io].pot = enew.pot;
}

void Energy::Normalization(const ParamModel& param_model, int jpop){
    int num_comp, num_part;
    num_comp = param_model.num_comp;
    num_part = param_model.num_part;
    int norm = jpop * num_comp * num_part;

    ewalk = eav.tot/jpop;
    epar.tot = eav.tot/double(norm);
    epar.kin = eav.kin/double(norm);
    epar.pot = eav.pot/double(norm);
}

void Energy::Average(){
    emean.tot = emean.tot + epar.tot;
    emean.kin = emean.kin + epar.kin;
    emean.pot = emean.pot + epar.pot;
}

void Energy::Print(int nwrite, int npopmean, const string & outDir){
    fstream outfile(outDir + "/Energy.dat", fstream::out| ios::app );
    outfile<<setprecision(18);
    outfile<< ntemps/nwrite <<" "<<emean.tot/float(nwrite)<<" "<<emean.kin/float(nwrite)<<" "<<emean.pot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";
    cout<< ntemps/nwrite <<" "<<emean.tot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";
    outfile.close();
}



void Energy::Calc(Locals& coord, Configuration& FMT, const ParamModel& param_mode){
    double Ekin_J1 = 0, Ekin_J2 = 0, Ekin_G = 0, Ekin_GJ1 = 0;
    double Ekin_GJ2 = 0, Ekin_J1J2 = 0, Epot = 0, Ekin_total = 0;
    double xi, xk, force_temp;
    double **FG, **FJ1, **FJ2;
    double sum = 0;
    int num_comp, num_part;
    double scat_length_bos, scat_length, width;

    num_comp = param_mode.num_comp;
    num_part = param_mode.num_part;
    scat_length_bos = param_mode.scat_lenght_bos;
    scat_length = param_mode.scat_lenght;
    width = param_mode.width;

    FG = new double*[num_comp];
    for(int i = 0; i < num_comp; i++) FG[i] = new double[num_part];

    FJ1 = new double*[num_comp];
    for(int i = 0; i < num_comp; i++) FJ1[i] = new double[num_part];

    FJ2 = new double*[num_comp];
    for(int i = 0; i < num_comp; i++) FJ2[i] = new double[num_part];

    for(int alpha = 0; alpha < num_comp; alpha++)
    {
        for(int i = 0; i < num_part; i++)
        {
            FG[alpha][i] = 0.0;
            FJ1[alpha][i] = 0.0;
            FJ2[alpha][i] = 0.0;
        }
    }

    for(int gamma = 0; gamma < num_comp; gamma++)
    {
        for(int k = 0; k < num_part; k++)
        {        
            xk = coord.metrop.GetParticleComp(gamma,k);

            FG[gamma][k] = -2.0*width*xk;

            for(int i = 0; i < num_part; i++)
            {
                if(i != k) {                    
                    xi = xMT.GetParticleComp(gamma,i);
                    FJ1[gamma][k] +=  ForcePartial(xk, xi, scat_length_bos);
                }
            }

            for(int alpha = 0; alpha < num_comp; alpha++)
            {
                if(alpha != gamma)
                {
                    for(int i  = 0; i < num_part; i++)
                    {
                        xi = coord.metrop.GetParticleComp(alpha,i);

                        FJ2[gamma][k] += ForcePartial(xk, xi, scat_length);
                    }
                }
            }
            force_temp = 2.0*(FG[gamma][k] + FJ1[gamma][k] + FJ2[gamma][k]);
            FMT.SetParticleComp(gamma,k,force_temp);
        }
    }


    for(int gamma = 0; gamma < num_comp; gamma++)
    {
        for(int k = 0; k < num_part; k++)
        {          
            xk = coord.metrop.GetParticleComp(gamma,k);

            Ekin_J1 = Ekin_J1 - 0.5*FJ1[gamma][k]*FJ1[gamma][k];

            for(int i = 0; i < num_part; i++)
            {
                if( i!= k){                   
                    xi = coord.metrop.GetParticleComp(gamma,i);

                    Ekin_J1 += EnergyPartial(xk, xi, scat_length_bos);}
            }

            Ekin_J2 = Ekin_J2 - 0.5 * FJ2[gamma][k]*FJ2[gamma][k];

            for(int alpha = 0; alpha < num_comp; alpha++)
            {
                if(alpha != gamma)
                {
                    for(int i = 0; i < num_part; i++)
                    {                       
                        xi = xMT.GetParticleComp(alpha, i);
                        Ekin_J2 += EnergyPartial(xk, xi, scat_length);
                    }
                }
            }

            Ekin_G += 0.5*(2.0*width - FG[gamma][k]*FG[gamma][k]);

            Ekin_GJ1 -= FG[gamma][k]*FJ1[gamma][k];

            Ekin_GJ2 -= FG[gamma][k]*FJ2[gamma][k];

            Ekin_J1J2 -= FJ1[gamma][k]*FJ2[gamma][k];

            Epot = Epot + 0.5*xk*xk;
        }
    }

    Ekin_total = Ekin_J1 + Ekin_J2 + Ekin_G + Ekin_GJ1 + Ekin_GJ2 + Ekin_J1J2;
    emetrop.tot = Ekin_total + Epot;

    for(int gamma = 0; gamma < num_comp; gamma++)
    {
        for(int k = 0; k < num_part; k++)
            force_temp = FMT.GetParticleComp(gamma,k);
            sum = sum + 0.25 * force_temp * force_temp;
    }

    emetrop.kin = 0.5 * sum;

    emetrop.pot = Epot;

    for(int ic = 0; ic < num_comp; ic++)
        delete [] FG[ic];
    delete [] FG;

    for(int ic = 0; ic < num_comp; ic++)
        delete [] FJ1[ic];
    delete [] FJ1;

    for(int ic = 0; ic < num_comp; ic++)
        delete [] FJ2[ic];
    delete [] FJ2;
}

double ForcePartial(double xk, double xi, double scat_length){
    double Fki = 0, temp;
    temp = fabs(xk-xi)-scat_length;
    Fki = temp*(xk-xi)/(fabs(temp)*fabs(temp)*fabs(xk-xi));
    return Fki;
}

double EnergyPartial(double xk, double xi, double scat_length){
    double Eki;
    Eki = 0.5/(fabs(fabs(xk-xi)-scat_length)*fabs(fabs(xk-xi)-scat_length));
    return Eki;
}
