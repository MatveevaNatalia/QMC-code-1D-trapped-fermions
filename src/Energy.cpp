#include "Energy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
using namespace std;

void EnergyState::SetZero(){
    tot = 0; pot = 0; kin = 0;
};

Energy::Energy():
    oldPage(dmnpop, EnergyState())
{
    mean.SetZero();
}

void Energy::Accept(){
    auxil = metrop;
}

void Energy::NotAccept(int ipop){
    auxil = oldPage[ipop];
}

void Energy::WalkerMatch(){
    newPage.push_back(auxil);
}

void Energy::WalkerCollect(int nsons){
    average.tot = average.tot + auxil.tot * nsons;
}

void Energy::PageSwap(){
     oldPage = newPage;
     newPage.clear();
}

void Energy::Normalization(const ParamModel& param_model, int jpop){
    int num_comp, num_part;
    num_comp = param_model.num_comp;
    num_part = param_model.num_part;
    int norm = jpop * num_comp * num_part;

    average_walk = average.tot/jpop;
    average_part.tot = average.tot/double(norm);
    average_part.kin = average.kin/double(norm);
    average_part.pot = average.pot/double(norm);
}


void Energy::Average(){
    mean.tot = mean.tot + average_part.tot;
    mean.kin = mean.kin + average_part.kin;
    mean.pot = mean.pot + average_part.pot;
}

void Energy::Print(int nwrite, int npopmean, int ntemps, const string & outDir){
    fstream outfile(outDir + "/Energy.dat", fstream::out| ios::app );
    outfile<<setprecision(18);
    outfile<< ntemps/nwrite <<" "<<mean.tot/float(nwrite)<<" "<<mean.kin/float(nwrite);
    outfile<<" "<<mean.pot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";
    cout<< ntemps/nwrite <<" "<<mean.tot/float(nwrite)<<" "<<float(npopmean)/float(nwrite)<<"\n";
    outfile.close();
}

void Energy::Calc(Locals& coord, Configuration& FMT, const ParamModel& param_mode){
    double Ekin_J1 = 0, Ekin_J2 = 0, Ekin_G = 0, Ekin_GJ1 = 0;
    double Ekin_GJ2 = 0, Ekin_J1J2 = 0, Epot = 0, Ekin_total = 0;
    double xi, xk, force_temp;
    double sum = 0;
    int num_comp, num_part;
    double scat_length_bos, scat_length, width;

    num_comp = param_mode.num_comp;
    num_part = param_mode.num_part;
    scat_length_bos = param_mode.scat_lenght_bos;
    scat_length = param_mode.scat_lenght;
    width = param_mode.width;

    vector<vector<double>> FG(num_comp, vector<double>(num_part, 0.));
    vector<vector<double>> FJ1(num_comp, vector<double>(num_part, 0.));
    vector<vector<double>> FJ2(num_comp, vector<double>(num_part, 0.));

    for(int gamma = 0; gamma < num_comp; gamma++)
    {
        for(int k = 0; k < num_part; k++)
        {
            xk = coord.metrop.GetParticleComp(gamma,k);

            FG[gamma][k] = -2.0*width*xk;

            for(int i = 0; i < num_part; i++)
            {
                if(i != k) {
                    xi = coord.metrop.GetParticleComp(gamma,i);
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

            Ekin_J1 -= 0.5*FJ1[gamma][k]*FJ1[gamma][k];

            for(int i = 0; i < num_part; i++)
            {
                if( i!= k){
                    xi = coord.metrop.GetParticleComp(gamma,i);

                    Ekin_J1 += EnergyPartial(xk, xi, scat_length_bos);}
            }

            Ekin_J2 -= 0.5 * FJ2[gamma][k]*FJ2[gamma][k];

            for(int alpha = 0; alpha < num_comp; alpha++)
            {
                if(alpha != gamma)
                {
                    for(int i = 0; i < num_part; i++)
                    {
                        xi = coord.metrop.GetParticleComp(alpha, i);
                        Ekin_J2 += EnergyPartial(xk, xi, scat_length);
                    }
                }
            }

            Ekin_G += 0.5*(2.0*width - FG[gamma][k]*FG[gamma][k]);

            Ekin_GJ1 -= FG[gamma][k]*FJ1[gamma][k];

            Ekin_GJ2 -= FG[gamma][k]*FJ2[gamma][k];

            Ekin_J1J2 -= FJ1[gamma][k]*FJ2[gamma][k];

            Epot += 0.5*xk*xk;
        }
    }

    Ekin_total = Ekin_J1 + Ekin_J2 + Ekin_G + Ekin_GJ1 + Ekin_GJ2 + Ekin_J1J2;
    metrop.tot = Ekin_total + Epot;

    for(int gamma = 0; gamma < num_comp; gamma++)
    {
        for(int k = 0; k < num_part; k++)
            force_temp = FMT.GetParticleComp(gamma,k);
            sum += 0.25 * force_temp * force_temp;
    }

    metrop.kin = 0.5 * sum;
    metrop.pot = Epot;
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

