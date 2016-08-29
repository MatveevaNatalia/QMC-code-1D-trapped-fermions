
#include "Energy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

//void Energy_calc(double **xaux, double **FT, double *Etot, double *Epot_trap, double *Ekin,  int num_comp, int np, double aB, double a, double width)
//void Energy_calc(double **xaux, double **FT, Energy& energy,  int num_comp, int np, double aB, double a, double width)


void Energy_calc(double **xMT, double** FMT, Energy& energy, ParamModel param_mode){
        double Ekin_J1 = 0, Ekin_J2 = 0, Ekin_G = 0, Ekin_GJ1 = 0;
        double Ekin_GJ2 = 0, Ekin_J1J2 = 0, Epot = 0, Ekin_total = 0;
        double xi, xk;
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

        //Forces (WITHOUT factor of 2)

        for(int gamma = 0; gamma < num_comp; gamma++)
        {
            for(int k = 0; k < num_part; k++)
            {
                xk = xMT[gamma][k];

                FG[gamma][k] = -2.0*width*xk;

                for(int i = 0; i < num_part; i++)
                {
                    if(i != k) {
                        xi = xMT[gamma][i];
                        FJ1[gamma][k] +=  ForcePartial(xk, xi, scat_length_bos);
                    }
                }

                for(int alpha = 0; alpha < num_comp; alpha++)
                {
                    if(alpha != gamma)
                    {
                        for(int i  = 0; i < num_part; i++)
                        {
                            xi = xMT[alpha][i];
                            FJ2[gamma][k] += ForcePartial(xk, xi, scat_length);
                        }
                    }
                }

                FMT[gamma][k] = 2.0*(FG[gamma][k] + FJ1[gamma][k] + FJ2[gamma][k]);

                //	cout<<"FG: "<<"\n";
                //	cout<<"gamma= "<<gamma<<" k= "<<k<<" "<<FG[gamma][k];
                //  cout<<" "<< FJ1[gamma][k]<<" "<<FJ2[gamma][k]<<" "<<FT[gamma][k]<<"\n";

            }
        }


        for(int gamma = 0; gamma < num_comp; gamma++)
        {
            for(int k = 0; k < num_part; k++)
            {
                xk = xMT[gamma][k];

                Ekin_J1 = Ekin_J1 - 0.5*FJ1[gamma][k]*FJ1[gamma][k];

                for(int i = 0; i < num_part; i++)
                {
                    if( i!= k){
                        xi = xMT[gamma][i];
                        Ekin_J1 += EnergyPartial(xk, xi, scat_length_bos);}
                }

                Ekin_J2 = Ekin_J2 - 0.5 * FJ2[gamma][k]*FJ2[gamma][k];

                for(int alpha = 0; alpha < num_comp; alpha++)
                {
                    if(alpha != gamma)
                    {
                        for(int i = 0; i < num_part; i++)
                        {
                            xi = xMT[alpha][i];
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
        energy.tot = Ekin_total + Epot;

        for(int gamma = 0; gamma < num_comp; gamma++)
        {
            for(int k = 0; k < num_part; k++)
                sum = sum + 0.25*FMT[gamma][k]*FMT[gamma][k];
        }

        energy.kin = 0.5*sum;

        energy.pot = Epot;

        //cout<<"Ekin_J1 = "<<Ekin_J1<<" Ekin_J2= "<<Ekin_J2<<" Ekin_G= "<<Ekin_G;
        //cout<<" Ekin_GJ1= "<< Ekin_GJ1<<" Ekin_GJ2= "<<Ekin_GJ2;
        //cout<<" Ekin_J1J2= "<<Ekin_J1J2<<" Epot_trap = "<<Epot<<"\n";

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
