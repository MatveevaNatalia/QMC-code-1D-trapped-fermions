#include "Normalization.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
using namespace std;

void Normalization(const char * name_file, const char * name_file_1, const char * name_file_2, int num, int np, int nc)
{

double * kk, *cor_fun;
double temp, sum, dk, norm;

kk = new double[num];
cor_fun = new double[num];
cor_fun_err = new double[num]; 

fstream filestr(name_file, fstream::in | fstream::out);
for (int i = 0; i< num; i++ )
{
	filestr >>setprecision(12)>> kk[i]>>cor_fun[i]>>cor_fun_err[i];
//	cout<<kk[i]<<" "<<mom_distr_av[i]<<"\n";
}
filestr.close();


dk = cor_fun[1]-cor_fun[0];
sum = 0;
for(int i = 0; i < num; i++)
{
	sum = sum + dk * cor_fun[i];
}

norm = (np*nc)/(2.0*sum);//Here normalization at the total number of particles; factor of 2 is because the integral must be from -inf to +inf

fstream outfile(name_file_1, fstream::out| ios::app );
for(int i = 0; i < num; i++)
{
    outfile<<kk[i]<<cor_fun[i]*norm<<cor_fun_err[i]*norm<<"\n";
}
outfile.close();		
cout<<" Normalization for "<<name_file<<" = "<<norm<<"sum= "<<sum<<"\n";

delete [] kk;
delete [] cor_fun;
delete [] cor_fun_err;
}
