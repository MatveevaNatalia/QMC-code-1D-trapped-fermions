#include "utils.h"
#include <iomanip>
#include <fstream>
using namespace std;

map<string, double> FillMap(const string & inFile, const string & outDir){
    map<string, double> paramMap;
    fstream filestr(inFile, fstream::in | fstream::out);
    fstream outfile(outDir+"/Output.dat", fstream::out| ios::app );
    outfile<<setprecision(18);
    while(!filestr.eof()){
        string str;
        double value;
        filestr >> str >> value;
        paramMap[str] = value;
        outfile<<str <<" "<<paramMap[str]<<endl;
    }
    outfile.close();
    return paramMap;
}

void PrintAcceptance(double accrate, int iblck, const string & outDir){
    fstream outfile(outDir + "/Accept.dat", fstream::out| ios::app );
    outfile<< iblck <<" "<<accrate<<"\n";
    outfile.close();
}


void Normalization(string name_file, string name_file_1, int num, int np, int nc)
{
    double * kk, *cor_fun, *cor_fun_err;
    double sum, dk, norm;

    kk = new double[num];
    cor_fun = new double[num];
    cor_fun_err = new double[num];

    fstream filestr(name_file, fstream::in | fstream::out);
    for (int i = 0; i< num; i++ )
    {
        filestr >>setprecision(12)>> kk[i]>>cor_fun[i]>>cor_fun_err[i];
    }
    filestr.close();


    dk = kk[1]-kk[0];
    sum = 0;
    for(int i = 0; i < num; i++)
    {
        if(cor_fun[i]>0)
            sum = sum + dk * cor_fun[i];
    }

    norm = (np*nc)/(2.0*sum);//Here normalization on the number of particles per component

    fstream outfile(name_file_1, fstream::out );
    for(int i = 0; i < num; i++)
    {
        outfile<<kk[i]<<" "<<cor_fun[i]*norm<<" "<<cor_fun_err[i]*norm<<"\n";
    }
    outfile.close();

    delete [] kk;
    delete [] cor_fun;
    delete [] cor_fun_err;
}
