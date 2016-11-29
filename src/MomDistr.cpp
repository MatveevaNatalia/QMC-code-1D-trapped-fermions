#include "MomDistr.h"
using namespace std;

MomentDistr::MomentDistr(const CorFunParam& param):
    average(param.num_points),
    auxil(param.num_points),
    oldPage(dmnpop, ConfigDistr(param.num_points))
{
    num_points = param.num_points;
    step = param.step;
}

void MomentDistr::SetZero(){
    average.setzero();
}

void MomentDistr::SetZeroAx(){
    auxil.setzero();
}

void MomentDistr::NotAccept(int ipop){
    auxil = oldPage[ipop];
}

void MomentDistr::WalkerMatch(){
    newPage.push_back(auxil);
}

void MomentDistr::WalkerCollect(int nsons){
    for(int i = 1; i < (num_points+1); i++)
        average[i] += nsons * auxil[i];
}

void MomentDistr::Normalization(const ParamModel& param, int norm){
    for(int i = 0; i < num_points;i++)
        average[i] = average[i]*param.num_part/(float(norm));
}

void MomentDistr::PrintDistr(const string& name_file)
{
    double coord_bin;
    fstream outfile(name_file, fstream::out| ios::app );
    for(int i = 1; i < (num_points+1);i++)
    {
        coord_bin = float(i) * step + step/2.0;
        outfile<<coord_bin<<" "<<average[i]<<"\n";
    }
    outfile.close();
}

void MomentDistr::PageSwap(){
     oldPage = newPage;
     newPage.clear();
}


