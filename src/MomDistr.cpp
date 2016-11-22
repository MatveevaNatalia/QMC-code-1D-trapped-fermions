#include "MomDistr.h"


using namespace std;

MomentDistr::MomentDistr(const CorFunParam& mom_distr){
    num_points = mom_distr.num_points;
    step = mom_distr.step;

    dnkup = new double[num_points];
    dnkupa = new double[num_points];

    dnkuplocal = new double*[num_points];
    for(int i = 0; i < num_points; i++) dnkuplocal[i] = new double[dmnpop];
}

void MomentDistr::SetZero(){
    setzero(dnkup);
}
void MomentDistr::SetZeroAx(){
    setzero(dnkupa);
}

void MomentDistr::NotAccept(int ipop){
    for(int i = 0; i < num_points; i++)
        dnkupa[i] = dnkuplocal[i][ipop];
}

void MomentDistr::WalkerMatch(int jpop){
    for(int i = 0; i < num_points; i++)
        dnkuplocal[i][jpop] = dnkupa[i];
}

void MomentDistr::WalkerCollect(int nsons){
    for(int i = 0; i < num_points; i++)
        dnkup[i] = dnkup[i] + dnkupa[i]*nsons;
}

void MomentDistr::Normalization(int np, int nkuppt){
    for(int ik = 0; ik < num_points;ik++)
        dnkup[ik] = dnkup[ik]*np/(float(nkuppt));
}


void MomentDistr::PrintDistr(const string& name_file)
{
    double bin;
    fstream outfile(name_file, fstream::out| ios::app );
    for(int ik = 0; ik < num_points;ik++)
    {

        bin = float(ik) * step;
        outfile<<bin<<" "<<dnkup[ik]<<"\n";      
    }
    outfile.close();
}


MomentDistr::~MomentDistr(){
    delete [] dnkup;
    delete [] dnkupa;

    for(int i = 0; i < num_points; i++)
        delete [] dnkuplocal[i];
    delete [] dnkuplocal;
}

