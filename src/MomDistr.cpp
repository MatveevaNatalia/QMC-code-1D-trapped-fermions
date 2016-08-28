#include "MomDistr.h"

MomentDistr::MomentDistr(int num_points){
    stored_num_points = num_points;
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
    for(int i = 0; i < stored_num_points; i++)
        dnkupa[i] = dnkuplocal[i][ipop];
}

void MomentDistr::WalkerMatch(int jpop){
    for(int i = 0; i < stored_num_points; i++)
    {
        dnkuplocal[i][jpop] = dnkupa[i];
    }
}

void MomentDistr::WalkerCollect(int nsons){
    for(int i = 0; i < stored_num_points; i++)
        dnkup[i] = dnkup[i] + dnkupa[i]*nsons;
}

void MomentDistr::Normalization(int np, int nkuppt){
    for(int ik = 0; ik < stored_num_points;ik++)
        dnkup[ik] = dnkup[ik]*np/(float(nkuppt));
}


MomentDistr::~MomentDistr(){
    delete [] dnkup;
    delete [] dnkupa;

    for(int i = 0; i < stored_num_points; i++)
        delete [] dnkuplocal[i];
    delete [] dnkuplocal;
}

