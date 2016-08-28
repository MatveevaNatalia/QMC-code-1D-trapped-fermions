#include "DistrR.h"


DistributionR::DistributionR(int num_points){
    stored_num_points = num_points;
    dr = new double[num_points+1];
    dra = new double[num_points+1];
    draMT = new double[num_points+1];

    drlocal = new double**[num_points+1];
    for(int i = 0; i < (num_points+1); i++){
        drlocal[i] = new double*[dmnpop];
        for(int j = 0; j < dmnpop; j++ ) drlocal[i][j] =  new double[2];
    }
}

void DistributionR::SetZero(){
    setzero(dr);
}

void DistributionR::SetZeroAx(){
    setzero(dra);
    setzero(draMT);
}

void DistributionR::Accept(){
    for(int i = 1; i < (stored_num_points+1); i++)
        dra[i] = draMT[i];
}

void DistributionR::NotAccept(int ipop, int in){
    for(int i = 1; i < (stored_num_points+1); i++)
        dra[i] = drlocal[i][ipop][in];
}

void DistributionR::WalkerMatch(int jpop, int io){
    for(int i = 1; i < (stored_num_points+1); i++)
        drlocal[i][jpop][io] = dra[i];
}

void DistributionR::WalkerCollect(int nsons){
    for(int i = 1; i < (stored_num_points+1); i++)
        dr[i] = dr[i] + nsons * dra[i];
}

void DistributionR::NormalizationGR(int ngr, int ncomp, int np, double step){
    double r1;
    for (int ir = 1; ir < (stored_num_points+1); ir++)
    {
        r1 = float(ir) * step;
        if (ngr > 0)
            dr[ir] = dr[ir]/float(ngr);
        else
            dr[ir] = 0.0;
        dr[ir] = dr[ir]/(ncomp*np*step);
    }
}

void DistributionR::NormalizationNR(int ngr, float step){
    double r1;
    for (int ir = 1; ir < (stored_num_points+1); ir++)
    {
        r1 = float(ir) * step;
        if (ngr > 0)
            dr[ir] = dr[ir]/float(ngr);
        else
            dr[ir] = 0.0;

        dr[ir] = dr[ir]/(2*step);
    }
}

DistributionR::~DistributionR(){
    delete[] dra;
    delete[] dr;
    delete[] draMT;
    for(int i = 0; i < (stored_num_points+1); i++)
    {
        for (int j = 0; j < dmnpop; j++) delete [] drlocal[i][j];
        delete [] drlocal[i];
    }
    delete [] drlocal;

}
