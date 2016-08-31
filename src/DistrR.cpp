#include "DistrR.h"


using namespace std;

DistributionR::DistributionR(const CorFunParam& pair_distr){
    num_points = pair_distr.num_points;
    step = pair_distr.step;
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
    for(int i = 1; i < (num_points+1); i++)
        dra[i] = draMT[i];
}

void DistributionR::NotAccept(int ipop, int in){
    for(int i = 1; i < (num_points+1); i++)
        dra[i] = drlocal[i][ipop][in];
}

void DistributionR::WalkerMatch(int jpop, int io){
    for(int i = 1; i < (num_points+1); i++)
        drlocal[i][jpop][io] = dra[i];
}

void DistributionR::WalkerCollect(int nsons){
    for(int i = 1; i < (num_points+1); i++)
        dr[i] = dr[i] + nsons * dra[i];
}

void DistributionR::NormalizationGR(int ngr, int ncomp, int np, double step){
    double r1;
    for (int ir = 1; ir < (num_points+1); ir++)
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
    for (int ir = 1; ir < (num_points+1); ir++)
    {
        r1 = float(ir) * step;
        if (ngr > 0)
            dr[ir] = dr[ir]/float(ngr);
        else
            dr[ir] = 0.0;

        dr[ir] = dr[ir]/(2*step);
    }
}

void DistributionR::PairDistrFirst(double **xMT, const ParamModel& param_model )
{
    double Lmax, deltax;
    int bin_number;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part ; i++)
    {
        for(int j = i+1; j < param_model.num_part; j++)
        {
            deltax = fabs(xMT[0][i] - xMT[0][j]);

            if(deltax < Lmax)
            {
                bin_number = deltax/step;
                draMT[bin_number] += 1;
            }

        }
    }
}

void DistributionR::PairDistrSecond(double **xMT, const ParamModel& param_model )
{
    double Lmax, deltax;
    int bin_number;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part; i++)
    {
        for(int j = i+1; j < param_model.num_part; j++)
        {
            deltax = fabs(xMT[1][i] - xMT[1][j]);

            if(deltax < Lmax)
            {
                bin_number = deltax/step;
                draMT[bin_number] += 1;
            }
        }
    }
}

void DistributionR::PairDistrCross(double **xMT, const ParamModel& param_model)
{
    double Lmax, deltax;
    int  bin_number;
    Lmax = step*num_points;

    for(int icomp = 1; icomp < param_model.num_comp; icomp++)
    {
        for(int i = 0; i < param_model.num_part; i++)
        {
            for(int j = 0; j < param_model.num_part; j++)
            {
                deltax = fabs(xMT[0][i] - xMT[icomp][j]);
                if(deltax < Lmax)
                {
                    bin_number = deltax/step;
                    draMT[bin_number] += 1;
                }
            }
        }
    }
}


void DistributionR::DensityFirst(double **xMT, const ParamModel& param_model)
{
    double Lmax;
    int bin_num;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part; i++)
    {
        if(fabs(xMT[0][i]) < Lmax)
        {
            bin_num = fabs(xMT[0][i])/step;
            draMT[bin_num] += 1;
        }
    }
}

void DistributionR::DensitySecond(double **xMT, const ParamModel& param_model)
{
    double Lmax;
    int bin_num;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part; i++)
    {
        if(fabs(xMT[1][i]) < Lmax)
        {
            bin_num = fabs(xMT[1][i])/step;
            draMT[bin_num] += 1;
        }
    }
}

void DistributionR::PrintDistr(const string& name_file)
{
    double coord_bin;
    fstream outfile(name_file, fstream::out| ios::app );
    for(int i = 1; i < (num_points+1);i++)
    {
        coord_bin = float(i) * step + step/2.0;
        outfile<<coord_bin<<" "<<dr[i]<<"\n";
        // cout<<r1<<" "<<gr_11[i]<<" "<<gr_12[i]<<"\n"; Impossible to do in the method?
    }
    outfile.close();
}




DistributionR::~DistributionR(){
    delete[] dra;
    delete[] dr;
    delete[] draMT;
    for(int i = 0; i < (num_points+1); i++)
    {
        for (int j = 0; j < dmnpop; j++) delete [] drlocal[i][j];
        delete [] drlocal[i];
    }
    delete [] drlocal;

}
