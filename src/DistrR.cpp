#include "DistrR.h"
using namespace std;

DistributionR::DistributionR(const CorFunParam& param):
    average(param.num_points),
    auxil(param.num_points),   
    oldPage(dmnpop, ConfigDistr(param.num_points))
{
    num_points = param.num_points;
    step = param.step;
}

void DistributionR::SetZero(){
    average.setzero();
}

void DistributionR::SetZeroAx(){
    auxil.setzero();
}

void DistributionR::NotAccept(int ipop){
    auxil = oldPage[ipop];
}

void DistributionR::WalkerMatch(){
    newPage.push_back(auxil);
}

void DistributionR::WalkerCollect(int nsons){
    for(int i = 1; i < (num_points+1); i++)
        average[i] += nsons * auxil[i];
}

void DistributionR::NormalizationGR(int nsons_total, const ParamModel& param){
    double norm;
    norm = param.num_comp * param.num_part * step;
    for (int i = 1; i < (num_points+1); i++)
    {
        if (nsons_total > 0)
            average[i] = average[i]/float(nsons_total);
        else
            average[i] = 0.;
        average[i] = average[i]/norm;
    }
}

void DistributionR::NormalizationNR(int nsons_total){
    for (int i = 1; i < (num_points+1); i++)
    {
        if (nsons_total > 0)
            average[i] = average[i]/float(nsons_total);
        else
            average[i] = 0.;

        average[i] = average[i]/(2*step);
    }
}

void DistributionR::PairDistrFirst(const Configuration& xMT, const ParamModel& param_model )
{
    double Lmax, deltax;
    int bin_number;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part ; i++)
    {
        for(int j = i+1; j < param_model.num_part; j++)
        {
            deltax = fabs(xMT.GetParticleComp(0,i) - xMT.GetParticleComp(0,j));
            if(deltax < Lmax)
            {
                bin_number = deltax/step;               
                auxil[bin_number] += 1;
            }

        }
    }
}

void DistributionR::PairDistrCross(const Configuration& xMT, const ParamModel& param_model)
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
                deltax = fabs(xMT.GetParticleComp(0,i) - xMT.GetParticleComp(icomp,j));
                if(deltax < Lmax)
                {
                    bin_number = deltax/step;                   
                    auxil[bin_number] += 1;
                }
            }
        }
    }
}

void DistributionR::DensityFirst(const Configuration& xMT, const ParamModel& param_model)
{
    double Lmax, x_modul;
    int bin_num;
    Lmax = step*num_points;

    for(int i = 0; i < param_model.num_part; i++)
    {
        x_modul = fabs(xMT.GetParticleComp(0,i));
        if(x_modul < Lmax)
        {
            bin_num = x_modul/step;           
            auxil[bin_num] += 1;
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
        outfile<<coord_bin<<" "<<average[i]<<"\n";
    }
    outfile.close();
}

void DistributionR::PageSwap(){
     oldPage = newPage;
     newPage.clear();
}

