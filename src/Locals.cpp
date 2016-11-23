#include "Locals.h"
#include "Algorithm.h"

using namespace std;

void Locals::ReadInitial(const string & startingConfig){

    fstream filestr(startingConfig, fstream::in | fstream::out);
    filestr>>setprecision(18);
    filestr >> num_walk;
    for (int ipop = 0; ipop < num_walk; ipop++ ){
        double etemp;
        filestr >> etemp;
        vector<vector<double>> vec_comp;
        for(int ic = 0; ic < num_comp;ic++){
            vector<double> vec_part;
            for(int ip = 0; ip < num_part; ip++){
                 double temp_coord;
                 filestr >> temp_coord;
                 vec_part.push_back(temp_coord);
            }
            vec_comp.push_back(vec_part);
        }
        Configuration config(vec_comp);
        PushOldConfig(config);
    }

}

void Locals::GenerateInitial(const string & startingConfig){
    long seed = -13156;

    fstream outfile(startingConfig, fstream::out );
    outfile << setprecision(18);
    outfile<<num_walk<<"\n";
    for(int ipop = 0; ipop < num_walk;ipop++){
        outfile<<ipop<<"\n";
        vector<vector<double>> vec_comp;
        for(int ic = 0; ic < num_comp; ic++){
            vector<double> vec_part;
            for(int ip = 0; ip < num_part; ip++){
                double temp_coord = ran2(&seed);
                outfile<<temp_coord<<"\n";
                vec_part.push_back(temp_coord);
            }
            vec_comp.push_back(vec_part);
        }
        Configuration config (vec_comp);
        PushOldConfig(config);
    }
    outfile.close();
}

void Locals::GaussianJump(int ntemps, int i_VMC, int ipop, Locals& force)
{
    double xgaus, x_old, force_old, x_new;
    vector<vector<double>> vec_comp;
    for (int ic = 0; ic < num_comp; ic++)
    {
        vector<double> vec_part;
        for(int ip = 0; ip < num_part; ip++)
        {

            x_old = oldPage[ipop].GetParticleComp(ic, ip);
            if (ntemps > 1)
                Gauss1D(&xgaus, alfa, &seed);
            else
                xgaus = 0.0;

            if (i_VMC == 1)
            {
                x_new = x_old + xgaus;
                vec_part.push_back(x_new);
            }
            else
            {
                force_old = force.oldPage[ipop].GetParticleComp(ic, ip);
                x_new = x_old + xgaus + force_old * step_jump;
                vec_part.push_back(x_new);
            }
        }
        vec_comp.push_back(vec_part);
    }
    Configuration config(vec_comp);
    metrop = config;
}

void Locals::Accept(){
    auxil = metrop;
}

void Locals::NotAccept(int ipop){
    auxil = oldPage[ipop];
}

void Locals::WalkerMatch(){
    newPage.push_back(auxil);
}

 void Locals::SetZero(){
    for (unsigned int i = 0; i < oldPage.size(); i++)
        oldPage[i].SetZero();
 }

 void Locals::PageSwap(){
     oldPage = newPage;
     newPage.clear();
 }

void Locals::PrintOld(){
    for (unsigned int i = 0; i < oldPage.size(); i++)
        oldPage[i].PrintConfig();
}

void Locals::PrintNew(){
    for (unsigned int i = 0; i < newPage.size(); i++)
        newPage[i].PrintConfig();
}

void Locals::PrintAll(const string & startingConfig, int in){
    fstream outfile(startingConfig, fstream::out );
    outfile<<num_walk<<"\n";
    outfile<<setprecision(18);
    for(int i = 0; i < num_walk;i++)
    {
        for(int ic = 0; ic < num_comp; ic++)
        {
            for(int ip = 0; ip < num_part; ip++ )

                {outfile<<oldPage[i].GetParticleComp(ic,ip)<<"\n"; }
        }
    }
    outfile.close();
}




