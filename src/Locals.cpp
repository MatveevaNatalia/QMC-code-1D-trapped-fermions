#include "Locals.h"
#include "Algorithm.h"

using namespace std;

Locals::Locals(ParamModel param_model){

    num_comp_saved = param_model.num_comp;
    num_part_saved = param_model.num_part;
    num_walk_saved = param_model.num_walk;
    seed_saved = param_model.seed;
    alfa_saved = param_model.alfa;
    step_jump = alfa_saved/4.0;

    auxil = new double*[num_comp_saved];
    for(int i = 0; i < num_comp_saved; i++) auxil[i] = new double[num_part_saved];

    total = new double***[num_comp_saved];
    for(int i = 0; i < num_comp_saved; i++)
    {
        total[i] = new double**[num_part_saved];
        for(int j = 0; j < num_part_saved; j++)
        {
            total[i][j] = new double*[dmnpop];
            for(int m = 0; m < dmnpop; m++)
                total[i][j][m] = new double[2];
        }
    }

    metrop = new double*[num_comp_saved];
    for(int i = 0; i < num_comp_saved; i++) metrop[i] = new double[num_part_saved];

}
void Locals::Locals::ReadInitial(const string & startingConfig){
    double etemp;
    fstream filestr(startingConfig, fstream::in | fstream::out);
    filestr>>setprecision(18);
    filestr >> num_walk_saved;
    for (int ipop = 0; ipop < num_walk_saved; ipop++ )
    {
        filestr >> etemp;
        for(int ic = 0; ic < num_comp_saved;ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++)
            {
                filestr >> total[ic][ip][ipop][0];
            }
        }

    }
    filestr.close();
}

void Locals::GenerateInitial()
{
    long seed = -13156;
    for (int ipop = 0; ipop < num_walk_saved; ipop++ )
    {
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++) {total[ic][ip][ipop][0] = ran2(&seed);}

        }
    }

    fstream outfile("../1D_Qt/Data/inprev.dat", fstream::out );
    outfile << setprecision(17);
    outfile<<num_walk_saved<<"\n";
    for(int ipop = 0; ipop < num_walk_saved;ipop++)
    {
        outfile<<ipop<<"\n";
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++)
                outfile<<total[ic][ip][ipop][0]<<"\n";
        }
    }
    outfile.close();
}

void Locals::GaussianJump(int ntemps, int in, int i_VMC, int ipop, double**** force_total){
    double xgaus;

    for (int ic = 0; ic < num_comp_saved; ic++)
    {
        for(int ip = 0; ip < num_part_saved; ip++)
        {
            if (ntemps > 1)
                Gauss1D(&xgaus, alfa_saved, &seed_saved);
            else
                xgaus = 0.0;

            if (i_VMC == 1){
                metrop[ic][ip] = total[ic][ip][ipop][in] + xgaus;

            }
            else
                metrop[ic][ip] = total[ic][ip][ipop][in] + xgaus + force_total[ic][ip][ipop][in] * step_jump;
        }
    }
}

void Locals::Accept(){
    for(int ic = 0; ic < num_comp_saved; ic++)
    {
        for(int ip = 0; ip < num_part_saved; ip++)
        {
            auxil[ic][ip] = metrop[ic][ip];
        }
    }
}

void Locals::NotAccept(int ipop, int in){
    for(int ic = 0; ic < num_comp_saved; ic++)
    {
        for(int ip = 0; ip < num_part_saved; ip++)
            auxil[ic][ip] = total[ic][ip][ipop][in];
    }
}

void Locals::WalkerMatch(int jpop, int io){
    for(int ic = 0; ic < num_comp_saved; ic++)
    {
        for(int ip = 0; ip < num_part_saved; ip++)
            total[ic][ip][jpop][io] = auxil[ic][ip];
    }
}

void Locals::SetZeroForceTotal(int in){
    for(int j = 0; j < num_walk_saved; j++ )
    {
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++)
            {
                total[ic][ip][j][in] = 0.0;
            }
        }
    }
}

void Locals::SetZeroForce(){
    for(int ic = 0; ic < num_comp_saved; ic++)
    {
        for(int ip = 0; ip < num_part_saved; ip++)
        {
            metrop[ic][ip] = 0.0;
        }
    }

}

/*void Locals::PrintConfig(const string& name_file, double** elocal_tot, int in)
{
    fstream outfile(name_file, fstream::out ); //Coordinate::Print
    outfile<<num_walk_saved<<"\n";
    outfile<<setprecision(18);
    for(int i = 0; i < num_walk_saved;i++)
    {
        outfile<<elocal_tot[i][in]<<"\n";
        for(int ic = 0; ic < num_comp_saved; ic++)
        {
            for(int ip = 0; ip < num_part_saved; ip++ )
                outfile<<total[ic][ip][i][in]<<"\n";
        }
    }
    outfile.close();
}*/





Locals::~Locals(){
    for(int i = 0; i < num_comp_saved; i++)
    {
        for(int j = 0; j < num_part_saved; j++)
        {
            for(int m = 0; m < dmnpop; m++)
                delete [] total[i][j][m];

            delete[] total[i][j];
        }
        delete [] total[i];
    }
    delete [] total;

    for(int i = 0; i < num_comp_saved; i++)
        delete [] auxil[i];
    delete [] auxil;


    for(int i = 0; i < num_comp_saved; i++)
        delete [] metrop[i];
    delete [] metrop;
}


