#ifndef LOCALS_H
#define LOCALS_H

#include "qmc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

class Configuration { // Class State
    vector<vector<double>> components;
public:

    Configuration(vector<vector<double>> _components):
        components(_components){}

    Configuration(int nComp, int npart):
        components(nComp, vector<double>(npart,0)){}

    double GetParticleComp(int nComp, int npart) const {
        return components[nComp][npart];    
    }

    void SetParticleComp(int nComp, int npart, double value){
        components[nComp][npart] = value;
    }

    void PrintConfig() const {
        for (int i = 0; i < components.size(); i++)
            for (int j = 0; j < components[i].size(); j++)
                cout << components[i][j] << endl;
    }

    void SetZero(){
       for (int i = 0; i < components.size(); i++)
           for (int j = 0; j < components[i].size(); j++)
               components[i][j] = 0.;
    }
};

class Locals{
private:
    int num_comp, num_part, num_walk;
    long seed;
    double alfa, step_jump;

public:
    vector<Configuration> oldPage; // index in
    vector<Configuration> newPage; // index io
    Configuration auxil, metrop;


    Locals(ParamModel& param):
        auxil(param.num_comp, param.num_part),
        metrop(param.num_comp, param.num_part),
        oldPage(param.num_walk, Configuration(param.num_comp,param.num_part)),
        newPage(param.num_walk, Configuration(param.num_comp,param.num_part))
    {
        num_comp = param.num_comp;
        num_part = param.num_part;
        num_walk = param.num_walk;
        seed = param.seed;
        alfa = param.alfa;
        step_jump = alfa/4.0;

    }

    void PushOldConfig(Configuration config){
        oldPage.push_back(config);
    }

    void PushNewConfig(Configuration config){
        newPage.push_back(config);
    }

    //double GetOldCoord( int iWalk, int iComp, int iPart) {
    //    return oldPage[iWalk].GetParticleComp(iComp, iPart);
    //}

    //double GetNewCoord( int iWalk, int iComp, int iPart) {
    //    return newPage[iWalk].GetParticleComp(iComp, iPart);
    //}

    void ReadInitial(const string & startingConfig);
    void GenerateInitial(const string & startingConfig);
    void GaussianJump(int ntemps, int i_VMC, int ipop, Locals& force);
    void Accept();
    void NotAccept(int ipop);
    void WalkerMatch();
    void SetZero();
    void PageSwap();
    void PrintOld();
    void PrintNew();
};

void Gauss1D(double * x, double alfa, long *kkk);
float ran2(long *idum);








/*class Configuration {
    vector<vector<double>> components;
public:
    double GetParticleComp(int npart, int nComp) {
        return components[npart][nComp];
    }
};

class Locals{
private:
    int num_comp, num_part, num_walk;
    long seed;
    double alfa, step_jump;
public:

    Configuration auxil, metrop;
    vector<Configuration> oldPage; // index in
    vector<Configuration> newPage; // index io

    //double ****total, **auxil, **metrop;


    Locals(ParamModel param_model);
    void ReadInitial(const string & startingConfig);

    void GenerateInitial();

   // void GaussianJump(int ntemps, int in, int i_VMC, int ipop, double ****FF);

    GaussianJump(int ntemps, int in, int i_VMC, int ipop, Locals& force);

    void Accept();

    void NotAccept(int ipop, int in);

    void WalkerMatch(int jpop, int io);

    void SetZeroForceTotal(int in);

    void SetZeroForce();

//    void PrintConfig(const string& name_file, double** elocal_tot, int in);


    ~Locals();

};*/










/*class Locals{
private:
    int num_comp_saved, num_part_saved, num_walk_saved;
    long seed_saved;
    double alfa_saved, step_jump;
public:
//    Configuration auxil, metrop;
//    vector<Configuration> oldPage;
//    vector<Configuration> newPage;

    double ****total, **auxil, **metrop;


    Locals(ParamModel param_model);
    void ReadInitial(const string & startingConfig);

    void GenerateInitial();

    void GaussianJump(int ntemps, int in, int i_VMC, int ipop, double ****FF);

    //GaussianJump(int ntemps, int in, int i_VMC, int ipop, Locals& force);
    // { force.oldPage[5].components}
    void Accept();

    void NotAccept(int ipop, int in);

    void WalkerMatch(int jpop, int io);

    void SetZeroForceTotal(int in);

    void SetZeroForce();

//    void PrintConfig(const string& name_file, double** elocal_tot, int in);


    ~Locals();

};*/

#endif // LOCALS_H
