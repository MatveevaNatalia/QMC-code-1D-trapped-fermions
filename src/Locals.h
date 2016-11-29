#ifndef LOCALS_H
#define LOCALS_H

#include "qmc.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

class Configuration {
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
        for (unsigned int i = 0; i < components.size(); i++)
            for (unsigned int j = 0; j < components[i].size(); j++)
                cout << components[i][j] << endl;
    }

    void SetZero(){
       for (unsigned int i = 0; i < components.size(); i++)
           for (unsigned int j = 0; j < components[i].size(); j++)
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

    void ReadInitial(const string & startingConfig);
    void GenerateInitial(const string & startingConfig);
    void GaussianJump(int ntemps, bool i_VMC, int ipop, Locals& force);
    void Accept();
    void NotAccept(int ipop);
    void WalkerMatch();
    void SetZero();
    void PageSwap();
    void PrintOld();
    void PrintNew();
    void PrintAll(const string & startingConfig);

};

void Gauss1D(double * x, double alfa, long *kkk);
float ran2(long *idum);

#endif // LOCALS_H
