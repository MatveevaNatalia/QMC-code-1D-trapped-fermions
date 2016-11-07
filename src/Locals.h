#ifndef LOCALS_H
#define LOCALS_H

#include "qmc.h"


#include <iostream>
#include <fstream>
#include <iomanip>


using namespace std;

/*class Configuration {
    vector<vector<double>> components;
public:
    double GetParticleComp(int npart, int nComp) {
        return components[npart][nComp];
    }
};*/

class Locals{
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

};

#endif // LOCALS_H
