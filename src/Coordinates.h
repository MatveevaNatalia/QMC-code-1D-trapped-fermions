#ifndef COORDINATES_H
#define COORDINATES_H

#include "qmc.h"
#include "Algorithm.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

class Coordinates{
private:
    int num_comp_saved, num_part_saved, num_walk_saved;
    long seed_saved;
    double alfa_saved, step_jump;
public:
    double ****x, **xaux, **xMT;
    Coordinates(ParamModel param_model);
    void ReadInitial(const string & startingConfig);

    void GenerateInitial();

    void GaussianJump(int ntemps, int in, int i_VMC, int ipop, double ****FF);

    void Accept();

    void NotAccept(int ipop, int in);

    void WalkerMatch(int jpop, int io);

   ~Coordinates();

};

#endif // COORDINATES_H
