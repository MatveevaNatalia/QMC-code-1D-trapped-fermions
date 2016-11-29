#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "qmc.h"
#include "Locals.h"
#include "Wave_fun.h"
#include "Energy.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

float ran2(long *);
void Initial(double ****, int, int, int);
void Gauss1D(double *, double, long *);
bool MetropolisDif(int ipop, ParamModel& param_model, WaveFunction& wave_func, Locals& coordinates, Locals& force, int& nprova, int ntemps, bool i_VMC);
int BranchingCalc(ParamModel& param_model, Energy& energy, bool accepta, int ntemps, int nacc, int nprova, int ipop);


#endif
