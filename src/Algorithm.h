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
void MetropolisDif(int ipop, ParamModel& param_model, WaveFunction& wave_func, Locals& coordinates, Locals& force, int& accepta, int& nprova, int ntemps, int in, int i_VMC);
int BranchingCalc(ParamModel& param_model, Energy& energy, int accepta, int ntemps,int nacc, int nprova);


#endif
