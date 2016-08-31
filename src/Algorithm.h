#ifndef ALGORITHM_H
#define ALGORITHM_H

#include "qmc.h"
#include "Locals.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

float ran2(long *);
void Initial(double ****, int, int, int);
void Gauss1D(double *, double, long *);
void MetropolisDif(int, int, int, double, double **, double **, double ****, double ****, double **, int *, int *, double *, int, long *, int, double, int);
//void MetropolisDif(int ipop, ParamModel param_model, double PsiTotal, double **flocal, Locals& coordinates, Locals& force, int& accepta, int& nprova, double& fvella, int ntemps, int in, int i_VMC);

//void MetropolisDif(int, ParamModel, double, double **, Locals&, Locals&, int&, int&, double&, int, int, int);


void BranchingCalc(int *, int, int, int, int, double, int, double, double, double, long *, int);



#endif
