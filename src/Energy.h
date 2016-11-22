#ifndef ENERGY_H
#define ENERGY_H
#include "qmc.h"
#include "Locals.h"
using namespace std;

//void Energy_calc(double **xMT, double** FMT, Energy& energy, const ParamModel& param_mode);
void Energy_calc(const Configuration& xMT, Configuration& FMT, Energy& energy, const ParamModel& param_mode);
double EnergyPartial(double xk, double xi, double scat_length);
double ForcePartial(double xk, double xi, double scat_length);

#endif
