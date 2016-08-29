#include "qmc.h"

using namespace std;
//void Energy_calc(double **, double **, double *, double *, double *, int, int, double, double, double);

//void Energy_calc(double **xaux, double **FT, Energy& energy,  int ncomp, int np, double aB, double a, double width);

void Energy_calc(double **xMT, double** FMT, Energy& energy, ParamModel param_mode);
double EnergyPartial(double xk, double xi, double scat_length);
double ForcePartial(double xk, double xi, double scat_length);
