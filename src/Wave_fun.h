#include "qmc.h"
#include "Locals.h"

using namespace std;
//double WaveFunction(ParamModel param_model, Coordinates coordinates, double* &PsiTotal);
double WaveFunction(ParamModel&, Locals&);

//double WaveFunction(int, int, double, double, double, double **, double *);
double WaveFunction_MC(int, int, double, double, double, double **, double *, int, double, int); 
