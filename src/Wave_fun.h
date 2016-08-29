#include "qmc.h"
#include "Locals.h"

using namespace std;
double WaveFunction(const ParamModel&, Locals&);
double WaveFunction_MC(int, int, double, double, double, double **, double *, int, double, int); 
