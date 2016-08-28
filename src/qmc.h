#ifndef QMC_H
#define QMC_H

#include <string>

using namespace std;

#define BIG  1e30
const int dmnpop = 300;
const double pi = 3.141592653589793;

struct ParamModel{
    int num_comp, num_part, num_walk;
    long seed;
    double scat_lenght, scat_lenght_bos, width, alfa;
};

struct CorFunParam{
    double step, max_value, num_points;
};


class Energy
{
public:
    double tot, pot, kin;
    void SetZero();
};


void run (const string & inFile, const string & startingConfig, const string & outDir );


#endif // QMC_H
