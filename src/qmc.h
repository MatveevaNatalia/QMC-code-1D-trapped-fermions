#ifndef QMC_H
#define QMC_H

#include <string>

using namespace std;

class Energy
{
public:
    double tot, pot, kin;
    void SetZero();
};


void run (const string & inFile, const string & startingConfig, const string & outDir );


#endif // QMC_H
