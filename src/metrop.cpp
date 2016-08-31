#include "metrop.h"

void MetropolisDif(int ipop, ParamModel param_model, double PsiTotal, double **flocal, Locals& coordinates, Locals& force, int& accepta, int& nprova, double& fvella, int ntemps, int in, int i_VMC)
{
    double fdif, QQ, DDF1, DDF2, DDS, dte = param_model.alfa/4;
    //long seed;
    //seed = param_model.seed;
    cout << "kkk= " << param_model.seed << "\n";
    if(ntemps == 1)
    {
        fvella = 0.0;
        QQ = 0.0;
    }
    else
    {
        fvella = flocal[ipop][in];
        QQ = 0.0;
        for(int inc = 0; inc < param_model.num_comp; inc++)
        {
            for(int ip = 0; ip < param_model.num_part; ip ++)
            {
                DDF1 = force.total[inc][ip][ipop][in] + force.metrop[inc][ip];
                DDF2 = force.total[inc][ip][ipop][in] - force.metrop[inc][ip];
                DDS = coordinates.metrop[inc][ip] - coordinates.total[inc][ip][ipop][in];
                QQ = QQ + 0.5*DDF1*(0.5*dte*DDF2-DDS);
            }
        }
    }

    if(i_VMC == 1)
        fdif = 2.0 * (PsiTotal - fvella); //For_VMC_calculations
    if(i_VMC == 0)
        fdif = 2.0 * (PsiTotal - fvella) + QQ; //SVMC and DMC, QQ because of drift jump

    accepta = 1;
    if(ntemps > 1)
    {
        nprova = nprova + 1;
        if (fdif >= 0)
            accepta = 1;
        else
        {
            accepta = 0;
            if(ran2(&param_model.seed)< exp(fdif))
            {
                cout << "kkk= " << param_model.seed <<" "<<ran2(&param_model.seed) << "\n";
                accepta = 1;
            }
        }
    }
}
