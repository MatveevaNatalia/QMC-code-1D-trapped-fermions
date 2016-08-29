using namespace std;
void PairDistribution_calc(double **, double *, double *, double *, int, double, int, int);
void CalculateFirst(double **xMT, double *graMT_11, int mgr, double dr, int ncomp, int np);
void CalculateSecond(double **xMT, double *graMT_22, int mgr, double dr, int ncomp, int np);
void CalculateCross(double **xMT, double *graMT_12, int mgr, double dr, int ncomp, int np);

