#ifndef CONFIGDISTR_H
#define CONFIGDISTR_H

#include <vector>
using namespace std;

class ConfigDistr{
private:
    vector<double> distribution;
public:
    ConfigDistr(int num_points):distribution(num_points,0.){}
    void setzero(){
        fill(distribution.begin(), distribution.end(), 0.);
    }

    double & operator[](int i){
        return distribution[i];
    }

    void SetValue(int i, double value){
        distribution[i] = value;
    }
};


#endif // CONFIGDISTR_H
