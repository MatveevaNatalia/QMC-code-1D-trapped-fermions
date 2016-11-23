#ifndef QMC_H
#define QMC_H

#include <string>
#include <map>
#include <iostream>
using namespace std;

#define BIG  1e30
const int dmnpop = 300;
const double pi = 3.141592653589793;

struct ParamModel{
    int num_comp, num_part, num_walk, nwalk_mean;
    long seed;
    double scat_lenght, scat_lenght_bos, width, alfa, Lmax;

    ParamModel(const map<string, double> & parametersMap){
        //num_comp = parametersMap["ncomp"]; can't be used with
        //'const' key word in declaration of the constructor, because []
        //is marked as non-const, and cannot be used in a const function

        //num_comp = parametersMap.find("ncomp")->second; //2003 standard
        num_comp = parametersMap.at("ncomp"); // 2011 standard
        num_part = parametersMap.at("np");
        width    = parametersMap.at("width");
        scat_lenght_bos = parametersMap.at("aB");
        scat_lenght = parametersMap.at("a");
        alfa = parametersMap.at("alfa");
        num_walk = parametersMap.at("num_walk");
        Lmax = parametersMap.at("num_walk");
        seed = parametersMap.at("seed");
        nwalk_mean = parametersMap.at("nwalk_mean");
    }

};


struct CorFunParam{
    double step, max_value;
    int num_points;
    CorFunParam(double _max_value, int _num_points){
        max_value = _max_value;
        num_points = _num_points;
        step = max_value / num_points;
    }

};

void run (const string & inFile, const string & startingConfig, const string & outDir );


#endif // QMC_H
