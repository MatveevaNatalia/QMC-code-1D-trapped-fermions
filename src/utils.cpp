#include "utils.h"
#include <iomanip>
#include <fstream>
using namespace std;

map<string, double> FillMap(const string & inFile, const string & outDir){
    map<string, double> paramMap;
    fstream filestr(inFile, fstream::in | fstream::out);
    fstream outfile(outDir+"/Output.dat", fstream::out| ios::app );
    outfile<<setprecision(18);
    while(!filestr.eof()){
        string str;
        double value;
        filestr >> str >> value;
        paramMap[str] = value;
        outfile<<str <<" "<<paramMap[str]<<endl;
    }
    outfile.close();

    return paramMap;
}
