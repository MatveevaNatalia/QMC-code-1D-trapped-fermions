#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>
#include "qmc.h"

void showhelpinfo(char *s);
using namespace std;

int main (int argc, char *argv[])
{
    std::string  inFile = "";
    std::string startingConfig = "";
    std::string outDir  = "";


    char tmp;
    /*if the program is ran witout options ,it will show the usgage and exit*/
    if(argc == 1)
    {
        showhelpinfo(argv[0]);
        exit(1);
    }

    while((tmp=getopt(argc,argv,"i:o:c:"))!=-1)
    {
        if(tmp == 'i') inFile = optarg;
            else if(tmp == 'c') startingConfig = optarg;
            else if(tmp == 'o') outDir = optarg;
    }

    if(inFile == "") {
        std::cout << "Error: input parameters file is required!\n";
        showhelpinfo(argv[0]);
        exit(1);
    }


    if(outDir == "") {
        std::cout << "Error: output folder name is required!\n";
        showhelpinfo(argv[0]);
        exit(1);
    }

    std::cout << "Starting the simulation with: \n";
    std::cout << " - input parameters file \"" << inFile << "\"\n";

    std::cout << " - output folder \"" << outDir << "\"\n";
    if(startingConfig == "")
        std::cout << " - no starting configuration was provided, generating new one\n";
    else
        std::cout << " - input starting configuration \"" << startingConfig << "\"\n";

    run (inFile, startingConfig, outDir);

    return 0;
}


/*funcion that show the help information*/
void showhelpinfo(char *s)
{
  cout<<"Usage:   "<<s<<" [-option] [argument]"<<endl;
  cout<<"option:  "<<"-i input file"<<endl;
  cout<<"         "<<"-o output directory"<<endl;
  cout<<"         "<<"-c starting configuration"<<endl;
}
