#include "Statistics.h"
#include <iostream>
#include <fstream>
#include <iomanip> 
#include <cmath>
using namespace std;

void Extrapolation(string name_VMC, string name_DMC, string name_extr, string name_extr1, int npt)
{

    double *rr, *fun_VMC, *fun_DMC, *err_VMC, *err_DMC, sigma;
    rr = new double[npt];
    fun_VMC = new double[npt];
    fun_DMC = new double[npt];
    err_VMC = new double[npt];
    err_DMC = new double[npt];


    fstream filestr(name_VMC, fstream::in | fstream::out);
    fstream filestr1(name_DMC, fstream::in | fstream::out);

    for (int i = 0; i< npt; i++ )
    {
        filestr >>setprecision(12)>> rr[i]>>fun_VMC[i]>>err_VMC[i];
        //        cout<<i<<" "<<rr[i]<<" "<<fun_VMC[i]<<" "<<err_VMC[i]<<"\n";
        filestr1 >>setprecision(12)>> rr[i]>>fun_DMC[i]>>err_DMC[i];
    }
    filestr.close();
    filestr1.close();

    fstream outfile(name_extr, fstream::out );
    fstream outfile1(name_extr1, fstream::out );

    for (int i = 0; i< npt; i++ )
    {
        sigma = sqrt(err_VMC[i]*err_VMC[i]+err_DMC[i]*err_DMC[i]);
        outfile<<setprecision(12)<< rr[i]<<" "<<fun_DMC[i]*fun_DMC[i]/fun_VMC[i]<<" "<<sigma<<"\n";
        outfile1<<setprecision(12)<< rr[i]<<" "<<2.*fun_DMC[i]-fun_VMC[i]<<" "<<sigma<<"\n";
    }

    outfile.close();
    outfile1.close();

    delete [] rr;
    delete [] fun_DMC;
    delete [] fun_VMC;
    delete [] err_VMC;
    delete [] err_DMC;

}

void Normalization(string name_file, string name_file_1, int num, int np, int nc)
{

    double * kk, *cor_fun, *cor_fun_err;
    double temp, sum, dk, norm;

    kk = new double[num];
    cor_fun = new double[num];
    cor_fun_err = new double[num];

    fstream filestr(name_file, fstream::in | fstream::out);
    for (int i = 0; i< num; i++ )
    {
        filestr >>setprecision(12)>> kk[i]>>cor_fun[i]>>cor_fun_err[i];
        //	cout<<kk[i]<<" "<<cor_fun[i]<<" "<<cor_fun_err[i]<<"\n";
    }
    filestr.close();


    dk = kk[1]-kk[0];
    sum = 0;
    for(int i = 0; i < num; i++)
    {
        if(cor_fun[i]>0)
            sum = sum + dk * cor_fun[i];
        //        cout<<i<<" "<<sum<<"\n";
    }

    norm = (np*nc)/(2.0*sum);//Here normalization on the number of particles per component

    fstream outfile(name_file_1, fstream::out );
    for(int i = 0; i < num; i++)
    {
        outfile<<kk[i]<<" "<<cor_fun[i]*norm<<" "<<cor_fun_err[i]*norm<<"\n";
    }
    outfile.close();
    //cout<<" Normalization for "<<name_file<<" = "<<norm<<"sum= "<<sum<<"\n";

    delete [] kk;
    delete [] cor_fun;
    delete [] cor_fun_err;
}

void Stat_Energy(string name_file, string name_file_1, string name_file_2, int ncol, int nel, int n_notused,  int nfull)
{

    int ntotal, nubloks = 0, ntemp;
    string str;
    double *etemp, *sum, *emig1, *emig2, *eblokE;

    etemp = new double[ncol];
    sum = new double[ncol];
    emig1 = new double[ncol];
    emig2 = new double[ncol];
    eblokE = new double[ncol];

    int icount = 0;

    fstream filestr(name_file, fstream::in | fstream::out);
    for(int j = 0; j < ncol; j++) {sum[j] = 0.0;}
    while(true)
    {
        bool stop = false;
        if(!(filestr >> ntemp)) {stop = true; break;} // to break here
        //    cout<<ntemp<<endl;
        for(int j = 0; j < ncol; j++) {filestr >> etemp[j];}
        icount = icount +1;
        if(icount >= n_notused)
        {
            for(int j = 0; j < ncol; j++) {sum[j] = sum[j] + etemp[j];}
        }
    }
    filestr.close();

    ntotal = icount;
    //cout<<sum[0]<<" "<<sum[1]<<" "<<icount<<endl;

    fstream outfile(name_file_1, fstream::out );

    fstream filestr1(name_file, fstream::in | fstream::out);

    for(int m = 0; m < ncol; m++) {emig1[m] = 0.0; emig2[m] = 0.0;}

    for(int i = 0; i < n_notused; i++)
    {
        filestr1 >> ntemp;
        for(int m = 0; m < ncol; m++) {filestr1 >> etemp[m];}
    }
    //cout <<nubloks<<" "<<nel<<"\n";

    //for(int i = 0; i < nubloks; i++)
    while(true)
    {
        bool stop = false;

        for(int m = 0; m < ncol; m++) {eblokE[m]  = 0.0;}

        ////for(int j = 0; j < nel; j++)
        {
            //	filestr1 >> ntemp;
            if(!(filestr1 >> ntemp)) { stop = true; break;}
            for(int m = 0; m < ncol; m++) {filestr1 >> etemp[m]; }
            for(int m = 0; m < ncol; m++){eblokE[m] = eblokE[m] + etemp[m]; }
        }

        if(stop) break;
        nubloks++;

        for(int m = 0; m < ncol; m++)
        {
            eblokE[m] = eblokE[m]/double(nel);
            emig1[m] = emig1[m] + eblokE[m];
            emig2[m] = emig2[m] + eblokE[m]*eblokE[m];
            //cout<<emig1[1]<<" "<<emig2[1]<<endl;
        }
    }
    filestr1.close();

    //cout<<nubloks<<endl;

    for(int m = 0; m < ncol; m++)
    {
        emig1[m] = emig1[m]/double(nubloks);
        emig2[m] = emig2[m]/double(nubloks);
        //cout<<emig1[0]<<" "<<emig2[0]<<endl;
        emig2[m] = sqrt((emig2[m] - emig1[m]*emig1[m])/double(nubloks));

        //cout<<emig1[0]*emig1[0]<<" "<<emig2[0]<<" "<<double(nubloks)<<endl;
    }


    for(int m = 0; m < ncol; m++)
    {
        outfile <<setprecision(18)<< sum[m]/double(ntotal - n_notused)<<" "<<emig2[m]<<" ";
        //cout<<sum[0]/double(ntotal - n_notused)<<" "<<emig2[0]<<"\n";

    }
    outfile.close();
    //-------------------------------------------------------------------//
    if(nfull == 1)
    {
        //fstream outfile1(name_file_2, fstream::out| ios::app );
        fstream outfile1(name_file_2, fstream::out );

        for(nel = 1; nel < int(double(ntotal - n_notused)/8); nel++)
        {
            nubloks = int(double(ntotal - n_notused)/nel);


            fstream filestr2(name_file, fstream::in | fstream::out);

            for(int m = 0; m < ncol; m++) {emig1[m] = 0.0; emig2[m] = 0.0; eblokE[m] = 0.0; etemp[m] = 0.0;}

            for(int i = 0; i < n_notused; i++)
            {
                filestr2 >> ntemp;
                for(int m = 0; m < ncol; m++) {filestr2 >> etemp[m];}
            }


            for(int i = 0; i < nubloks; i++)
            {
                for(int m = 0; m < ncol; m++) {eblokE[m]  = 0.0;}

                for(int j = 0; j < nel; j++)
                {
                    filestr2 >> ntemp;
                    for(int m = 0; m < ncol; m++) {filestr2 >> etemp[m];}
                    for(int m = 0; m < ncol; m++){eblokE[m] = eblokE[m] + etemp[m];}
                }

                for(int m = 0; m < ncol; m++)
                {
                    eblokE[m] = eblokE[m]/double(nel);
                    emig1[m] = emig1[m] + eblokE[m];
                    emig2[m] = emig2[m] + eblokE[m]*eblokE[m];
                }
            }
            filestr2.close();

            for(int m = 0; m < ncol; m++)
            {
                emig1[m] = emig1[m]/double(nubloks);
                emig2[m] = emig2[m]/double(nubloks);
                emig2[m] = sqrt((emig2[m] - emig1[m]*emig1[m])/double(nubloks));
            }

            //		cout<<nel<<" "<<nubloks<<" ";
            outfile1<<nel<<" "<<nubloks<<" ";

            for(int m = 0; m < ncol; m++)
            {
                outfile1<<" "<<emig2[m]<<" ";
                //			cout<<" "<<emig2[m]<<" ";

            }
            //		cout<<"\n";
            outfile1<<"\n";

        }
        outfile1.close();
    }

    delete [] sum;
    delete [] etemp;
    delete [] emig1;
    delete [] emig2;
    delete [] eblokE;
}


//_______________________________________________________________________________________________//

void statistics_cor(string name_file, string name_file_1, int npt, int ncomp, int nfora, int nel)
{

    int ntemp;

    int nublocks = 0;
    //	nublocks = (nblck - nfora)/nel;

    double **func, **func2, **func_block, *rr, *func_block_temp;
    double rtemp, ftemp;

    func = new double*[ncomp];
    for(int i = 0; i < ncomp; i++) func[i] = new double[npt];

    func_block_temp = new double [ncomp];

    func2 = new double*[ncomp];
    for(int i = 0; i < ncomp; i++) func2[i] = new double[npt];

    func_block = new double*[ncomp];
    for(int i = 0; i < ncomp; i++) func_block[i] = new double[npt];

    rr = new double[npt];


    for(int i = 0; i < ncomp; i++)
    {
        for(int j = 0; j < npt; j++)
        {
            func[i][j] = 0.0;
            func2[i][j] = 0.0;
        }
    }

    fstream filestr(name_file, fstream::in | fstream::out);


    for(int i = 0; i < nfora; i++)
    {
        for (int j = 0; j < npt; j++)
        {
            filestr >> rtemp;
            for(int icomp = 0; icomp < ncomp; icomp ++)
                filestr >> ftemp;
        }
    }



    //	for(int iblck = 0; iblck < nublocks; iblck++)
    while(true)
    {
        bool stop = false;

        for(int i = 0; i < ncomp; i++)
        {
            for(int j = 0; j < npt; j++)
                func_block[i][j] = 0.0;
        }

        for(int i = 0; i < nel; i++)
        {
            for(int j = 0; j < npt; j++)
            {

                if(!(filestr >> rr[j])){ stop = true; break; }

                for(int icomp = 0; icomp < ncomp; icomp++)
                {
                    filestr >> ftemp;
                    func_block[icomp][j] = func_block[icomp][j]+ftemp;
                }

            }
            if(stop) break;
        }
        if(stop) break;

        nublocks++;

        for(int j = 0; j < npt; j++)
        {
            for(int icomp = 0; icomp < ncomp; icomp++)
            {
                func_block[icomp][j] = func_block[icomp][j]/double(nel);
                //            cout<<func_block[icomp][j]<<" ";
                func[icomp][j] = func[icomp][j] + func_block[icomp][j];
                func2[icomp][j] = func2[icomp][j] + func_block[icomp][j]*func_block[icomp][j];

            }
            //		cout<<endl;
        }
    }

    cout<<"nublocks= "<<nublocks<<endl;

    fstream outfile(name_file_1, fstream::out);

    double sum1, sum2;
    for(int j = 0; j < npt; j++)
    {
        outfile << rr[j]<<" ";
        sum1 = 0.0;
        sum2 = 0.0;
        for(int icomp = 0; icomp < ncomp; icomp++)
        {
            func[icomp][j] = func[icomp][j] / double(nublocks);
            func2[icomp][j] = func2[icomp][j] / double(nublocks);
            func2[icomp][j]=sqrt((func2[icomp][j]-func[icomp][j]*func[icomp][j])/double(nublocks));
            outfile <<func[icomp][j]<<" "<<func2[icomp][j]<<" ";
            sum1 = sum1 + func[icomp][j];
            sum2 = sum2 + func2[icomp][j]*func2[icomp][j];
        }
        //	if(i_single == 0){outfile << sum1/double(ncomp)<<" "<<sqrt(sum2);}
        outfile<<"\n";
    }

    filestr.close();
    outfile.close();


    delete [] rr;

    for(int i = 0; i < ncomp; i++)
        delete[] func[i];
    delete [] func;

    for(int i = 0; i < ncomp; i++)
        delete[] func2[i];
    delete [] func2;

    for(int i = 0; i < ncomp; i++)
        delete[] func_block[i];
    delete [] func_block;

    delete [] func_block_temp;


}

