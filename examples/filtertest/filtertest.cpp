#include <iostream>
#include <cstdlib>
#include <random>
#include <fstream>
#include <typeinfo>
#include "/home/vd/rls/source/RLS_Estimation_Object.h"
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>

using namespace RLS;
using namespace std;

void testFF(const double* Y, double ff, int len ,double init_covar)

{  	
	typedef Matrix< float, Dynamic, 1 > Vec;
    Vec RegF=Vec::Identity(2,1);
   	RLS_Estimator<double> Test_Gen(2, ff ,init_covar);  // (Number of factors ,Forgeting factor, Initial Covarience function)

	//Print output of parameters from the generalized class in .txt file
    ofstream out;
    out.open ("build/FF.txt");

	for (int i = 0; i < len; i++) {
		RegF(0)=pow(i,0);
		RegF(1)=pow(i,1);
        Test_Gen.update_par(RegF,Y[i]); // Update parameters in respect to Input and Output
		out << Test_Gen.getEstimatedOutput(RegF)<<'\t' << (Test_Gen.getEstimatedParameters()(0)) << "\t" << (Test_Gen.getEstimatedParameters()(1)) << endl;	
		}
   		out.close();
    }

void testWindow(const double* Y, int len,double ff, int win ,double init_covar)
{
    typedef Matrix< float, Dynamic, 1 > Vec;
    Vec RegW=Vec::Identity(2,1);
    BlockRLS<double> Test_Block (2, ff , win, init_covar);

	//Check output of parameters from the generalized class
   	ofstream out;
	out.open ("build/RecWin.txt");

    for (int i = 0; i < len; i++) { //for each sample i
		//update polinomial factors
		//cout<<"i="<<i<<endl;
		RegW(0)=pow(i,0);
		RegW(1)=pow(i,1);
		Test_Block.update_par(RegW,Y[i]); 
		out <<Test_Block.getEstimatedOutput(RegW)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	
		//cout << Test_Block.getEstimatedOutput(RegW)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	
		//cout<<"\nRegW=\n"<<RegW<<endl;
	}
   	out.close();
}
   

int main() {

	system("rm build/FF.txt");
	system("rm build/Signal.txt");
	system("rm build/RecWin.txt");

	//Initializing parameters//
	int len = 1000; //Length of Array of Output/Input-Iteration function
	double init_covar = 80;
	double ff = 0.97; 
	int window_size =80;

	//Initializing random signal
	double Y[len]; 
	ofstream SignalFile;
	SignalFile.open ("build/Signal.txt");
	for (int i = 0; i < 2*len/8; i++) {
		Y[i]=2.0 + 0.1*(2.0*rand()/RAND_MAX-1.0);
		SignalFile<<Y[i]<<'\n';
	}
	for (int i = 2*len/8; i <3*len/8; i++) {

		Y[i]=0.70 + 0.1*(2.0*rand()/RAND_MAX-1.0);
		SignalFile<<Y[i]<<'\n';
	}
	for (int i = 3*len/8; i < len; i++) {

		Y[i]=1.3 + 0.1*(2.0*rand()/RAND_MAX-1.0);
		SignalFile<<Y[i]<<'\n';
	}
	SignalFile.close();
	
	//Call Estimating Functions
	testWindow( Y, len,ff, window_size , init_covar);
	testFF( Y,ff, len, init_covar);
return 0;

}

