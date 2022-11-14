#include <iostream>
#include <cstdlib>
#include <random>
#include <fstream>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <Eigen/Dense>
#include <random>

#include "../../source/RLS_Estimation_Object.h"

using namespace RLS;
using namespace std;
using namespace Eigen;

#define real_num float
typedef Matrix <real_num , Dynamic, Dynamic > Mat;
typedef Matrix< real_num, Dynamic, 1 > Vec;

void Forgeting_Factor_Method(const real_num* Y, real_num ff, int len , real_num init_covar)

{  	
	typedef Matrix< real_num, Dynamic, 1 > Vec;
    Vec RegF=Vec::Identity(2,1);
   	RLS_Estimator<real_num> Test_Gen(2, ff ,init_covar);  // (Number of factors ,Forgeting factor, Initial Covarience function)

	//Print output of parameters from the generalized class in .txt file
    ofstream out;
    out.open ("FF.txt");

	for (int i = 0; i < len; i++) {
		RegF(1)=i;
        Test_Gen.update_par(RegF,Y[i]); // Update parameters in respect to Input and Output
		out << Test_Gen.getEstimatedOutput(RegF)<<'\t' << (Test_Gen.getEstimatedParameters()(0)) << "\t" << (Test_Gen.getEstimatedParameters()(1)) << endl;	
		}
   		out.close();
}

void Rec_Window_Method(const real_num* Y, int len,real_num ff, int win ,real_num init_covar)
{
    typedef Matrix< real_num, Dynamic, 1 > Vec;
    Vec RegW=Vec::Identity(2,1);
    BlockRLS<real_num> Test_Block (2, 1 , win, init_covar);

	//Check output of parameters from the generalized class
   	ofstream out;
	out.open ("RecWin.txt");

    for (int i = 0; i < len; i++) { //for each sample i
		//update polinomial factors
		RegW(1)=i;
		Test_Block.update_par(RegW,Y[i],i); 
		out <<Test_Block.getEstimatedOutput(RegW)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	

	}
   	out.close();
}

void Aug_Cholesky_Method(int win, real_num init_covar, Vec Y, int start, int l )
{
	ofstream out;
	out.open ("aug_chol.txt");

	Augmented_Cholesky_RLS_Estimator<real_num> Test_AugChol (2 , win);
	Vec RegW=Vec::Identity(2,1);
    // Update Cholesky matrix and compute Theta again
    for(int j=1; j<l; j++){ 
		RegW(1)=j;
   		Test_AugChol.update_par(RegW, Y[j] );
		out << Test_AugChol.getEstimatedOutput(RegW)<<'\t' << (Test_AugChol.getEstimatedParameters()(0)) << "\t" << (Test_AugChol.getEstimatedParameters()(1)) << endl;	
	}
   	out.close();
}

int main() {
	system("rm FF.txt");
	system("rm Signal.txt");
	system("rm RecWin.txt");
	system("rm aug_chol.txt");

	//Initializing parameters//
	int len = 1000; //Length of Array of Output/Input-Iteration function
	real_num init_covar = 80;
	real_num ff = 0.97; 
	int window_size =50;

	int   l=3500,p=2000, m=1500, n=800;
    float l1=0.8, l2=6.3, l3=3.4, l4=4;
	len=l;
	//Initializing random signal
	Vec Y1=Vec::Zero(len); 
	real_num Y2[len]; 
	ofstream SignalFile;
	SignalFile.open ("Signal.txt");
	
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(-0.5, 0.5);



	for (int i = 0; i < n; i++) {
        Y1[i]= l1 + dis(gen);
        Y2[i]= l1 + dis(gen);
		SignalFile<<Y1[i]<<'\n';
	}
	for (int i = n; i <m; i++) {
        Y1[i]= l2 + dis(gen);
        Y2[i]= l2 + dis(gen);
		SignalFile<<Y1[i]<<'\n';
	}
	for (int i = m; i < p; i++) {
        Y1[i]= l3 + dis(gen);
        Y2[i]= l3 + dis(gen);
		SignalFile<<Y1[i]<<'\n';
	}
	for (int i = p; i < l; i++) {
        Y1[i]= l4 + dis(gen);
        Y2[i]= l4 + dis(gen);
		SignalFile<<Y1[i]<<'\n';
	}

	SignalFile.close();
	//const float *Signal_Y=Y;
	//Call Estimating Functions
	Rec_Window_Method( Y2, len,ff, window_size , init_covar);
	Forgeting_Factor_Method( Y2,ff, len, init_covar);
	Aug_Cholesky_Method(window_size, init_covar, Y1, 10 , len );
	
return 0;

}

