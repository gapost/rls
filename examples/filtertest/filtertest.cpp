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
		RegF(0)=pow(i,0);
		RegF(1)=pow(i,1);
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
		RegW(0)=pow(i,0);
		RegW(1)=pow(i,1);
		Test_Block.update_par(RegW,Y[i],i); 
		out <<Test_Block.getEstimatedOutput(RegW)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	

	}
   	out.close();
}

void Aug_Cholesky_Method(int win, real_num init_covar, Vec Y, int start, int l )
{
	ofstream out;
	out.open ("aug_chol.txt");

	Mat Phi(start,2);
	for (int i =0; i<start; i++){
		Phi(i,0)=1;
		Phi(i,1)=i;
	}
    
	// Matrix used to update llt
    Vec v_up=Vec::Zero(3,1);
	Vec v_down=Vec::Zero(3,1);
    v_up(0,0)=1;   
	v_down(0,0)=1;

	Augmented_Cholesky_RLS_Estimator<real_num> Test_Block (2 , win);

    // Update Cholesky matrix and compute Theta again
    for(int j=1; j<l; j++){ 
        if(j>win){
            v_up(1,0)=j;
            v_up(2,0)=Y(j);
            v_down(1,0)=j-win;
            v_down(2,0)=Y(j-win);
			
        }
        else{
            v_up(1,0)=j;
            v_up(2,0)=Y(j);
        }
		//cout<<v_up.topLeftCorner(2, 1)<<"\n\n\n";
   		Test_Block.update_par(Phi(j,2), Y(j) );
		Vec v=v_up.topLeftCorner(2, 1);
		out << Test_Block.getEstimatedOutput(v)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	
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
	int window_size =80;

	int   l=2500, m=1500, n=800;
    float l1=3.8, l2=6.3, l3=3.4;
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
	for (int i = m; i < l; i++) {
        Y1[i]= l3 + dis(gen);
        Y2[i]= l3 + dis(gen);
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

