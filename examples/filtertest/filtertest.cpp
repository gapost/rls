#include <iostream>
#include <cstdlib>
#include <random>
#include <fstream>
#include <typeinfo>
#include "/home/tinaa/Documents/Intership/rls/source/RLS_Estimation_Object.h"
#include <stdio.h>
#include <stdlib.h>

using namespace RLS;
using namespace arma;
using namespace std;
using namespace std;

void testFF(const double* Y, int len,int factors ,double init_covar)

{  	
    vec Reg(factors, fill::eye);

   	RLS_Estimator<double> Test_Gen(factors,0.99,init_covar);  // (Number of factors ,Forgeting factor, Initial Covarience function)
    // Update of Parameters with every new data// 
	for (int i = 0; i < len; i++) { //for each sample i
		for (int j=0; j<factors; j++){	//update polinomial factors
			Reg(j)=pow(i,j);
		}
		Test_Gen.update_par(Reg,Y[i]); 
	}

	//Print output of parameters from the generalized class in .txt file
    ofstream out;
    out.open ("build/FF.txt");

	for (int i = 0; i < len; i++) {
		for (int j=0; j<factors; j++){	//update polinomial factors
			Reg(j)=pow(i,j);
		}
        Test_Gen.update_par(Reg,Y[i]); // Update parameters in respect to Input and Output

		out << Test_Gen.getEstimatedOutput(Reg)<<'\t' << (Test_Gen.getEstimatedParameters()(0)) << "\t" << (Test_Gen.getEstimatedParameters()(1)) << endl;	}

    out.close();
    }

void testWindow(const double* Y, int len,int factors ,double init_covar)
    
    {   
        vec Reg(factors, fill::eye);
        
        BlockRLS<double> Test_Block(factors,1,50, init_covar);
        	for (int i = 0; i < len; i++) { //for each sample i
		for (int j=0; j<factors; j++){	//update polinomial factors
			Reg(j)=pow(i,j);
		}
		Test_Block.update_par(Reg,Y[i]); 
	}

	//Check output of parameters from the generalized class
 	ofstream a0;
    ofstream a1;
    ofstream out;

    out.open ("build/RecWin.txt");

	for (int i = 0; i < len; i++) {
		for (int j=0; j<factors; j++){	//update polinomial factors
			Reg(j)=pow(i,j);
		}
        Test_Block.update_par(Reg,Y[i]); // Update parameters in respect to Input and Output

		out << Test_Block.getEstimatedOutput(Reg)<<'\t' << (Test_Block.getEstimatedParameters()(0)) << "\t" << (Test_Block.getEstimatedParameters()(1)) << endl;	}
    out.close();

    }

int main() {
	
	//Initializing parameters//
	const int len = 1500; //Length of Array of Output/Input-Iteration function
	double init_covar = 100.;//??
	int factors=2;

	//Generate Data for Validation
	double Y[len]; 
	ofstream SignalFile;
	SignalFile.open ("build/Signal.txt");

	for (int i = 0; i < 2*len/5; i++) {
		Y[i]=1.0 + 0.1*(2.0*rand()/RAND_MAX-1.0);
		SignalFile<<Y[i]<<'\n';
	}

	for (int i = 2*len/5; i < len; i++) {

		Y[i]=5.0 + 0.1*(2.0*rand()/RAND_MAX-1.0);
		SignalFile<<Y[i]<<'\n';
	}
	SignalFile.close();

    testFF( Y, len,factors, init_covar);
	testWindow( Y, len,factors, init_covar);

return 0;
}