#include <cmath>
#include <cstdlib>
#include <random>
#include <fstream>
#include "C:\Users\nicks\rls\source\RLS_Estimation_Object.h"

using namespace RLS;
using namespace arma;
using namespace std;



int main() {
	//Initializing parameters//
	const int len = 10000; //Length of Array of Output/Input
	double init_covar = 100000.;
	//Generate Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };
	

	double counter = 0.;
	const double mean = 0.0;
	const double stddev = 10;
	std::default_random_engine generator;
	std::normal_distribution<double> dist(mean, stddev);


	for (int i = 0; i < 50; i++) {
		Y[i] = counter; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}
	for (int i = 50; i < 100; i++) {
		Y[i] = 50.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}

	for (int i = 100; i < 150; i++) {
		Y[i] = counter + 50.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}

	for (int i = 150; i < 200; i++) {
		Y[i] = 20.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}

	for (int i = 200; i < 250; i++) {
		Y[i] = (i-198.)*(i-198.) + 2.*(i-198.) + 1.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}

	// Add Gaussian noise and Write Output//
	ofstream outFile0;
	outFile0.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test_Output.txt");
	for (int i = 0; i < 250; i++) {
		Y[i] = Y[i] + dist(generator);
		outFile0 << Y[i] << " ";
	}
	outFile0.close();
	

	

	//Test the "Set" and "Get" Functions in the class//
	RLS_Estimator<double, 2> Test_Alg(1., 10000);
	//Do some updates to notice changes//
	for (int i = 0; i < 10; i++) {
		Test_Alg.update_par(Y[i]); // Update parameters 
	}
	cout << "Testing Functions : " <<  endl;
	cout << "Parameter function : " << endl;
	Test_Alg.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Alg.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Alg.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Alg.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Change some parameters with "Set" and test them again//
	Test_Alg.setLambda(0.5);
	Test_Alg.setCovariance(5000);
	cout << "New forgetting factor : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "New initial covariancce : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Test "reset" function//
	cout << "Testing reset function : " << endl; 
	Test_Alg.reset();
	cout << "Parameter function : " << endl;
	Test_Alg.getEstimatedParameters().print();
	cout << '\n';
	cout << "Covariance function : " << endl;
	Test_Alg.getCovarianceMat().print();
	cout << '\n';
	cout << "Gain function : " << endl;
	Test_Alg.getGains().print();
	cout << '\n';
	cout << "Iteration function : " << Test_Alg.getIterations() << endl;
	cout << '\n';
	cout << "Forgetting factor function : " << Test_Alg.getLambda() << endl;
	cout << '\n';
	cout << "Initial Covariance function : " << Test_Alg.getCovar() << endl;
	cout << '\n';

	//Initialize some different parameters to check differences//
	//Test 1 : 2 Parameters and forgetting factor 1

	RLS_Estimator<double, 2> Alg_1(1., init_covar);

	cout << "First Test : We expect parameters for each area to have errors as every single input is taken in regard" << endl;
	cout << "and doesn't accurately represent the area." << endl;

	for (int i = 0; i < 200; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_1.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i <<" : " << endl;
		Alg_1.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		
	}
	
	//Test 2 : 2 Parameters and forgetting factor 0.9

	RLS_Estimator<double, 2> Alg_2(0.9, init_covar);

	cout << "Second Test : We expect parameters to be time-variant and change depending on the output" << endl;
	cout << "In our case , we should see that it resembles the polynomial we put as input." << endl;

	//Writing Data from Parameters//
	ofstream outFile1;
	ofstream outFile2;
	outFile1.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test2_Param_a0.txt");
	outFile2.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test2_Param_a1.txt");

	for (int i = 0; i < 250; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_2.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_2.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		outFile1 << Alg_2.getEstimatedParameters().row(0) << " ";
		outFile2 << Alg_2.getEstimatedParameters().row(1) << " ";
	}
	outFile1.close();
	outFile2.close();

	//Test 3 : 3 Parameters and forgetting factor 0.9

	RLS_Estimator<double, 3> Alg_3(0.9, init_covar);

	cout << "Third Test : We are estimating everything as a 2nd order Polynomial"<< endl;

	//Writing Data from Parameters//
	ofstream outFile3;
	ofstream outFile4;
	ofstream outFile5;
	outFile3.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test3_Param_a0.txt");
	outFile4.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test3_Param_a1.txt");
	outFile5.open("C:/Users/nicks/rls/MATLAB/TXT-Files/Test3_Param_a2.txt");
	for (int i = 0; i < 250; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_3.update_par(Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_3.getEstimatedParameters().print(); //Print Matrix 
		cout << '\n';
		outFile3 << Alg_3.getEstimatedParameters().row(0) << " ";
		outFile4 << Alg_3.getEstimatedParameters().row(1) << " ";
		outFile5 << Alg_3.getEstimatedParameters().row(2) << " ";
	}
	outFile3.close();
	outFile4.close();
	outFile5.close();

	cout << "Takes much longer to converge with this method to the correct input,if at all. " << endl;
	cout << "Need to look for big shifts in parameters,to validate the change,ie 50 ->6500 at that moment";
	cout << " or 0.0002 -> -5.8";


	return 0;
}