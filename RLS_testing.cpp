#include <cmath>
#include <cstdlib>
#include <random>
#include "C:\armadillo-11.2.1\examples\Header.h"

using namespace RLS;
using namespace arma;



int main() {
	//Initializing parameters//
	const int len = 10000; //Length of Array of Output/Input
	double init_covar = 100000.;
	//Generate Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };
	

	double counter = 0.;
	const double mean = 0.0;
	const double stddev = 0.1;
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
		Y[i] = pow(i-198,2) + 2*(i-198) + 1.; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}
	// Add Gaussian noise
	/*for (int i = 0; i < 200; i++) {
		Y[i] = Y[i] + dist(generator);
	}*/

	//Initialize some different parameters to check differences//
	//Test 1 : 2 Parameters and forgetting factor 1

	RLS_Estimator<double, 2> Alg_1(1., init_covar);

	cout << "First Test : We expect parameters for each area to be wrong as every single input it taken in regard" << endl;
	cout << "and doesn't accurately represent the area." << endl;

	for (int i = 0; i < 200; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_1.update_par(X[i], Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i <<" : " << endl;
		Alg_1.estimatedParameters().print(); //Print Matrix 
		cout << '\n';
	}

	//Test 2 : 2 Parameters and forgetting factor 0.5

	RLS_Estimator<double, 2> Alg_2(0.5, init_covar);

	cout << "Second Test : We expect parameters to be time-variant and change depending on the output" << endl;
	cout << "In our case , we should see that it resembles the polynomial we put as input." << endl;

	for (int i = 0; i < 200; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_2.update_par(X[i], Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_2.estimatedParameters().print(); //Print Matrix 
		cout << '\n';
	}
	//Test 3 : 4 Parameters and forgetting factor 0.6

	RLS_Estimator<double, 3> Alg_3(0.5, init_covar);

	cout << "Third Test : We are estimating everything as a 2nd order Polynomial"<< endl;

	for (int i = 0; i < 250; i++) {
		cout << "Here is output: " << Y[i] << endl;
		Alg_3.update_par(X[i], Y[i]); // Update parameters in respect to Input and Output
		cout << "Here are the estimated parameters at time " << i << " : " << endl;
		Alg_3.estimatedParameters().print(); //Print Matrix 
		cout << '\n';
	}

	cout << "Takes much longer to converge with this method to the correct input,if at all" << endl;
	cout << "Need to look for big shifts in parameters,to validate the change,ie 50 ->6500 at that moment";
	cout << "or 0.0002 -> -5.8";


	return 0;
}