/*#include <cstdlib>
#include <random>
#include <fstream>
#include "C:\armadillo-11.2.1\examples\RLS_Estimation_Object.h"

using namespace RLS;
using namespace arma;
using namespace std;

int main(){
	const int len = 10; //Length of Array of Output/Input
	double init_covar = 100000.;
	//Generate Data for Validation
	double Y[4*len] = { 0. };
	double X[4*len] = { 0. };

	double counter = 0.;
	const double mean = 0.0;
	const double stddev = 5;
	std::default_random_engine generator;
	std::normal_distribution<double> dist(mean, stddev);

	Y[0] = 0;
	for (int i = 1; i < 2*len; i++) {
		Y[i] = Y[i-1] - 1; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}
	for (int i = 2*len; i < 4 * len; i++) {
		Y[i] = Y[i-1] + 1; //OUTPUTS OF SYSTEM
		X[i] = counter; //INPUTS OF SYSTEM
		counter += 1.;
	}
	vec Reg(3, fill::eye);
	BlockRLS<double> Test_Block(3,1, len, 10000.);
	for (int i = 0; i < 4 * len; i++) {
		Reg(1) = i;
		Reg(2) = i * i;
		cout << "True output: " << Y[i] << " |";
		cout << "Estimated output: " << Test_Block.getEstimatedOutput(Reg) << " |";
		Test_Block.update_par(Reg, Y[i]);
		cout << "Estimated parameters at time " << i << " : " << " |";
		trans(Test_Block.getEstimatedParameters()).print(); //Print Matrix 
		cout << '\n';
	}
}*/