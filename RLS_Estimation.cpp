#include <cmath>
#include <cstdlib>
#include <random>
#include "C:\armadillo-11.2.1\examples\Header.h"

using namespace RLS;
using namespace arma;

/*int main() {
	//INITIALIZATION OF PARAMETERS//
	
	const int len = 10000; //Length of Array of Output/Input
	const int np = 2; //Number of parameters of system to be calculated
	double l = 0.65; //Forgetting factor

	//Generate Random Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };
	double theta[np] = { 0. };

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

	// Add Gaussian noise
	for (int i = 0; i < 200; i++) {
		Y[i] = Y[i] + dist(generator);
	}

	double P[np * np] = { 0. };
	for (int i = 0; i < np * np; i += (np + 1)) {
		P[i] = 10000.; //Starting P large to declare indifference
	}


	mat P_matrix(P, np, np);
	mat Y_matrix(Y, len, 1);
	mat X_matrix(X, len, 1);

	//END OF INITIALIZATION OF PARAMETERS//
	//INITIALIZATION OF PHI//
	mat K1;
	mat K;
	bool check;
	mat theta_matrix(theta, np, 1);
	double phi[np] = { 0. };

	for (int i = 0; i < 200; i++) {

		phi[0] = 1;
		phi[1] = X[i];

		mat phi_matrix(phi, 1, np);
		mat phi_matrix_T(phi, np, 1);

		//END OF INITIALIZATION OF PHI AND PHI TRANSPOSE//

		//CALCULATION OF GAIN VECTOR 

		K1 = phi_matrix * (1. / l) * P_matrix * phi_matrix_T + 1.;
		K1 = inv(K1);
		cout << "Here is gains : " << endl; 
		K1.print();
		K = (1. / l) * P_matrix * phi_matrix_T * K1;


		//CALCULATION OF PARAMETERS//

		theta_matrix = theta_matrix + K * (Y[i] - phi_matrix * theta_matrix);

		//CALCULACTION OF NEW COVARIANCE MATRIX//
		P_matrix = ((1. / l) * P_matrix - K * phi_matrix * P_matrix);

		cout << "Iteration number " << i << endl;
		theta_matrix.print();
		cout << "Here is output:" << Y[i] << endl;
		cout << "\n" << endl;
	}


	return 0;
}*/

int main(){
	const int len = 10000; //Length of Array of Output/Input
	const int np = 2; //Number of parameters of system to be calculated
	double l = 0.65; //Forgetting factor
	
	

	//Generate Random Data for Validation
	double Y[len] = { 0. };
	double X[len] = { 0. };
	double theta[np] = { 0. };

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

	// Add Gaussian noise
	for (int  i = 0; i < 200; i++) {
		Y[i] = Y[i] + dist(generator);
	}
	double init_covar = 100000.;
	RLS_Estimator<double, 2> a(l, init_covar);
	Mat<double> A(5, 5, fill::zeros);
	
	
	
	for (int i = 0; i < 200; i++) {
		a.update_par(X[i], Y[i]);
		cout << "Here are the estimated parameters: " << endl;
		a.estimatedParameters().print();
		cout << '\n';
	}
	



	return 0;
}