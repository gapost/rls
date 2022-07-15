#ifndef RLS_Estimation
#define RLS_Estimation

//Directories//
#include <armadillo>
#include <iostream>
#include <vector>

using namespace arma;


namespace RLS {


	template <typename T, const int N>
	class RLS_Estimator {
	public:
		typedef Mat<T> Type_Mat;


	private:
		const int np; //Number of Parameters - Order of Polynomial
		double lambda; //Forgetting Factor
		double init_covar; //Initial Covariance (Preferably large to declare indiffirence at first iterations//
		Type_Mat phi_matrix;
		Type_Mat theta_matrix; //Matrix of Parameters//
		Type_Mat P_matrix; //Covariance Matrix//
		Type_Mat K; //Gain Vector//
		Type_Mat temp;
		unsigned long long num_update; //Number of updates//

	public: //Set Defaults Values at Construction of Object//
		RLS_Estimator(double lam, double init)
			:np(N), lambda(1.0),
			init_covar(init),
			theta_matrix(Type_Mat(N,1,fill::zeros)),
			P_matrix(Type_Mat(N,N,fill::eye)),
			phi_matrix(Type_Mat(1, N,fill::zeros)),
			K(Type_Mat(N,1,fill::zeros)),
			temp(Type_Mat(N,1,fill::zeros)),
			num_update(0) {
			setLambda(lam);
			setCovariance(init);
			P_matrix = P_matrix * init;
		}

		// Update of Parameters with New data (data)
		void update_par(T data) {

			phi_matrix.col(0) = 1.;
			for (int i = 1; i < N; i++) {
				phi_matrix.col(i) = phi_matrix.col(i - 1) * num_update;
			}
			temp = P_matrix * trans(phi_matrix);
			K = temp / (dot(phi_matrix,temp) + lambda);


			//CALCULATION OF  NEW PARAMETERS//

			theta_matrix += K * (data - phi_matrix * theta_matrix); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//CALCULATION OF NEW COVARIANCE MATRIX//
			temp = trans(temp);
	
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < i+1; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					if (j != i) {
						P_matrix(j, i) = P_matrix(i, j);
					}
				}
			}
			P_matrix /= lambda;
	
			num_update += 1;
		};
		//"Set" Functions//
		void setLambda(double lam) {
			if ((lam > 0) && (lam <= 1.0)) {
				lambda = lam;
			
			}
			else {
				throw std::invalid_argument(
					"Forgetting factor must be  (0 < lambda <= 1).");
			}
		};

		void setCovariance(double init) {
			if (init > 0.0) {
				init_covar = init;
			}
			else {
				throw std::invalid_argument(
					"Covariance must be positive.");
			}
		};

		//"Get" Functions//
		const Type_Mat& getEstimatedParameters() const noexcept { return theta_matrix; }
		const Type_Mat& getCovarianceMat() const noexcept { return P_matrix; }
		const Type_Mat& getGains() const noexcept { return K; }
		int getIterations() const noexcept { return num_update;  }
		double getEstimatedOutput() const noexcept { return dot(phi_matrix, theta_matrix); }
		const double getLambda() const noexcept { return lambda; }
		const double getCovar() const noexcept { return init_covar; }

		void reset() noexcept {
			theta_matrix = Type_Mat(1, N, fill::zeros);
			P_matrix = Type_Mat(N, N,fill::eye) * init_covar;
			K = Type_Mat(1, N, fill::zeros);
			phi_matrix = Type_Mat(1, N, fill::zeros);
			num_update = 0;
		};
	};

};

#endif






























