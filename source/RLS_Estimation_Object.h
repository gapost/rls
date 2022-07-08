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
		unsigned long long num_update; //Number of updates//

	public: //Set Defaults Values at Construction of Object//
		RLS_Estimator(double lam, double init)
			:np(N), lambda(1.0),
			init_covar(init),
			theta_matrix(Type_Mat(N,1,fill::zeros)),
			P_matrix(Type_Mat(N,N,fill::eye)),
			phi_matrix(Type_Mat(1, N,fill::zeros)),
			K(Type_Mat(N,1,fill::zeros)),
			num_update(0) {
			setLambda(lam);
			setCovariance(init);
			P_matrix = P_matrix * init;
		}

		// Update of Parameters at Time (Last) and with New data (data)
		void update_par(T data) {

			phi_matrix.col(0) = 1.;
			for (int i = 1; i < N; i++) {
				phi_matrix.col(i) = phi_matrix.col(i - 1) * num_update;
			}

			//mat phi_matrix(phi, 1, N);

			K = (1. / lambda) * P_matrix * trans(phi_matrix) * (phi_matrix * (1. / lambda) * P_matrix * trans(phi_matrix) + 1.).i();
			//K = K.i();
			//K = (1. / lambda) * P_matrix * phi_matrix_T * K;


			//CALCULATION OF  NEW PARAMETERS//

			theta_matrix = theta_matrix + K * (data - phi_matrix * theta_matrix); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//CALCULATION OF NEW COVARIANCE MATRIX//
			P_matrix = ((1. / lambda) * P_matrix - K * phi_matrix * P_matrix);

			num_update += 1;
		};

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

		const Type_Mat& estimatedParameters() const noexcept { return theta_matrix; }
		const Type_Mat& CovarianceMat() const noexcept { return P_matrix; }

		void reset() noexcept {
			theta_matrix = Type_Mat.zeros(1, np);
			P_matrix = Type_Mat.eye(np, np) * init_covar;
			err_ = 0.0;
			num_update = 0;
		};
	};

};

#endif






























