#ifndef RLS_Estimation
#define RLS_Estimation

//Directories//
#include <armadillo>
#include <iostream>
#include <vector>
using namespace std;
using namespace arma;


namespace RLS {


	template <typename T, int N>
	class RLS_Estimator {
	public:
		typedef Mat<T>(N, N) Two_D;
		typedef Mat<T>(N, 1) One_D;


	private:
		int np; //Number of Parameters - Order of Polynomial
		double lambda; //Forgetting Factor
		double init_covar; //Initial Covariance (Preferably large to declare indiffirence at first iterations//
		One_D theta_matrix; //Matrix of Parameters//
		Two_D P_matrix; //Covariance Matrix//
		double K; //Gain //
		unsigned long long num_update; //Number of updates//

	public: //Set Defaults Values at Construction of Object//
		RLS_Estimator(double lam, double init)
			:np(N), lambda(1.0),
			init_covar(init),
			theta_matrix(One_D.zeros()),
			P_matrix(Two_D.eye()),
			K(One_D.zeros()),
			num_update(0) {
			setLambda(lam);
			setCovariance(init);
			P_matrix = P_matrix * init;
		}

		// Update of Parameters at Time (Last) and with New data (data)
		void update_par(double time, T data) {
			mat phi_matrix(1, np, fill::zeros);
			phi_matrix.col(0) = 1;

			for (int i = 1; i < np; i++) {
				phi_matrix.col(i) = pow(time,i);
			}

			phi_matrix_T = phi_matrix.t();

			K = phi_matrix * (1. / l) * P_matrix * phi_matrix_T + 1.;
			K = inv(K1);
			K = (1. / l) * P_matrix * phi_matrix_T * K1;


			//CALCULATION OF  NEW PARAMETERS//

			theta_matrix = theta_matrix + K * (data - phi_matrix * theta_matrix);

			//CALCULACTION OF NEW COVARIANCE MATRIX//
			P_matrix = ((1. / l) * P_matrix - K * phi_matrix * P_matrix);
		};

		void setLambda(double lam) {
			if ((lam > 0) && (lam <= 1.0)) {
				lam_ = lam;
				lam_inv_ = 1.0 / lam_;
			}
			else {
				throw std::invalid_argument(
					"Invalid forgetting factor (0 < lambda <= 1).");
			}
		};

		void setCovariance(double init) {
			if (init > 0.0) {
				init_covar = init;
			}
			else {
				throw std::invalid_argument(
					"Invalid covariance matrix gain factor (delta > 0).");
			}
		};

		const One_D& estimatedParameters() const noexcept { return theta_matrix; };
		const Two_D& CovarianceMat() const noexcept { return P_matrix; };

		void reset() noexcept {
			theta_matrix = One_D.zeros(np);
			P_matrix = Two_D.eye(np, np) * init_covar;
			K = One_D.zeros(np);
			err_ = 0.0;
			num_update = 0;
		}
	};

};

#endif






























