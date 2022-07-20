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
		typedef Col<T> Type_Vec;

	protected:
		const int np; //Number of Parameters - Order of Polynomial
		double lambda; //Forgetting Factor
		double init_covar; //Initial Covariance (Preferably large to declare indiffirence at first iterations//
		double cost;
		double error;
		Type_Vec phi;
		Type_Vec theta; //Matrix of Parameters//
		Type_Mat P_matrix; //Covariance Matrix//
		Type_Vec K; //Gain Vector//
		Type_Vec temp;
		unsigned long long num_update; //Number of updates//

	public: //Set Defaults Values at Construction of Object//
		RLS_Estimator(double lam, double init)
			:np(N), lambda(1.0),
			init_covar(init),
			theta(Type_Vec(N,fill::zeros)),
			P_matrix(Type_Mat(N,N,fill::eye)),
			phi(Type_Vec(N,fill::zeros)),
			K(Type_Vec(N,fill::zeros)),
			temp(Type_Vec(N,fill::zeros)),
			num_update(0),
			cost(0),
			error(0){
			setLambda(lam);
			setCovariance(init);
			P_matrix = P_matrix * init;
		}

		// Update of Parameters with New data (data)
		void update_par(vec& x, T data) {

			for (int i = 0; i < N; i++) {
				phi(i) = x(i);
			}

			temp = P_matrix * phi;

			K = temp / (dot(phi, temp) + lambda);

			//CALCULATION OF  NEW PARAMETERS//
			error = data - dot(phi, theta);
			cost = lambda * cost + error * error;
			
			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//CALCULATION OF NEW COVARIANCE MATRIX//

			for (int i = 0; i < N; i++) {
				P_matrix(i, i) -= K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					//Matrix is symmetric - assign values for less computations
					P_matrix(j, i) = P_matrix(i, j);

				}
			};
			
			num_update += 1; //Update number of iterations
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
		const Type_Vec& getEstimatedParameters() const noexcept { return theta; }
		const Type_Mat& getCovarianceMat() const noexcept { return P_matrix; }
		const Type_Vec& getGains() const noexcept { return K; }
		int getIterations() const noexcept { return num_update;  }
		double getEstimatedOutput() const noexcept { return dot(phi, theta); }
		double getCost() const noexcept { return cost; }
		double getError() const noexcept { return error; }
		const double getLambda() const noexcept { return lambda; }
		const double getCovar() const noexcept { return init_covar; }

		void reset() noexcept {
			theta = Type_Vec(N, fill::zeros);
			P_matrix = Type_Mat(N, N,fill::eye) * init_covar;
			K = Type_Vec(N, fill::zeros);
			phi = Type_Vec(N, fill::zeros);
			temp = Type_Vec(N,fill::zeros);
			num_update = 0;
		};
	};
	template <typename T, const int N>
	class RLS_Estimator_Poly : public RLS_Estimator<T, N> {
	public:
		typedef Mat<T> Type_Mat;

	public:
		RLS_Estimator_Poly(double lam, double init)
			: RLS_Estimator<T, N>(lam, init) {}
		void update_par(T data) {
			
			phi(0) = 1.;
			for (int i = 1; i < N; i++) {
				phi(i) = phi(i - 1) * num_update;
			}
	
			temp = P_matrix * phi;
	
			K = temp / (dot(phi, temp) + lambda);

			error = data - dot(phi, theta);

			//CALCULATION OF  NEW PARAMETERS//

			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//CALCULATION OF NEW COVARIANCE MATRIX//

			for (int i = 0; i < N; i++) {
				P_matrix(i, i) -= K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					//Matrix is symmetric - assign values for less computations
					P_matrix(j, i) = P_matrix(i, j);

				}
			};
			cost *= lambda;
			cost += error * error;
			num_update += 1; //Update number of iterations
		};
		

	};
};
#endif






























