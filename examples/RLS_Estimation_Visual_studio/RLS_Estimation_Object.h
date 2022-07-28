#ifndef RLS_Estimation
#define RLS_Estimation

//Directories//
#include <armadillo>
#include <iostream>
#include <vector>
#include <stdexcept>

using namespace arma;


namespace RLS {

	template <typename T, const int N>
	class RLS_Estimator {
	public:
		static_assert(std::is_arithmetic_v<T>,"The estimator doesn't support string inputs.");
		typedef Mat<T> Type_Mat;
		typedef Col<T> Type_Vec;

	protected:
		const int np; //Number of Parameters - Order of Polynomial
		double lambda; //Forgetting Factor
		double init_covar; //Initial Covariance (Preferably large to declare indiffirence at first iterations//
		double cost; //Cost function 
		double error; //Error value
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
		void update_par(vec &x, T data) {
			//Pass regressors given in phi matrix
			if (x.n_rows != phi.n_rows) {
				throw std::invalid_argument("Matrix of Regressors needs to be change same length as the parameters");
			}
			for (int i = 0; i < N; i++) {
				phi(i) = x(i);
			}

			temp = P_matrix * phi;
		
			//Set gain vector
			K = temp / (dot(phi, temp) + lambda);

			//Calculate error and cost function values
			error = data - dot(phi, theta);
			cost = lambda * cost + error * error;

			//Calculation of new parameters//
			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//Calculation of new covariance//
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
		double getEstimatedOutput(vec& x) const noexcept { return dot(x, theta); }
		double getCost() const noexcept { return cost; }
		double getError() const noexcept { return error; }
		const double getLambda() const noexcept { return lambda; }
		const double getCovar() const noexcept { return init_covar; }
		
		//Reset Function
		void reset() noexcept {
			theta = Type_Vec(N, fill::zeros);
			P_matrix = Type_Mat(N, N,fill::eye) * init_covar;
			K = Type_Vec(N, fill::zeros);
			phi = Type_Vec(N, fill::zeros);
			temp = Type_Vec(N,fill::zeros);
			cost = 0;
			error = 0;
			num_update = 0;
		};
	};

	//Subclass of RLS_Estimator that implemenets a polynomial estimation of N parameters in regards with time//
	template <typename T, const int N>
	class PolyRLS : public RLS_Estimator<T, N> {
	public:
		PolyRLS(double lam, double init)
			: RLS_Estimator<T, N>(lam, init) {}
		void update_par(T data) {

			//Set regressors with regards to updates//
			phi(0) = 1.;
			for (int i = 1; i < N; i++) {
				phi(i) = phi(i - 1) * num_update;
			}

			temp = P_matrix * phi;
			
			//Set gain vector
			K = temp / (dot(phi, temp) + lambda);

			//Calculate error and cost function values
			error = data - dot(phi, theta);
			cost = lambda * cost + error * error;

			//Calculation of new parameters//
			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
			
			//Calculation of new covariance//

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

		}
		double getEstimatedOutput() const noexcept { return dot(phi, theta); }
	};
	template <typename T, const int N>
	class BlockRLS: public RLS_Estimator<T, N> {
	private: 
		int check;
		int window;
		bool remove;
		Type_Mat pin;
		Type_Vec pout;
		
	public:
		BlockRLS(double lam,int win, double init)
			: RLS_Estimator<T, N>(lam, init),
			  window(win),
			  check(0),
			  pout(Type_Vec(win+1, fill::zeros)),
			  pin(Type_Mat(N, win+1,fill::zeros)),
			  remove(false){}
		void update_par(vec& x, T data){

			for (int i = 0; i < N; i++) {
				phi(i) = x(i);
			}

			pin = shift(pin, -1, 1);
			pout = shift(pout, -1);

			pin.col(window) = phi;
			
			pout(window) = data;

			check += 1; //Raise check for when array is full
			temp = P_matrix * phi;
			if (check > window ) {
				remove = true;
			}

			K = temp / (dot(phi, temp) + lambda);

			//CALCULATION OF  NEW PARAMETERS//
			error = data - dot(phi, theta);

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
			}
			num_update += 1; //Update number of iterations
			if (remove) {
				temp = P_matrix * pin.col(0);

				K = temp / (lambda - dot(pin.col(0), temp));

				error = pout(0) - dot(pin.col(0), theta);
				//CALCULATION OF  NEW PARAMETERS//

				//cost = lambda * cost + error * error;

				theta -= K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
				for (int i = 0; i < N; i++) {
					P_matrix(i, i) += K(i) * temp(i);
					P_matrix(i, i) /= lambda;
					for (int j = 0; j < i; j++) {
						P_matrix(i, j) += K(i) * temp(j);
						P_matrix(i, j) /= lambda;
						//Matrix is symmetric - assign values for less computations
						P_matrix(j, i) = P_matrix(i, j);

					}
				}

			}
			
			
		};
		//Get Function
		const Type_Vec& getpout() const noexcept { return pout; }
		const Type_Mat& getpin() const noexcept { return pin; }
		//Reset Function
		void reset() noexcept {
			theta = Type_Vec(N, fill::zeros);
			P_matrix = Type_Mat(N, N, fill::eye) * init_covar;
			pin = Type_Mat(N, window+1, fill::zeros);
			pout = Type_Vec(window+1, fill::zeros);
			K = Type_Vec(N, fill::zeros);
			phi = Type_Vec(N, fill::zeros);
			temp = Type_Vec(N, fill::zeros);
			lambda = 1.;
			cost = 0;
			check = 0;
			remove = false;
			error = 0;
			num_update = 0;
		};


	};
};
#endif






























