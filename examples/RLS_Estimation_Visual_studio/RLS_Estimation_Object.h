#ifndef RLS_Estimation
#define RLS_Estimation

//Directories//
#include <armadillo>
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>

using namespace arma;


namespace RLS {

    template <typename T>
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
		Type_Vec theta; //Matrix of Parameters//
		Type_Mat P_matrix; //Covariance Matrix//
		Type_Vec K; //Gain Vector//
        Type_Vec temp; // temporary
		unsigned long long num_update; //Number of updates//

	public: //Set Defaults Values at Construction of Object//
        RLS_Estimator(int n, double lam, double init)
            : np(n), lambda(1.0),
			init_covar(init),
            theta(n,fill::zeros),
            P_matrix(n,n,fill::eye),
            K(n,fill::zeros),
            temp(n,fill::zeros),
			num_update(0),
			cost(0),
            error(0)
        {
			setLambda(lam);
			setCovariance(init);
			P_matrix = P_matrix * init;
		}

		// Update of Parameters with New data (data)
        void update_par(const Type_Vec& phi, T data) {

			temp = P_matrix * phi;
		
			//Set gain vector
			K = temp / (dot(phi, temp) + lambda);

			//Calculate error and cost function values
			error = data - dot(phi, theta);
			cost = lambda * cost + error * error;

			//Calculation of new parameters//
			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//Calculation of new covariance//
            for (int i = 0; i < np; i++) {
				P_matrix(i, i) -= K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					//Matrix is symmetric - assign values for less computations
					P_matrix(j, i) = P_matrix(i, j);

				}
			};
			
            num_update++; //Update number of iterations
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
        double getLambda() const noexcept { return lambda; }
        double getCovar() const noexcept { return init_covar; }
		
		//Reset Function
		void reset() noexcept {
            theta.fill(fill::zeros);
            P_matrix.fill(fill::eye);
            P_matrix *= init_covar;
            K.fill(fill::zeros);
			cost = 0;
			error = 0;
			num_update = 0;
		};
	};

    template <typename T>
    class BlockRLS: public RLS_Estimator<T> {
    public:
        int window; // block size
        typedef typename RLS_Estimator<T>::Type_Mat Type_Mat;
        typedef typename RLS_Estimator<T>::Type_Vec Type_Vec;
        Type_Mat pin;
        Type_Vec pout;
        using RLS_Estimator<T>::np;
        using RLS_Estimator<T>::init_covar;
        using RLS_Estimator<T>::lambda; //Forgetting Factor
        using RLS_Estimator<T>::cost; //Cost function
        using RLS_Estimator<T>::error; //Error value
        using RLS_Estimator<T>::theta; //Matrix of Parameters
        using RLS_Estimator<T>::P_matrix; //Covariance Matrix
        using RLS_Estimator<T>::K; //Gain Vector//
        using RLS_Estimator<T>::temp;
        using RLS_Estimator<T>::num_update; //Number of updates
	public:
        BlockRLS(int n, double lam, int win, double init)
            : RLS_Estimator<T>(n, lam, init),
			window(win),
            pout(win + 1, fill::zeros),
            pin(n, win + 1, fill::zeros)
		 {
            // initialize regressors according to Zhang 2004
            for (int i = 0; i < n; i++)
                pin(i, win - n + i) = 1. / sqrt(init);

		}
        void update_par(const Type_Vec& phi, T data)
        {

			pin = shift(pin, -1, 1);
			pout = shift(pout, -1);

			pin.col(window) = phi;
			
			pout(window) = data;

			temp = P_matrix * phi;
			

			K = temp / (dot(phi, temp) + lambda);

			//CALCULATION OF  NEW PARAMETERS//
			error = data - dot(phi, theta);

			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			//CALCULATION OF NEW COVARIANCE MATRIX//
			
            for (int i = 0; i < np; i++) {
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
	
			temp = P_matrix * pin.col(0);

			K = temp / (lambda - dot(pin.col(0), temp));

			error = pout(0) - dot(pin.col(0), theta);
			//CALCULATION OF  NEW PARAMETERS//

			//cost = lambda * cost + error * error;

            theta -= K * error; //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
            for (int i = 0; i < np; i++) {
				P_matrix(i, i) += K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) += K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					//Matrix is symmetric - assign values for less computations
					P_matrix(j, i) = P_matrix(i, j);
				}
			}

		};
		//Get Function
        const typename RLS_Estimator<T>::Type_Vec& getpout() const noexcept { return pout; }
        const typename RLS_Estimator<T>::Type_Mat& getpin() const noexcept { return pin; }
        int getWindow() const { return window; }

		//Reset Function
		void reset() noexcept {
            RLS_Estimator<T>::reset();
            pin.fill(fill::zeros);
            for (int i = 0; i < np; i++) {
                pin(i, window - np + i) = 1./sqrt(init_covar);
			}
            pout.fill(fill::zeros);
		};


	};

    //Subclass of RLS_Estimator that implemenets a polynomial estimation
    template <typename T, class Estimator>
    class PolyRLS : public Estimator
    {
    };

    template<typename T>
    class PolyRLS<T, RLS_Estimator<T>> : public RLS_Estimator<T>
    {
        typename RLS_Estimator<T>::Type_Vec phi;
      public:
        PolyRLS(int n, double lam, double init) :
            RLS_Estimator<T>(n, lam, init),
            phi(n, fill::ones)
        {}
        void update_par(T data) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
            }

            RLS_Estimator<T>::update_par(phi, data);

        }
        double getEstimatedOutput() const noexcept { return dot(phi, this->theta); }
    };

    template<typename T>
    class PolyRLS<T, BlockRLS<T>> : public BlockRLS<T>
    {
        typename BlockRLS<T>::Type_Vec phi;
      public:
        PolyRLS(int n, double lam, int w, double init) :
            BlockRLS<T>(n, lam, w, init),
            phi(n, fill::ones)
        {}
        void update_par(T data) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
            }

            BlockRLS<T>::update_par(phi, data);

        }
        double getEstimatedOutput() const noexcept { return dot(phi, this->theta); }
    };



} // namespace RLS

#endif






























