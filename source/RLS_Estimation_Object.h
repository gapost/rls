#ifndef RLS_Estimation
#define RLS_Estimation

//Directories//
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;


namespace RLS {

    template <typename T>
	class RLS_Estimator {
	public:
		static_assert(std::is_arithmetic_v<T>,"The estimator doesn't support string inputs.");
		typedef Matrix <T , Dynamic, Dynamic > Type_Mat;
		typedef Matrix< T, Dynamic, 1 > Type_Vec;

	protected:
		int np; //Number of Parameters - Order of Polynomial
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
        RLS_Estimator(int n, double ff, double init) // n: number of factors, ff: forgrting factor, init: Initial Covarience function
            : np(n), lambda(1.0),
			init_covar(init),
            theta(n),			// new parameters
            P_matrix(n,n),		// [1 0 ; 0 1]
            K(n),				// gain vector
            temp(n),			 
			num_update(0),		// Update number of iterations
			cost(0),			// Cost function value
            error(0)			// Error function value
        {   
            theta.setZero();
			K.setZero();
			P_matrix.setIdentity();
			temp.setZero();
			setLambda(ff);
			setCovariance(init);
			P_matrix = P_matrix * init;		// [init 0 ; 0 init]

		}

		// Update of Parameters with New data (data)
        void update_par(const Type_Vec& phi, T data) {
			
			//cout << "\nphi:\n" << phi << endl;

			temp = P_matrix * phi;
		
			//Set gain vector
			K = temp / (phi.dot( temp) + lambda);

			//Calculate error and cost function values
			error = data - phi.dot(theta);
			cost = lambda * cost + error * error;

			//Calculation of new parameters//
			theta += K * (error); //Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
			//cout<<"\nafto to ypoxthonio theta=\n"<<theta<<endl;
			//cout<<"\nafto to ypoxthonio error=\n"<<error<<endl;
			//cout<<"\nafto to ypoxthonio phi=\n"<<phi<<endl;
			//cout << "\nnew parameters:\n" << theta <<endl;
			
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
			}
			
            num_update++; //Update number of iterations
		};
		//"Set" Functions//
		void setLambda(double ff) {
			if ((ff > 0) && (ff <= 1.0)) {
				lambda = ff;
			
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
		const Type_Vec & getEstimatedParameters() const noexcept { return theta; }
		const Type_Mat & getCovarianceMat() const noexcept { return P_matrix; }
		const Type_Vec & getGains() const noexcept { return K; }
		int getIterations() const noexcept { return num_update;  }
		double getEstimatedOutput(Type_Vec& x) const noexcept { 
			//cout<<"\nx=\n"<<x<<"\ntheta=\n"<<theta<<"\n it returns:"<<x.dot(theta)<<endl;
			return x.dot(theta); } //x=Reg!, theta=new parameters
		double getCost() const noexcept { return cost; }
		double getError() const noexcept { return error; }
        double getLambda() const noexcept { return lambda; }
        double getCovar() const noexcept { return init_covar; }
		
		//Reset Function
		void reset() noexcept {
            theta.setZero();
            P_matrix.setIdentity();
            P_matrix *= init_covar;
            K.setZero();
			cost = 0;
			error = 0;
			num_update = 0;
		};
	};

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        using RLS_Estimator<T>::lambda;		 //Forgetting Factor
        using RLS_Estimator<T>::cost;		 //Cost function
        using RLS_Estimator<T>::error;		 //Error value
        using RLS_Estimator<T>::theta;		 //Matrix of Parameters
        using RLS_Estimator<T>::P_matrix;	 //Covariance Matrix
        using RLS_Estimator<T>::K;			 //Gain Vector
        using RLS_Estimator<T>::temp;
        using RLS_Estimator<T>::num_update;	 //Number of updates
	public:
        BlockRLS(int n, double ff, int win, double init) // 
            : RLS_Estimator<T>( n, ff, init), //why???
			window(win),
            pout(win + 1),
            pin(n, win + 1)
		 {
            // initialize regressors according to Zhang 2004
            for (int i = 0; i < n; i++)
                pin(i, win - n + i) = 1. / sqrt(init);

		}
        void update_par(const Type_Vec& phi, T data)
        {   
			typedef Matrix <float , Dynamic, Dynamic > Type_Mat;


            for (int i=0; i<window; i++){
                pin.col(i)=pin.col(i+1);
            }
            pin.col(window)=phi;

            for (int i=0; i<window; i++){
			pout(i)=pout(i+1);
			}
			pout(window)=data;
			
            temp.setZero();
			temp = P_matrix * phi;
			if ((phi.dot(temp) + lambda)==0){
			//cout<<"\nwhat if its zero:\n"<<(phi.dot(temp) + lambda)<<endl;
			}			
			K = (temp) / (phi.dot(temp) + lambda); 	// Gain Vector
			error = data - theta.dot(phi);			// Error Value

			//CALCULATION OF  NEW PARAMETERS//
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

			num_update += 1; 							//Update number of iterations
			temp = P_matrix * pin.col(0);
			K = temp / (lambda - temp.dot(pin.col(0)));

			error = pout(0) - theta.dot(pin.col(0));
			//CALCULATION OF  NEW PARAMETERS//
			//cost = lambda * cost + error * error;

            theta -= K * error; 						//Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
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
            pin.setZero();
            for (int i = 0; i < np; i++) {
                pin(i, window - np + i) = 1./sqrt(init_covar);
			}
            pout.setZero();
		};


	};

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    //Subclass of RLS_Estimator that implemenets a polynomial estimation
/* template <typename T, class Estimator>
    class PolyRLS : public Estimator
    {
    };

  template<typename T>
    class PolyRLS<T, RLS_Estimator<T>> : public RLS_Estimator<T>
    {
        typename RLS_Estimator<T>::Type_Vec phi;
      public:
        PolyRLS(int n, double ff, double init) : // n: number of factors, lam: forgrting factor, init: Initial Covarience function
            RLS_Estimator<T>(n, ff, init),
            phi(n) 
            {	
				cout << "yes" << endl;

                for (int i=0; i<n; i++)
			    phi(i)=1;
			}

        void update_par(T data) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
				
				cout << "this:" << this << endl;
            
			}

            RLS_Estimator<T>::update_par(phi, data);
        
			cout << "phi dot theta =" << dot(phi, this->theta) << endl;

        }
		double getEstimatedOutput() const noexcept { 
			
			cout << "this:" << this << endl;

			return phi.dot(this->theta); }
    };    
	
	template<typename T>
    class PolyRLS<T, BlockRLS<T>> : public BlockRLS<T>
    {
        typename BlockRLS<T>::Type_Vec phi;
      public:
        PolyRLS(int n, int w, double init) :
            BlockRLS<T>(n, w, init), // n: number of factors, w: window size, init: Initial Covarience function
            phi(n)
            {
                for (int i=0; i<n;i++)
			    phi(i)=1;
			}
        void update_par(T data) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
            }

            BlockRLS<T>::update_par(phi, data);

        }
        double getEstimatedOutput() const noexcept { 
			
			cout << "phi dot theta =" << dot(phi, this->theta) << endl;

			return phi.dot(this->theta); }
    };
	 */

} // namespace RLS

#endif