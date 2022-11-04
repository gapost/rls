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

    template <typename real_num>
	class RLS_Estimator {
	public:
		static_assert(std::is_arithmetic_v<real_num>,"The estimator doesn't support string inputs.");
		typedef Matrix <real_num , Dynamic, Dynamic > Type_Mat;
		typedef Matrix< real_num, Dynamic, 1 > Type_Vec;

	protected:
		int np; 										// Number of Parameters - Order of Polynomial
		real_num lambda; 								// Forgetting Factor
		real_num init_covar; 							// Initial Covariance (Preferably large to declare indiffirence at first iterations//
		real_num cost; 									// Cost function 
		real_num error; 								// Error value
		Type_Vec theta; 								// Matrix of Parameters//
		Type_Mat P_matrix;								// Covariance Matrix//
		Type_Vec K; 									// Gain Vector//
        Type_Vec temp; 									// temporary
		unsigned long long num_update; 					// Number of updates//

	public: 											//Set Defaults Values at Construction of Object//
        RLS_Estimator(int n, real_num ff, real_num init)// n: number of factors, ff: forgrting factor, init: Initial Covarience function
            : np(n), lambda(1.0),
			init_covar(init),
            theta(n),									// new parameters
            P_matrix(n,n),								// [1 0 ; 0 1]
            K(n),										// gain vector
            temp(n),			 
			num_update(0),								// Update number of iterations
			cost(0),									// Cost function value
            error(0)									// Error function value
        {   
            theta.setZero();
			K.setZero();
			P_matrix.setIdentity();
			temp.setZero();
			setLambda(ff);
			setCovariance(init);
			P_matrix = P_matrix * init;					// [init 0 ; 0 init]

		}

		// Update of Parameters with New data (data)
        void update_par(const Type_Vec& phi, real_num data) {
			
			temp = P_matrix * phi;
		
			//Set gain vector
			K = temp / (phi.dot( temp) + lambda);

			//Calculate error and cost function values
			error = data - phi.dot(theta);
			cost = lambda * cost + error * error;

			//Calculation of new parameters//
			theta += K * (error); 						// Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
			
			//Calculation of new covariance//
            for (int i = 0; i < np; i++) {
				P_matrix(i, i) -= K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					P_matrix(j, i) = P_matrix(i, j);	// Matrix is symmetric - assign values for less computations

				}
			}
			
            num_update++; 								// Update number of iterations
		};
		//"Set" Functions//
		void setLambda(real_num ff) {
			if ((ff > 0) && (ff <= 1.0)) {
				lambda = ff;
			
			}
			else {
				throw std::invalid_argument(
					"Forgetting factor must be  (0 < lambda <= 1).");
			}
		};

		void setCovariance(real_num init) {
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
		real_num getEstimatedOutput(Type_Vec& x) const noexcept { 
			return x.dot(theta); } 						// x=time vector (t^0 t^1 .. t^n), theta=new parameters
		real_num getCost() const noexcept { return cost; }
		real_num getError() const noexcept { return error; }
        real_num getLambda() const noexcept { return lambda; }
        real_num getCovar() const noexcept { return init_covar; }
		
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

    template <typename real_num>
    class BlockRLS: public RLS_Estimator<real_num> {
    public:
        int window; 									// block size
        typedef typename RLS_Estimator<real_num>::Type_Mat Type_Mat;
        typedef typename RLS_Estimator<real_num>::Type_Vec Type_Vec;
        Type_Mat pin;
        Type_Vec pout;
        using RLS_Estimator<real_num>::np;
        using RLS_Estimator<real_num>::init_covar;
        using RLS_Estimator<real_num>::lambda;		 	// Forgetting Factor
        using RLS_Estimator<real_num>::cost;		 	// Cost function
        using RLS_Estimator<real_num>::error;		 	// Error value
        using RLS_Estimator<real_num>::theta;		 	// Matrix of Parameters
        using RLS_Estimator<real_num>::P_matrix;		// Covariance Matrix
        using RLS_Estimator<real_num>::K;			 	// Gain Vector
        using RLS_Estimator<real_num>::temp;
        using RLS_Estimator<real_num>::num_update;	 	// Number of updates
	public:
        BlockRLS(int n, real_num ff, int win, real_num init) // 
            : RLS_Estimator<real_num>( n, ff, init), 	// why???
			window(win),
            pout(win + 1),
            pin(n, win + 1)
		 {
            // initialize regressors according to Zhang 2004
            for (int i = 0; i < n; i++)
                pin(i, win - n + i) = 1. / sqrt(init);

		}
        void update_par(const Type_Vec& phi, real_num data, int ADD) // what is phi? And why do I use it to update my matrix?
        {   
			typedef Matrix <real_num , Dynamic, Dynamic > Type_Mat;
			
			// shift elements 

			if (ADD<window){							// if the window is not filled
				for (int i=0; i<ADD; i++){
					pin.col(i)=pin.col(i+1);
				}
				for (int i=0; i<ADD; i++){
					pout(i)=pout(i+1);
				}
				pout(ADD)=data;							// add data in the end of the matrix
	            pin.col(ADD)=phi;						// add phi in the end of the matrix
			}
			else{
				for (int i=0; i<window; i++){
					pin.col(i)=pin.col(i+1);
				}
				for (int i=0; i<window; i++){
					pout(i)=pout(i+1);
				}
				pout(window)=data;						// add data in the end of the matrix
			    pin.col(window)=phi;					// add phi in the end of the matrix
			}
			// end of shift 

            temp.setZero();
			temp = P_matrix * phi;
			
			// UPDATE MATRIX //
			K = temp / (lambda + phi.dot(temp)); 		// Gain Vector
			error = data - theta.dot(phi);				// Error Value
			theta += K * (error); 						// Output is in ascending order , ie: a0 + a1*t + a2*t^2.....

			
            for (int i = 0; i < np; i++) {
				P_matrix(i, i) -= K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) -= K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					P_matrix(j, i) = P_matrix(i, j);	// Matrix is symmetric - assign values for less computations
				}
			}

			// DOWNDATE MATRIX //
			num_update += 1; 							// Update number of iterations
			temp = P_matrix * pin.col(0);				// This line take pout into consideration
			K = temp / (lambda - temp.dot(pin.col(0))); // This line take pout into consideration

			error = pout(0) - theta.dot(pin.col(0));	// This line take pout and pin int consideration
			//CALCULATION OF  NEW PARAMETERS//
			//cost = lambda * cost + error * error;

            theta -= K * error; 						// Output is in ascending order , ie: a0 + a1*t + a2*t^2.....
            for (int i = 0; i < np; i++) {
				P_matrix(i, i) += K(i) * temp(i);
				P_matrix(i, i) /= lambda;
				for (int j = 0; j < i; j++) {
					P_matrix(i, j) += K(i) * temp(j);
					P_matrix(i, j) /= lambda;
					P_matrix(j, i) = P_matrix(i, j);	// Matrix is symmetric - assign values for less computations
				}
			}

		};

		//Get Function
        const typename RLS_Estimator<real_num>::Type_Vec& getpout() const noexcept { return pout; }
        const typename RLS_Estimator<real_num>::Type_Mat& getpin() const noexcept { return pin; }
        int getWindow() const { return window; }

		//Reset Function
		void reset() noexcept {
            RLS_Estimator<real_num>::reset();
            pin.setZero();
            for (int i = 0; i < np; i++) {
                pin(i, window - np + i) = 1./sqrt(init_covar);
			}
            pout.setZero();
		};


	};

//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    template <typename real_num>
    class Augmented_Cholesky_RLS_Estimator: public RLS_Estimator<real_num> {
    public:
        int window; // block size
		typedef Matrix <real_num , Dynamic, Dynamic > Type_Mat;
		typedef Matrix< real_num, Dynamic, 1 > Type_Vec;
		LLT<Type_Mat> llt;
        using RLS_Estimator<real_num>::np;
        using RLS_Estimator<real_num>::init_covar;
        using RLS_Estimator<real_num>::cost;		 	// Cost function
        using RLS_Estimator<real_num>::error;		 	// Error value
		using RLS_Estimator<real_num>::lambda;		 	// Forgetting Factor
        using RLS_Estimator<real_num>::theta;		 	// Matrix of Parameters
        using RLS_Estimator<real_num>::P_matrix;		// Covariance Matrix
        using RLS_Estimator<real_num>::K;			 	// Gain Vector
        using RLS_Estimator<real_num>::temp;
        using RLS_Estimator<real_num>::num_update;	 	// Number of updates
	public:
        Augmented_Cholesky_RLS_Estimator(int n, int win, real_num init, int start, Type_Mat Phi, Type_Vec Y) // n: number of factors, 
														// n: number of factors, init: Initial Covarience function
        : RLS_Estimator<real_num>( n, 1, init), 	// why???
			window(win)
		{   
            theta.setZero();
			K.setZero();
			P_matrix.setIdentity();
			temp.setZero();
			//setLambda(ff);
			//setCovariance(init);
			P_matrix = P_matrix * init;					// [init 0 ; 0 init]
			
			Type_Mat Phi_Aug(start,3);
			Phi_Aug.setZero();
			Phi_Aug.topRightCorner(start,1)=Y;
			Phi_Aug.topLeftCorner(start,2)=Phi;
			Type_Mat A_Aug=Phi_Aug.adjoint()*Phi_Aug;
			llt.compute(A_Aug);
			//L_Aug =llt.matrixL();
		}

        void update_par(Type_Vec v_up,Type_Vec v_down, int ADD )
        {   
			//typedef Matrix <real_num , Dynamic, Dynamic > Type_Mat;
			//Matrix <float , Dynamic, Dynamic > L2 = llt.matrixL();
			llt.rankUpdate(v_up, 1);   
			if(window<ADD){
				llt.rankUpdate(v_down,-1);
			}
		    Type_Mat L_Aug =llt.matrixL();
			theta=(L_Aug.bottomLeftCorner(1,2)).adjoint();
		   	Matrix <float , Dynamic, Dynamic > L=(L_Aug.topLeftCorner(2,2)); //gives error if I initialize it with Type_Mat or real_num
		    (L.adjoint()).triangularView< Upper>().solveInPlace(theta); 



		};

		//Get Function
        //const typename RLS_Estimator<real_num>::Type_Vec& getpout() const noexcept { return pout; }
        //const typename RLS_Estimator<real_num>::Type_Mat& getpin() const noexcept { return pin; }
        int getWindow() const { return window; }

		//Reset Function
		void reset() noexcept {/*
            RLS_Estimator<real_num>::reset();
            pin.setZero();
            for (int i = 0; i < np; i++) {
                pin(i, window - np + i) = 1./sqrt(init_covar);
			}
            pout.setZero();*/
		};


	};
//-----------------------------------------------------------------------------------------------------------------------------------------
    //Subclass of RLS_Estimator that implemenets a polynomial estimation
/* template <typename T, class Estimator>
    class PolyRLS : public Estimator
    {
    };

  template<typename T>
    class PolyRLS<real_num, RLS_Estimator<real_num>> : public RLS_Estimator<real_num>
    {
        typename RLS_Estimator<real_num>::Type_Vec phi;
      public:
        PolyRLS(int n, real_num ff, real_num init) : // n: number of factors, lam: forgrting factor, init: Initial Covarience function
            RLS_Estimator<real_num>(n, ff, init),
            phi(n) 
            {	
				cout << "yes" << endl;

                for (int i=0; i<n; i++)
			    phi(i)=1;
			}

        void update_par(real_numdata) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
				
				cout << "this:" << this << endl;
            
			}

            RLS_Estimator<real_num>::update_par(phi, data);
        
			cout << "phi dot theta =" << dot(phi, this->theta) << endl;

        }
		real_num getEstimatedOutput() const noexcept { 
			
			cout << "this:" << this << endl;

			return phi.dot(this->theta); }
    };    
	
	template<typename T>
    class PolyRLS<real_num, BlockRLS<real_num>> : public BlockRLS<T>
    {
        typename BlockRLS<T>::Type_Vec phi;
      public:
        PolyRLS(int n, int w, real_num init) :
            BlockRLS<T>(n, w, init), // n: number of factors, w: window size, init: Initial Covarience function
            phi(n)
            {
                for (int i=0; i<n;i++)
			    phi(i)=1;
			}
        void update_par(real_numdata) {

            //Set regressors with regards to updates//
            // phi(0) = 1.; this is set in constructor
            for (int i = 1; i < this->np; i++) {
                phi(i) = phi(i - 1) * this->num_update;
            }

            BlockRLS<real_num>::update_par(phi, data);

        }
        real_num getEstimatedOutput() const noexcept { 
			
			cout << "phi dot theta =" << dot(phi, this->theta) << endl;

			return phi.dot(this->theta); }
    };
	 */

} // namespace RLS

#endif
