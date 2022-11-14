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

		typedef Matrix <real_num , Dynamic, Dynamic > Type_Mat;
		typedef Matrix< real_num, Dynamic, 1 > Type_Vec;
		int window;										// block size
		int counter;
		int N;
		Type_Mat Phi_Aug;
		//Type_Vec YY;
		Matrix <float , Dynamic, Dynamic > L;

		Type_Mat v_up;
		Type_Mat v_down;
		Type_Mat L_Aug;
		Type_Mat A_Aug;
		//Type_Mat phi_save;
		//Matrix <float , Dynamic, Dynamic > L;

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
        Augmented_Cholesky_RLS_Estimator(int n, int win)
														// n: number of factors, init: Initial Covarience function
        : RLS_Estimator<real_num>( n, 1, 1), 	
			window(win),
			N(n),
			v_up(N+1,1),
			v_down(N+1,1),
			L_Aug(N+1,N+1),
			A_Aug(window,N+1),
			Phi_Aug(window,N+1)
		{   
			reset();
		}

        void update_par( Type_Vec phi, real_num Y )
        {   
			if (counter < N+1) {
				for (int i=0; i<N; i++){
					Phi_Aug(counter,i)=phi(i);
					}
				Phi_Aug(counter,N)=Y;
			}
			else if (counter == N+1) {
				for (int i=0; i<N; i++){
					Phi_Aug(counter,i)=phi(i);
				}
				Phi_Aug(counter,N)=Y;

				A_Aug=Phi_Aug.adjoint()*Phi_Aug;
				llt.compute(A_Aug);
				L_Aug =llt.matrixL();
				theta=(L_Aug.bottomLeftCorner(1,2)).adjoint();
				L=(L_Aug.topLeftCorner(2,2)); 			//gives error if I initialize it with Type_Mat or real_num
				(L.adjoint()).triangularView< Upper>().solveInPlace(theta); 

			}
			else if ( counter < window) {

				for (int i=0; i<N; i++){
					Phi_Aug(counter,i)=phi(i);
				}
				Phi_Aug(counter,N)=Y;

				for (int i=0; i<N; i++){
					v_up(i,0)=phi(i);
				}
				v_up(N,0)=Y;
				llt=update_llt(llt, v_up); 
				L_Aug =llt.matrixL();
				
				theta=(L_Aug.bottomLeftCorner(1,2)).adjoint();
				L=(L_Aug.topLeftCorner(2,2)); 			//gives error if I initialize it with Type_Mat or real_num
				(L.adjoint()).triangularView< Upper>().solveInPlace(theta); 
			}
			else {
				
				for (int i=0; i<N; i++){
					v_up(i,0)=phi(i);
				}
				v_up(N,0)=Y;

				for (int i =0; i<N; i++){
					v_down(i,0)=Phi_Aug(0,i);
		     	}
				v_down(N,0)=Phi_Aug(0,N);

				llt=up_down_date_llt(llt, v_up, v_down); 
				
				for (int i=0; i<window-1; i++){
					Phi_Aug.row(i)=Phi_Aug.row(i+1);
				}
				for (int i=0; i<N; i++){
					Phi_Aug(window-1,i)=phi(i);
				}
				Phi_Aug(window-1,N)=Y;

				L_Aug =llt.matrixL();
				theta=(L_Aug.bottomLeftCorner(1,2)).adjoint();
				L=(L_Aug.topLeftCorner(2,2)); 			//gives error if I initialize it with Type_Mat or real_num
				(L.adjoint()).triangularView<Eigen::Upper>().solveInPlace(theta); 

			}
		counter+=1;
		};
		LLT<Type_Mat> update_llt(LLT<Type_Mat> llt, const Type_Vec v ){ llt.rankUpdate(v, 1); return llt; }
		LLT<Type_Mat> up_down_date_llt(LLT<Type_Mat> llt, const Type_Vec& v_up, const Type_Vec& v_down ){ 
			llt.rankUpdate(v_up, +1); 
			llt.rankUpdate(v_down, -1);
			return llt; }




		//Get Function
        //const typename RLS_Estimator<real_num>::Type_Vec& getpout() const noexcept { return pout; }
        //const typename RLS_Estimator<real_num>::Type_Mat& getpin() const noexcept { return pin; }
        int getWindow() const { return window; }

		//Reset Function
			void reset() noexcept {
			counter=0;
            theta.setZero();
			//K.setZero();
			//P_matrix.setIdentity();
			//temp.setZero();
			//P_matrix = P_matrix * init;				// [init 0 ; 0 init]
			Phi_Aug.setZero();
			v_up.setZero();
			v_down.setZero();
			L_Aug.setZero();
			A_Aug.setZero();
			//phi_save.setZero(window,N);
			L.setZero(2,2);
			v_up(0,0)=1;
			v_down(0,0)=1;
			//YY.setZero(window);


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
