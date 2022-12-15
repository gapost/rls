#ifndef RLS_Estimation
#define RLS_Estimation

// Directories//
#include <iostream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <Eigen/Dense>

namespace RLS
{
	/**
	 * @brief Base class for RLS estimators
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template < typename T >
	class RlsEstimatorBase
	{
	public:
		typedef T ScalarType;
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixType;
		typedef typename Eigen::Matrix<T, Eigen::Dynamic, 1> VectorType;

	protected:
		int np_;					// number of parameters
		T cost_;					// current cost funtion value
		VectorType theta_;			// current parameter vector
		unsigned long long iter_; 	// # iterations
	public:
		explicit RlsEstimatorBase(int n) : np_(n), theta_(n)
		{
			reset();
		}
		int np() const { return np_; }
		const VectorType& estimatedPar() const { return theta_; }
		T cost() const { return cost_; }
		void setPar(const VectorType &v) { theta_ = v; }
		const unsigned long long & iter() const { return iter_; }
		T estimatedOutput(const VectorType &phi) const { return phi.dot(theta_); }
		void reset() {
			iter_ = 0;
			cost_ = 0;
			theta_.setZero();
		}
	};

	/**
	 * @brief Exp weighted RLS estimator with forgetting factor
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T>
	class ExpWeightedRLS : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		T ff_;				// Forgetting Factor
		MatrixType P;		// Covariance Matrix
		VectorType k;		// Gain Vector
		VectorType temp;	// temporary buffer
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

	public:															
		explicit ExpWeightedRLS(int n, T ff = 0.98, T init_covar = 1.e6) 
		: RlsEstimatorBase<T>(n), ff_(ff),
		P(n, n), k(n), temp(n)
		{
			reset(init_covar);
			setff(ff);
		}

		// Update of Parameters with New data (data)
		void update(const VectorType& phi, const T& data)
		{
			temp = P * phi;

			// Set gain vector
			k = temp / (phi.dot(temp) + ff_);

			// Calculate error and cost function values
			T error = data - phi.dot(theta_);
			cost_ = ff_ * cost_ + error * error;

			// Calculation of new parameters//
			theta_ += k * error; 

			// Calculation of new covariance//
			for (int i = 0; i < np_; i++)
			{
				P(i, i) -= k(i) * temp(i);
				P(i, i) /= ff_;
				for (int j = 0; j < i; j++)
				{
					P(i, j) -= k(i) * temp(j);
					P(i, j) /= ff_;
					P(j, i) = P(i, j); // Matrix is symmetric - assign values for less computations
				}
			}

			iter_++; // Update number of iterations
		};
		//"Set" Functions//
		int setff(T ff)
		{
			if ((ff > 0) && (ff <= 1.0))
			{
				ff_ = ff;
				return 1;
			}
			return 0; 
		};
		T ff() { return ff_; }

		//"Get" Functions//
		const MatrixType &covar() const noexcept { return P; }
		const VectorType &gain()  const noexcept { return k; }

		// Reset Function
		void reset(T init_covar = 1) noexcept
		{
			iter_ = 0;
			theta_.setZero();
			P.setIdentity();
			P *= init_covar;
			k.setZero();
			cost_ = 0;			
		};
	};

	//----------------------------------------------------------------------------//

	/**
	 * @brief RLS estimator with forgetting factor
	 * 
	 * Updating based on P = Q*Q^T decomposition where Q is lower triangular
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T>
	class ExpWeightedRLS2 : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		T ff_;				// Forgetting Factor
		T invsqrtff_;		// 1/sqrt(ff)
		MatrixType Q;		// Square-root Covariance Matrix
		VectorType k;		// Gain Vector
		VectorType f;	// temporary buffers
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

	public:															
		explicit ExpWeightedRLS2(int n, T ff = 0.98, T init_covar = 1.e6) 
		: RlsEstimatorBase<T>(n), ff_(ff),
		Q(n, n), k(n), f(n)
		{
			reset(init_covar);
			setff(ff);
		}

		// Update of Parameters with New data (data)
		void update(const VectorType& phi, const T& data)
		{
			// algorithm from Ljung & Soederstoem (1987) p. 328 
			// f = Q.template triangularView<Eigen::Lower>().adjoint() * phi;
			// T b = ff_ + f.squaredNorm();
			// k = Q * f;
			T b = ff_;
			T error = data;
			for (int i = 0; i < np_; ++i)
			{
				f(i) = 0;
				for (int j = i; j < np_; ++j) f(i) += Q(j,i)*phi(j); // f = Q^T*phi
				b += f(i)*f(i); // b = ff - f^T*f
				k(i) = 0;
				for (int j = 0; j <= i; j++) k(i) += Q(i,j)*f(j); // k = Q*f
				error -= phi(i)*theta_(i);
			}

			
			// Calculate error and cost function values
			// T error = data - phi.dot(theta_);
			cost_ = ff_ * cost_ + error * error;

			// Calculation of new parameters//
			// theta_ += k * (error/b); 

			// Calculation of new covariance //
			// Q is lower triangular !!
			T a = 1/(b + std::sqrt(b*ff_));
			error /= b;
			for (int i = 0; i < np_; i++)
			{
				theta_(i) += k(i)*error;
				Q(i, i) -= a * k(i) * f(i);
				Q(i, i) *= invsqrtff_;
				for (int j = 0; j < i; j++)
				{
					Q(i, j) -= a * k(i) * f(j);
					Q(i, j) *= invsqrtff_;
				}
			}

			iter_++; // Update number of iterations
		};
		//"Set" Functions//
		int setff(T ff)
		{
			if ((ff > 0) && (ff <= 1.0))
			{
				ff_ = ff;
				invsqrtff_ = 1/std::sqrt(ff_);
				return 1;
			}
			return 0; 
		};
		T ff() { return ff_; }

		//"Get" Functions//
		const MatrixType &covar() const noexcept { return Q; }
		const VectorType &gain()  const noexcept { return k; }

		// Reset Function
		void reset(T init_covar = 1) noexcept
		{
			iter_ = 0;
			theta_.setZero();
			Q.setIdentity();
			Q *= std::sqrt(init_covar);
			k.setZero();
			cost_ = 0;			
		};
	};
	
	//----------------------------------------------------------------------------//

	/**
	 * @brief Ring buffer object
	 *
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T>
	class RingBuffer
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;
		typedef Eigen::Block<const MatrixType, Eigen::Dynamic, 1, Eigen::RowMajor> BlockType;
	public:
		/**
		 * @brief Construct a new Ring Buffer object
		 *
		 * @param n : number of column vectors in buffer
		 * @param w : width of each vector
		 */
                RingBuffer(int n, int w) : buff_(w, n), i_(0)
		{
			buff_.setZero();
		}

		int size() const { return buff_.cols(); }
		int width() const { return buff_.rows(); }

		void resize(int n, int w) {
                        buff_ = MatrixType(w,n);
			reset();
		}

		const MatrixType& buff() const { return buff_; }

		template <class V>
		void push(const V &v)
		{
			int n = v.size();
			n = std::min(width(), n);
			for (int j = 0; j < n; ++j)
				buff_(j, i_) = v(j);
			i_++;
            i_ %= size();
		}

		// overloads for scalar values
		void push(const long double &v) { push_scalar(v); }
		void push(const double &v) { push_scalar(v); }
		void push(const float &v) { push_scalar(v); }

		const BlockType first() const {
			return buff_.col(i_); // i_ points to the first block in the buffer
		}
		const BlockType last() const {
			return buff_.col(i_ ? i_ - 1 : size() - 1); // i_-1 points to the last
		}
		const BlockType	at(int i) const
		{
			int j = (i_ + i) % size();
			return buff_.col(j); 
		}

		void reset()
		{
			buff_.setZero();
            i_ = 0;
		}

	private:
		MatrixType buff_;
		int i_;

		template <typename U>
		void push_scalar(const U &v)
		{
			buff_(0, i_) = v;
			i_++;
            i_ %= size();
		}
	};

	//----------------------------------------------------------------------------//

	/**
	 * @brief Block RLS estimator 
	 * 
	 * RLS estimation is employed on a fixed-size, rolling block of data
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T>
	class BlockRLS : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		RingBuffer<T> phi_buff_, data_buff_; // window storage buffers
		Eigen::LLT< MatrixType > llt;
		MatrixType P;		// Covariance Matrix
		VectorType k;		// Gain Vector
		VectorType temp;	// temporary buffer
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

	public:
		explicit BlockRLS(int n, int w = 10, T init_covar = 1.e6) 
                : RlsEstimatorBase<T>(n), phi_buff_(w, n), data_buff_(1, n),
                  P(n,n), k(n), temp(n)
		{
			reset(init_covar);
		}
		void update(const VectorType& phi, const T& data)
		{
			// UPDATE //
			temp = llt.matrixU() * phi;
			T d = 1./(1. + temp.squaredNorm());
			temp = llt.matrixL() * temp;
			k = temp * d; // Gain Vector
			T error = data - theta_.dot(phi);	// Error Value
			theta_ += k * error;			

			llt.rankUpdate(temp, -d);

			cost_ += error*error*d*(1-d);
			error = data - theta_.dot(phi);
			cost_ += error*error;

			// DOWNDATE //
            temp = llt.matrixU() * phi_buff_.first();
			d = temp.squaredNorm();
			if (d != 1.) {
				temp = llt.matrixL() * temp;
				d = 1./(1. - d);  
				k = temp * d; 
				error = data_buff_.first()(0) - theta_.dot(phi_buff_.first()); 
				theta_ -= k * error; 

				llt.rankUpdate(temp, d);
				
				cost_ -= error*error*d*(1-d);
				error = data_buff_.first()(0) - theta_.dot(phi_buff_.first());
				cost_ -= error*error;
			}

			// store phi, data in window buffers
			phi_buff_.push(phi);
			data_buff_.push(data);

			iter_++; // Update number of iterations
		};

		int size() const { return phi_buff_.size(); }
		void setSize(int w) {
			phi_buff_.resize(w, np_);
			data_buff_.resize(w, 1);
			VectorType th = theta_;
			T c = cost_;
			MatrixType P0 = P;
			reset();
			theta_ = th;
			P = P0;
			cost_ = c;
		}

		// Reset Function
		void reset(T init_covar = 1.e6) noexcept
		{
			iter_ = 0;
			theta_.setZero();
			P.setIdentity();
			P *= init_covar;
			k.setZero();
			cost_ = 0;	
			llt.compute(P);		

			phi_buff_.reset();
			data_buff_.reset();

			// initialize regressors according to Zhang 2004
			VectorType p(np_);
			p.setZero();
			int i = 0;
			for (; i < size() - np_; ++i) phi_buff_.push(p);
			for (i = 0; i < np_; ++i)
			{
				p.setZero();
				p(i) = 1 / std::sqrt(init_covar);
				phi_buff_.push(p);
			}
		};
	};

	//----------------------------------------------------------------------------//

	/**
	 * @brief Block RLS estimator with Cholesky decomposition
	 * 
	 * Block RLS utilizing up-/down-dating of the Cholesky decomposition of the
	 * normal equation matrix.
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T>
	class CholeskyBlockRls : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		RingBuffer<T> phi_buff_; // block storage buffer
		VectorType aug_phi;
		Eigen::LLT< MatrixType > llt;
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

	public:
		explicit CholeskyBlockRls(int n, int w = 10)
		: RlsEstimatorBase<T>(n), phi_buff_(w, n+1), aug_phi(n+1)
		{
			reset();
		}

		void update(const VectorType& phi, const T& y)
		{
			aug_phi.head(np_) = phi;
			aug_phi(np_) = y;

        	if (iter_ >= size()) { // buffer is full
				// update
				llt.rankUpdate(aug_phi);
				// downdate
				llt.rankUpdate(phi_buff_.first(), -1);
				solve();
				phi_buff_.push(aug_phi);
			} else {
				phi_buff_.push(aug_phi);
				if (iter_ < np_) {} // very few data 
				else if (iter_ == np_) {
					// init LLT
					const MatrixType& b = phi_buff_.buff();
					llt.compute((b.leftCols(np_+1)) * (b.leftCols(np_+1).adjoint()));
					solve();
				} else {
					// update only
					llt.rankUpdate(aug_phi);
					solve();
				}
			}

			iter_++; // Update number of iterations
		};

		int size() const { return phi_buff_.size(); }
		void setSize(int w) {
			phi_buff_.resize(w, np_+1);
			VectorType th = theta_;
			T c = cost_;
			reset();
			theta_ = th;
			cost_ = c;
		}


		// Reset Function
		void reset() noexcept
		{
			RlsEstimatorBase<T>::reset();
			phi_buff_.reset();
		};

	private:
		void solve()
		{
			const MatrixType& L =  llt.matrixLLT();
			cost_ = L(np_, np_)*L(np_, np_);
			theta_ = L.bottomLeftCorner(1, np_).adjoint();
			L.adjoint().topLeftCorner(np_, np_).template triangularView<Eigen::Upper>().template solveInPlace(theta_);
		}
	};

	//-----------------------------------------------------------------//

	// Subclass of RLS_Estimator that implemenets a polynomial estimation

	/**
	 * @brief Polynomial RLS estimator object 
	 * 
	 * A data sequence $y_t$, $t=0,1,2,...$ is fitted to a polynomial
	 * $y(t) = \theta_0 + \theta_1\cdot t + theta_2\cdot t^2 + ...$
	 * using one of the RLS estimators
	 * 
	 * @tparam EstimatorType_ The type of estimator used
	 */
	template < typename EstimatorType_ >
	class PolyRLS : public EstimatorType_
	{
	public:
		typedef EstimatorType_ EstimatorType;
		using VectorType = typename EstimatorType::VectorType;
		using ScalarType = typename EstimatorType::ScalarType;
	protected:
		using EstimatorType::iter_;
		using EstimatorType::np_;
		using EstimatorType::theta_;
		VectorType phi_;
	public:
		explicit PolyRLS(int n) : EstimatorType(n), phi_(n)
		{
			for (int i = 0; i < n; i++)	phi_(i) = 1;
		}

		void update(const ScalarType &data)
		{
			// Set regressors with regards to updates//
			//  phi(0) = 1.; this is set in constructor
			for (int i = 1; i < np_; i++)
				phi_(i) = phi_(i - 1) * iter_;

			EstimatorType::update(phi_, data);
		}
		ScalarType estimatedOutput() const noexcept
		{
			return EstimatorType::estimatedOutput(phi_);
		}
		ScalarType estimatedRate() const noexcept
		{
			ScalarType r = 0.;
			for (int i = 1; i < np_; i++)
					r += i * phi_(i - 1) * theta_(i);
			return r;
		}
		VectorType phi() const { return phi_; }
	};

} // namespace RLS

#endif
