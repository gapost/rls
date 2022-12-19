#ifndef RLS_H_
#define RLS_H_

#include <vector>
#include <cmath>
#include <cstring>

namespace RLS
{
	enum UpdateType : int { CovarianceUpdate=1, SquareRootUpdate=2 };

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
		typedef typename std::vector<T> MatrixType;
		typedef typename std::vector<T> VectorType;

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
		T estimatedOutput(const VectorType &phi) const { 
			T v = 0;
			for(int i=0; i<np_; ++i) v += theta_[i]*phi[i];
			return v; 
		}
		void reset() {
			iter_ = 0;
			cost_ = 0;
			std::fill(theta_.begin(), theta_.end(), T(0) );
		}
		// UpdateType updateType() const { return (UpdateType)Upd_; }
	};

	template<typename T, class Estimator_>
	struct rls_upd_impl;

	/**
	 * @brief Exp weighted RLS estimator with forgetting factor
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T, int Upd_ = CovarianceUpdate>
	class ExpWeightedRLS : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		T ff_;				// Forgetting Factor
		T invsqrtff_;		// 1/sqrt(ff)
		MatrixType P;		// Covariance Matrix
		VectorType k;		// Gain Vector
		VectorType u;	    // temp Vector
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

		typedef ExpWeightedRLS<T, Upd_> MyT_;
		friend struct rls_upd_impl<T, MyT_>;

	public:															
		explicit ExpWeightedRLS(int n, T ff = 0.98, T init_covar = 1) 
		: RlsEstimatorBase<T>(n), ff_(ff), P(n*n), k(n), u(n)
		{
			reset(init_covar);
			setff(ff);
		}

		// Update of Parameters with New data (data)
		template<typename Vector_>
		void update(const Vector_& phi, const T& data)
		{
			rls_upd_impl<T, MyT_>::update(*this, &phi[0], data);
			// Update number of iterations
			iter_++; 
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
		const MatrixType &covar() const noexcept { return P; }
		const VectorType &gain()  const noexcept { return k; }

		// Reset Function
		void reset(T init_covar = 1) noexcept
		{
			iter_ = 0;
			cost_ = 0;
			std::fill(theta_.begin(), theta_.end(), T(0) );
			std::fill(k.begin(), k.end(), T(0) );
			std::fill(P.begin(), P.end(), T(0) );
			for(int i=0; i<np_; ++i) 
				P[i*(np_+1)] = Upd_==CovarianceUpdate ? init_covar : std::sqrt(init_covar); // P(i,i)=1			
		};
	};

	template<typename T>
	struct rls_upd_impl<T, ExpWeightedRLS<T, CovarianceUpdate> >
	{
		typedef ExpWeightedRLS<T, CovarianceUpdate> EstimatorType;
		
		static void update(EstimatorType& E, const T* phi, const T& data) {
			const T& ff_ = E.ff_;				// Forgetting Factor
			typename EstimatorType::MatrixType& P = E.P;		// Covariance Matrix
			typename EstimatorType::VectorType& k = E.k;		// Gain Vector
			typename EstimatorType::VectorType& u = E.u;	    // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			T b = ff_, error = data;
			for(int i=0; i<np_; ++i)
			{
				// u = P * phi
				T* p = P.data() + i*np_; // P(i,1)
				u[i] = *p * phi[0];
				for(int j=1; j<np_; ++j) {
					p++;
					u[i] += *p * phi[j];
				}
				// b = ff + phi*u
				b += u[i]*phi[i];
				// e = y - theta*phi
				error -= theta_[i]*phi[i];
			}
			for(int i=0; i<np_; ++i) {
				// Set gain vector
				// k = u / (phi.dot(u) + ff_);
				k[i] = u[i] / b;
				// Calculation of new parameters//
				theta_[i] += k[i]*error;
				// Calculation of new covariance//
				T* p = P.data() + i*(np_ + 1); // P(i,i)
				*p -= k[i] * u[i];
				*p /= ff_;
				for(int j=i-1; j>=0; --j) {
					p--;
					*p -= k[i] * u[j]; // P(i,j). j<i
					*p /= ff_;
					P[i+j*np_] = *p;  // P(j,i) = P(i,j)					
				}
			}
			// calc cost function
			cost_ = ff_ * cost_ + error * error;
		}
	};

	template<typename T>
	struct rls_upd_impl<T, ExpWeightedRLS<T, SquareRootUpdate> >
	{
		typedef ExpWeightedRLS<T, SquareRootUpdate> EstimatorType;
		
		static void update(EstimatorType& E, const T* phi, const T& data) {
			const T& ff_ = E.ff_;				// Forgetting Factor
			const T& invsqrtff_ = E.invsqrtff_; // 1/sqrt(ff)
			typename EstimatorType::MatrixType& Q = E.P; // square root Covariance Matrix
			typename EstimatorType::VectorType& k = E.k; // Gain Vector
			typename EstimatorType::VectorType& u = E.u; // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			// algorithm from Ljung & Soederstoem (1987) p. 328 
			T b = ff_, error = data;
			for (int i = 0; i < np_; ++i) {				
				// u = Q^T*phi
				T* q = Q.data() + i; // Q(0,i)
				u[i] = *q * phi[0];
				for (int j = 1; j < np_; ++j) {
					q += np_;
					u[i] += *q * phi[j]; // Q(j,i)*phi(j)
				}
				// b = ff + u^T*u
				b += u[i]*u[i]; 				
				// e = y - theta*phi
				error -= theta_[i]*phi[i];
			}
			// k = Q*u
			for (int i = 0; i < np_; ++i) {
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] = *q * u[0];
				for (int j = 1; j < np_; ++j) {
					q++;
					k[i] += *q * u[j]; 
				}
			}
			// calc cost function
			cost_ = ff_ * cost_ + error * error;			
			// Calculation of new covariance and theta//
			T a = 1/(b + std::sqrt(b*ff_));
			error /= b;
			for (int i = 0; i < np_; i++)
			{
				// Update parameters//
				theta_[i] += k[i]*error;
				// Update square root of covariance //
				// Q = (Q - a*k*u^T) / sqrt(ff)
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] *= a;
				*q -= k[i] * u[0];
				*q *= invsqrtff_;
				for (int j = 1; j < np_; ++j) {
					q++;
					*q -= k[i] * u[j]; 
					*q *= invsqrtff_;
				}
			}
		}
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
		/**
		 * @brief Construct a new Ring Buffer object
		 *
		 * @param n : number of vectors in buffer
		 * @param w : size of each vector
		 */
        RingBuffer(int n, int w) : sz_(n), w_(w), buff_(n*w)
		{
			reset();
		}

		int size() const { return sz_; }
		int width() const { return w_; }

		void resize(int n, int w) {
			sz_ = n;
			w_ = w;
            buff_.resize(n*w);
			reset();
		}

		const T* data() const { return buff_.data(); }

		void push(const T* v, int n) {
			std::memcpy(buff_.data() + i_*w_, v, std::min(width(), n)*sizeof(T));
			i_++;
            i_ %= size();
		}

		void push(const T& v) {
			buff_[i_*w_] = v;
			i_++;
            i_ %= size();
		}

		const T* first() const {
			return data() + i_*w_; // i_ points to the first element in the buffer
		}
		const T* last() const {
			return data() + (i_ ? i_ - 1 : size() - 1)*w_; // i_-1 points to the last
		}
		const T* at(int i) const
		{
			int j = (i_ + i) % size();
			return data() + j*w_; 
		}

		void reset()
		{
			std::fill(buff_.begin(), buff_.end(), T(0));
            i_ = 0;
		}

	private:
		int sz_;
		int w_;
		std::vector<T> buff_;
		int i_;		
	};

//----------------------------------------------------------------------------//

	/**
	 * @brief Block RLS estimator 
	 * 
	 * RLS estimation is employed on a fixed-size, rolling block of data
	 * 
	 * @tparam T Real number type (float, double, long double)
	 */
	template <typename T, int Upd_ = CovarianceUpdate>
	class BlockRLS : public RlsEstimatorBase<T>
	{
	public:
		using MatrixType = typename RlsEstimatorBase<T>::MatrixType;
		using VectorType = typename RlsEstimatorBase<T>::VectorType;

	protected:
		RingBuffer<T> phi_buff_, data_buff_; // window storage buffers
		MatrixType P;		// Covariance Matrix
		VectorType k;		// Gain Vector
		VectorType u;		// temporary buffer
		using RlsEstimatorBase<T>::theta_;
		using RlsEstimatorBase<T>::cost_;
		using RlsEstimatorBase<T>::np_;
		using RlsEstimatorBase<T>::iter_;

		int failedDowndate_;

		typedef BlockRLS<T, Upd_> MyT_;
		friend struct rls_upd_impl<T, MyT_>;

	public:
		explicit BlockRLS(int n, int w = 10, T init_covar = 1) 
                : RlsEstimatorBase<T>(n), phi_buff_(w, n), data_buff_(1, n),
                  P(n*n), k(n), u(n)
		{
			reset(init_covar);
			setSize(w);
		}
		template<typename Vector_>
		void update(const Vector_& phi, const T& data)
		{
			rls_upd_impl<T, MyT_>::update(*this, &phi[0], data);

			failedDowndate_ = 0;
			if (rls_upd_impl<T, MyT_>::downdate(*this, 
				phi_buff_.first(), *(data_buff_.first())) < 0) {
				// TODO
				// What happens if downdate fails?
				// Which means that the updated P is not positive-definite
				failedDowndate_ = 1;
			}

			// store phi, data in window buffers
			phi_buff_.push(phi.data(),np_);
			data_buff_.push(data);

			iter_++; // Update number of iterations
		};

		int failedDowndate() const { return failedDowndate_; }
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
		void reset(T init_covar = 1) noexcept
		{
			iter_ = 0;
			cost_ = 0;
			std::fill(theta_.begin(), theta_.end(), T(0) );
			std::fill(k.begin(), k.end(), T(0) );
			std::fill(P.begin(), P.end(), T(0) );
			for(int i=0; i<np_; ++i) 
				P[i*(np_+1)] = Upd_==CovarianceUpdate ? 
				               init_covar : std::sqrt(init_covar); 			

			phi_buff_.reset();
			data_buff_.reset();

			// initialize regressors according to Zhang 2004
			VectorType p(np_, T(0));
			int i = 0;
			for (; i < size() - np_; ++i) phi_buff_.push(p.data(),np_);
			for (i = 0; i < np_; ++i)
			{
				VectorType p1(np_, T(0));
				p1[i] = 1 / std::sqrt(init_covar);
				phi_buff_.push(p1.data(),np_);
			}
		};
	};

	template<typename T>
	struct rls_upd_impl<T, BlockRLS<T, CovarianceUpdate> >
	{
		typedef BlockRLS<T, CovarianceUpdate> EstimatorType;
		
		static void update(EstimatorType& E, const T* phi, const T& data) {
			typename EstimatorType::MatrixType& P = E.P;		// Covariance Matrix
			typename EstimatorType::VectorType& k = E.k;		// Gain Vector
			typename EstimatorType::VectorType& u = E.u;	    // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			T b(1), error(data), error2(data);
			for(int i=0; i<np_; ++i)
			{
				// u = P * phi
				T* p = P.data() + i*np_; // P(i,1)
				u[i] = *p * phi[0];
				for(int j=1; j<np_; ++j) {
					p++;
					u[i] += *p * phi[j];
				}
				// b = 1 + phi*u
				b += u[i]*phi[i];
				// e = y - theta*phi
				error -= theta_[i]*phi[i];
			}
			b = 1/b;
			cost_ += error*error*b*(1-b);
			for(int i=0; i<np_; ++i) {
				// Set gain vector
				// k = u / (phi.dot(u) + 1);
				k[i] = u[i] * b;
				// Calculation of new parameters//
				theta_[i] += k[i]*error;
				// e = y - theta*phi
				error2 -= theta_[i]*phi[i];
				// Calculation of new covariance//
				T* p = P.data() + i*(np_ + 1); // P(i,i)
				*p -= k[i] * u[i];
				for(int j=i-1; j>=0; --j) {
					p--;
					*p -= k[i] * u[j]; // P(i,j). j<i
					P[i+j*np_] = *p;  // P(j,i) = P(i,j)					
				}
			}
			cost_ += error2*error2;
		}

		static int downdate(EstimatorType& E, const T* phn, const T& data) {
			typename EstimatorType::MatrixType& P = E.P;		// Covariance Matrix
			typename EstimatorType::VectorType& k = E.k;		// Gain Vector
			typename EstimatorType::VectorType& u = E.u;	    // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			T b(1), error(data), error2(data);
			for(int i=0; i<np_; ++i)
			{
				// u = P * phi
				T* p = P.data() + i*np_; // P(i,1)
				u[i] = *p * phn[0];
				for(int j=1; j<np_; ++j) {
					p++;
					u[i] += *p * phn[j];
				}
				// b = 1 - phi*u
				b -= u[i]*phn[i];
				// e = y - theta*phi
				error -= theta_[i]*phn[i];
			}
			if (b==T(0)) return -1;
			b = 1/b;
			cost_ -= error*error*b*(1-b);
			for(int i=0; i<np_; ++i) {
				// Set gain vector
				// k = u / (phi.dot(u) + 1);
				k[i] = u[i] * b;
				// Calculation of new parameters//
				theta_[i] -= k[i]*error;
				// e = y - theta*phi
				error2 -= theta_[i]*phn[i];
				// Calculation of new covariance//
				T* p = P.data() + i*(np_ + 1); // P(i,i)
				*p += k[i] * u[i];
				for(int j=i-1; j>=0; --j) {
					p--;
					*p += k[i] * u[j]; // P(i,j). j<i
					P[i+j*np_] = *p;  // P(j,i) = P(i,j)					
				}
			}
			cost_ -= error2*error2;
			return 0;
		}
	};

	template<typename T>
	struct rls_upd_impl<T, BlockRLS<T, SquareRootUpdate> >
	{
		typedef BlockRLS<T, SquareRootUpdate> EstimatorType;
		
		static void update(EstimatorType& E, const T* phi, const T& data) {
			typename EstimatorType::MatrixType& Q = E.P;		// Square roor Covariance Matrix
			typename EstimatorType::VectorType& k = E.k;		// Gain Vector
			typename EstimatorType::VectorType& u = E.u;	    // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			T b(1), error(data), error2(data);
			for(int i=0; i<np_; ++i)
			{
				// u = Q^T*phi
				T* q = Q.data() + i; // Q(0,i)
				u[i] = *q * phi[0];
				for (int j = 1; j < np_; ++j) {
					q += np_;
					u[i] += *q * phi[j]; // Q(j,i)*phi(j)
				}
				// b = 1 + u^T*u
				b += u[i]*u[i];
				// e = y - theta*phi
				error -= theta_[i]*phi[i];
			}
			// k = Q*u
			for (int i = 0; i < np_; ++i) {
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] = *q * u[0];
				for (int j = 1; j < np_; ++j) {
					q++;
					k[i] += *q * u[j]; 
				}
			}
			b = 1/b;
			cost_ += error*error*b*(1-b);
			T a = b/(1 + std::sqrt(b));
			error *= b;
			for(int i=0; i<np_; ++i) {
				// Calculation of new parameters//
				theta_[i] += k[i]*error;
				// e = y - theta*phi
				error2 -= theta_[i]*phi[i];
				// Update square root of covariance //
				// Q = (Q - a*k*u^T)
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] *= a;
				*q -= k[i] * u[0];
				for (int j = 1; j < np_; ++j) {
					q++;
					*q -= k[i] * u[j]; 
				}
			}
			cost_ += error2*error2;
		}

		static int downdate(EstimatorType& E, const T* phi, const T& data) {
			typename EstimatorType::MatrixType& Q = E.P;		// Square roor Covariance Matrix
			typename EstimatorType::VectorType& k = E.k;		// Gain Vector
			typename EstimatorType::VectorType& u = E.u;	    // temp Vector
			typename EstimatorType::VectorType& theta_ = E.theta_;
			T& cost_ = E.cost_;
			const int& np_ = E.np_;

			T b(1), error(data), error2(data);
			for(int i=0; i<np_; ++i)
			{
				// u = Q^T*phi
				T* q = Q.data() + i; // Q(0,i)
				u[i] = *q * phi[0];
				for (int j = 1; j < np_; ++j) {
					q += np_;
					u[i] += *q * phi[j]; // Q(j,i)*phi(j)
				}
				// b = 1 - u^T*u
				b -= u[i]*u[i];
				// e = y - theta*phi
				error -= theta_[i]*phi[i];
			}
			if (b <= T(0)) return -1;
			// k = Q*u
			for (int i = 0; i < np_; ++i) {
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] = *q * u[0];
				for (int j = 1; j < np_; ++j) {
					q++;
					k[i] += *q * u[j]; 
				}
			}
			b = 1/b;
			cost_ -= error*error*b*(1-b);
			T a = b/(1 + std::sqrt(b));
			error *= b;
			for(int i=0; i<np_; ++i) {
				// Calculation of new parameters//
				theta_[i] -= k[i]*error;
				// e = y - theta*phi
				error2 -= theta_[i]*phi[i];
				// Update square root of covariance //
				// Q = (Q - a*k*u^T)
				T* q = Q.data() + i*np_; // Q(i,1)
				k[i] *= a;
				*q += k[i] * u[0];
				for (int j = 1; j < np_; ++j) {
					q++;
					*q += k[i] * u[j]; 
				}
			}
			cost_ -= error2*error2;
			return 0;
		}
	};

//----------------------------------------------------------------------------//



} // namespace RLS

#endif
