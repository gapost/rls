#ifndef POLYRLS_H_
#define POLYRLS_H_

#include "RLS.h"

namespace RLS {

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
class poly_rls_impl : public EstimatorType_
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
    explicit poly_rls_impl(int n) : EstimatorType(n), phi_(n)
    {
        for (int i = 0; i < n; i++)	phi_[i] = 1;
    }

    void update(const ScalarType &data)
    {
        // Set regressors with regards to updates//
        // phi(0) = 1.; this is set in constructor
        for (int i = 1; i < np_; i++)
            phi_[i] = phi_[i-1] * iter_;

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
                r += i * phi_[i - 1] * theta_[i];
        return r;
    }
    VectorType phi() const { return phi_; }
    ScalarType t() const { return phi_[1]; }
};

template< typename T, int Upd_, typename EstimatorType_ >
class PolyRLS;

template< typename T, int Upd_ >
class PolyRLS<T, Upd_, ExpWeightedRLS<T, Upd_> > {    
public:
    typedef ExpWeightedRLS<T, Upd_> EstimatorType;		
    using VectorType = typename EstimatorType::VectorType;
    using ScalarType = typename EstimatorType::ScalarType;
private:
    typedef poly_rls_impl<EstimatorType> PolyRlsType;
    PolyRlsType P1, P2;
    PolyRlsType* p[2];   
    int W_;
public:
    explicit PolyRLS(int n) : P1(n), P2(n)
    {
        setff(0.98);
    }

    void update(const ScalarType &data)
    {
        if (p[0]->iter() == 2*W_) {
            // reset p[0] and swap
            p[0]->reset();
            PolyRlsType* d = p[0];
            p[0] = p[1];
            p[1] = d;
        }
        P1.update(data);
        P2.update(data);
    }
    ScalarType estimatedOutput() const noexcept
    {
        return p[0]->estimatedOutput();
    }
    ScalarType estimatedRate() const noexcept
    {
        return p[0]->estimatedRate();
    }
    VectorType phi() const { return p[0].phi(); }
    const VectorType& estimatedPar() const { return p[0]->estimatedPar(); }
    T cost() const { return p[0]->cost(); }
    void reset() {
        P1.reset(); P2.reset();
        for(int i=0; i<W_; ++i) P1.update(0);
        p[0] = &P1; p[1] = &P2;
    }
    T ff() const { return P1.ff(); }
    void setff(T ff) {
        if ((ff > 0) && (ff <= 1.0)) {
            P1.setff(ff); P2.setff(ff);
            W_ = (int) (1/(1-P1.ff()));
            W_ *= 20;
            reset();
        }
    }
    T t() const { return p[0]->t(); }
};

template< typename T, int Upd_ >
class PolyRLS<T, Upd_, BlockRLS<T, Upd_> > {
public:
    typedef BlockRLS<T, Upd_> EstimatorType;		
    using VectorType = typename EstimatorType::VectorType;
    using ScalarType = typename EstimatorType::ScalarType;
private:
    typedef poly_rls_impl<EstimatorType> PolyRlsType;
    PolyRlsType P1, P2;
    PolyRlsType* p[2];   
public:
    explicit PolyRLS(int n) : P1(n), P2(n)
    {
        reset();
    }

    void update(const ScalarType &data)
    {
        if (p[0]->iter() == 2*P1.size()) {
            // reset p[0] and swap
            p[0]->reset();
            PolyRlsType* d = p[0];
            p[0] = p[1];
            p[1] = d;
        }
        P1.update(data);
        P2.update(data);
    }
    ScalarType estimatedOutput() const noexcept
    {
        return p[0]->estimatedOutput();
    }
    ScalarType estimatedRate() const noexcept
    {
        return p[0]->estimatedRate();
    }
    VectorType phi() const { return p[0].phi(); }
    const VectorType& estimatedPar() const { return p[0]->estimatedPar(); }
    T cost() const { return p[0]->cost(); }
    void reset() {
        P1.reset(); P2.reset();
        for(int i=0; i<P1.size(); ++i) P1.update(0);
        p[0] = &P1; p[1] = &P2;
    }
    int size() const { return P1.size(); }
    void setSize(int w) {
        P1.setSize(w); P2.setSize(w);
        reset();
    }
    T t() const { return p[0]->t(); }
};

} // namespace RLS

#endif