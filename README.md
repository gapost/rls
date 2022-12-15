# RLS

Recursive least squares in C++

## About

Recursive least squares (RLS) refers to algorithms that recursively find the coefficients that minimize a weighted linear least squares cost function. 

It is assumed that a measured time-dependent quantity $y(t)$ follows a linear model plus noise

$$ y(t) = \theta^T \cdot \phi(t) + e(t) $$

where $\theta$ is a parameter vector, $\phi(t)$ is a vector of known functions of time and $e(t)$ is white noise, $\langle e(t) \rangle = 0$, $\langle e(t)e(t') \rangle = \sigma^2\delta(t-t')$.

A least squares cost function is defined as

$$ J(\theta) = \sum_t {w(t) \left[ y(t) - \theta^T \cdot \phi(t) \right]^2} $$

where $w(t)$ is a weighting function.

RLS algorithms seek to estimate $\theta$ that minimizes $J$ and recursively refine this estimation as new data become available. This is acheived with a very elegant numerical technique which is given here very briefly. 

The solution of the least squares problem (for $w(t)=1$) is given by the normal equations

$$ \left[ \Phi^T \cdot \Phi \right] \theta = A\cdot \theta = \Phi^T\cdot Y \rightarrow \theta = A^{-1}\cdot \Phi^T\cdot Y$$

where $\Phi = [\phi(1) \; \phi(2) \dots ]^T$, $Y = [y(1) \; y(2) \dots ]^T$

When a new data point becomes available the matrix $A$ is modified as

$$ A' = A + \phi \cdot \phi^T $$
 
 The inverse of $A'$, $P'=A'^{-1}$  can be obtained from $P=A^{-1}$ using the **Wilkinson** matrix formula:

 $$ P' = (A + \phi \cdot \phi^T)^{-1} = P - \frac{P\cdot\phi\cdot\phi^T\cdot P}{1 + \phi^T\cdot P\cdot \phi} $$

 Thus, in RLS algorithms we keep a copy of $P$ which is updated by the above formula. Note that $P$ is the covariance matrix of $\theta$.

 The parameters are updated by 

$$ e(t) = y(t) - \theta^T(t-1)\cdot \phi(t) $$

$$ \theta' = \theta + k\, e(t) $$

where $k$ is the "gain" which is defined below.

## Implemented Algorithms

### 1. Exponentially weighted RLS

In this algorithm the weighting function is

$$ w(t) = \lambda ^ {n-t}, \quad t\leq n$$

where $0 < \lambda \leq 1$ is the "forgetting factor". The data are exponentially weighted so that older data are "forgotten". The effective   forgetting time constant $\tau$ can be found from

$$ \lambda = e^{-1/\tau} \approx 1 - 1/\tau $$

The recursion is given by

> $$ e(t) = y(t) - \theta^T(t-1)\cdot \phi(t) $$
> $$ u = P(t-1)\cdot \phi(t) $$
> $$ k = u / \left[ \lambda + \phi^T(t)\cdot u\right] $$
> $$ \theta(t) = \theta(t-1) + k\, e(t) $$
> $$ P(t) = \left[P(t-1) - k \cdot u^T\right]/\lambda $$
> $$ J(t) = \lambda J(t-1) + e^2(t) $$

The algorithm is still valid for $\lambda=1$ but with infinite memory length ($\tau \to \infty$)

### 2. Exponentially weighted RLS with square-root algorithm

The square root algorithms (initially due to Potter (1963)) utilize the fact that $A$ and $P$ are symmetric, positive definite matrices. This can be understood e.g. by the fact that $J(\theta)$ close to the minimum can be expanded to 2nd order in $\theta$
$$ J = J_0 + \theta^{-1} P^{-1} \theta $$
thus the matrix $A=P^{-1}$ must be positive definite so that there is a minimum.

Such matrices can be Cholesky-decomposed as $P=Q\cdot Q^T$ where $Q$ is lower triangular.
$Q$ is also called the square root of $P$, thus the name of the algorithm.
The covariance update can be turned into an update of $Q$:

$$ P' = Q'\cdot Q'^T = Q\left[ I - \alpha u\cdot u^T \right] Q^T$$

where $u=Q^T \cdot\phi$ and $\alpha = (1 + u^T\cdot u )^{-1}$.

It can be easily shown that

$$ \left[ I - \alpha u\cdot u^T \right]  = \left[ I - \sigma u\cdot u^T \right]^2 $$

with $\sigma = \alpha / \left[ 1 + \sqrt{1 + \alpha u^T\cdot u}\right]$. Finally

$$ Q' = Q \left[ I - \sigma u\cdot u^T \right]$$

The benefit of $Q$ updating is that it ensures that $P$ retains symmetry and positive-definetness.

The RLS with recursion for $Q$ becomes

> $$ e(t) = y(t) - \theta^T(t-1)\cdot \phi(t) $$
 
> $$ u = Q^T(t-1) \cdot \phi(t) $$

> $$ \beta = \lambda + u^T \cdot u $$

> $$ \alpha = 1 / \left[\beta + \sqrt{\beta\,\lambda}\right] $$

> $$ k = Q(t-1) \cdot u $$

> $$ \theta(t) = \theta(t-1) + k\, [e(t)/\beta] $$

> $$ Q(t) = \left[Q(t-1) - \alpha \, k\cdot u^T\right]/\sqrt{\lambda} $$

> $$ J(t) = \lambda J(t-1) + e^2(t) $$

This algorithm is taken from Ljung & Soederstroem (1987) "Theory & Practice of Recursive Identification", p. 328

Note that this algorithm does not make $Q$ lower triangular. However it ensures that $P=Q\cdot Q^T$.

**TODO:** We are updating only the lower triangular Q and it is working. Why??

### 3. Block RLS

The rectangular moving block RLS uses the last $N$ points to estimate the parameters.
In each recursion the estimate is first *updated* with the new data point and then *downdated* by removing the oldest point.

The updating sequence is given by

> $$ e(t) = y(t) - \theta^T(t-1)\cdot \phi(t) $$
> $$ u = P(t-1)\cdot \phi(t) $$
> $$ \beta = \left[ 1 + \phi(t)^T \cdot u \right]^{-1} $$
> $$ k = \beta \, u $$
> $$ \bar{\theta}(t) = \theta(t-1) + k\, e(t) $$
> $$ \bar{P}(t) = \left[P(t-1) - k \cdot u^T\right] $$
> $$ e'(t) = y(t) - \theta(t)^T\cdot \phi(t) $$
> $$ \bar{J}(t) = J(t-1) + e^2(t)\,\beta\,(1-\beta) + e'^2(t) $$

The down-dating sequence is

> $$ e(t-N) = y(t-N) - \bar{\theta}^T(t)\cdot \phi(t-N) $$
> $$ u = P(t)\cdot \phi(t-N) $$
> $$ \beta = \left[ 1 - \phi(t-N)^T \cdot u \right]^{-1} $$
> $$ k = \beta\, u $$
> $$ \theta(t) = \theta(t) - k\, e(t-N) $$
> $$ P(t) = \left[P(t-1) + k \cdot u^T\right] $$
> $$ e'(t-N) = y(t) - \theta(t)^T\cdot \phi(t-N) $$
> $$ J(t) = \bar{J}(t) - e^2(t-N)\,\beta\,(1-\beta) - e'^2(t-N) $$

### 3. Block RLS with square root update

Similar to the above, the block RLS algorithm can be formulated with update of the square root of $P$.



