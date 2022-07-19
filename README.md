# RLS

Recursive least squares in C++

# About

Recursive least squares (RLS) is an adaptive filter algorithm that recursively finds the coefficients that minimize a weighted linear least squares cost function relating to the input signals. This approach is in contrast to other algorithms such as the least mean squares (LMS) that aim to reduce the mean square error. In the derivation of the RLS, the input signals are considered deterministic, while for the LMS and similar algorithm they are considered stochastic. Compared to most of its competitors, the RLS exhibits extremely fast convergence. However, this benefit comes at the cost of high computational complexity.

The algorithm works as follows:
$`N= number of parameters`$, 
$`λ= forgetting factor`$, 
$`Φ= \begin{bmatrix}a1(n)\\a2(n)\\...\end{bmatrix}`$
, $`Φ(n)`$ is the matrix of regressors used to calculate the parameters.

$`P= \begin{bmatrix} P(0) & 0 & ..\\0 & P(0) & ....\\...\end{bmatrix}`$
,$`P`$ is the covariance matrix and $`P(0)`$ is an initial large value to declare indifference.

$`e(n)`$ is the error from new data with the previous parameters.

$`d(n)`$ is the new data.

$`g(n)`$ is the gain vector used for correction.

$`w(n)`$ is the vector of parameters that are calculated.

## Recursion:

$`e(n) = d(n) -Φ^T(n)w(n-1)`$ , calculating new error.

$`g(n) = P(n-1)Φ(n) / ( λ + Φ^T(n)P(n-1)Φ(n) )`$, calculating new gain vector.

$`P(n) = (λ^-1)(P(n-1) -g(n)Φ^T(n)P(n-1))`$, calculate new covariance matrix.

$`w(n) = w(n-1) + e(n)g(n)`$, calculate new parameters, estimate and repeat.





