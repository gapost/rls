# RLS

Recursive least squares in C++

# About

Recursive least squares (RLS) is an adaptive filter algorithm that recursively finds the coefficients that minimize a weighted linear least squares cost function relating to the input signals. This approach is in contrast to other algorithms such as the least mean squares (LMS) that aim to reduce the mean square error. In the derivation of the RLS, the input signals are considered deterministic, while for the LMS and similar algorithm they are considered stochastic. Compared to most of its competitors, the RLS exhibits extremely fast convergence. However, this benefit comes at the cost of high computational complexity.

## Exponential:
The exponential algorithm works as follows:
$`N= number  of  parameters`$, 

$`λ= forgetting  factor`$, 

$`Φ= \begin{bmatrix}r_1(n)\\r_2(n)\\...\end{bmatrix}`$ , is the matrix of regressors used to calculate the parameters.

$`P= \begin{bmatrix} P(0) & 0 & ..\\0 & P(0) & ....\\...\end{bmatrix}`$ , is the covariance matrix and $`P(0)`$ is an initial large value to declare indifference.

$`e(n)`$ is the error from new data with the previous parameters.

$`d(n)`$ is the new data.

$`g(n)`$ is the gain vector used for correction.

$`w(n)`$ is the vector of parameters that are calculated.

## Recursion:

$`e(n) = d(n) - Φ^T(n)w(n-1)`$ , calculating new error.

$`g(n) = P(n-1)Φ(n) / ( λ + Φ^T(n)P(n-1)Φ(n) )`$, calculating new gain vector.

$`P(n) = (1/λ)(P(n-1) -g(n)Φ^T(n)P(n-1))`$, calculate new covariance matrix.

$`w(n) = w(n-1) + e(n)g(n)`$, calculate new parameters, estimate and repeat.

The algorithm typically uses values for the forgetting factor between 0.98 and 1. Its 
purpose is to accurately estimate the next value of the output depending on the input
of the regressors by calculating the correct parameters. Changing the forgetting factor
changes how much the algorithm takes into account past values of the output.

## Rectangular Window:

The rectangular window approach uses the N last points to estimated the future values.
It uses the same algorithm as the exponential, with a forgetting factor of 1, in addition to removing 
observations from the estimated parameters.

## Recursion:

$`e_(i+N,i) = d_(i+N) - (Φ_(i+N))^T*w_(i+N-1,i)`$ , calculating new error.

$`g_(i+N,i) = P_(i+N-1,i)Φ_(i+N) / ( 1 + (Φ_(i+N))^T*P_(i+N-1,i)Φ_(i+N) )`$, calculating new gain vector.

$`P_(i+N,i) = (P_(i+N-1,i) - g_(i+N,i)*(Φ_(i+N))^T*P_(i+N-1,i))`$, calculate new covariance matrix.

$`w_(i+N,i) = w_(i+N-1,i) + e_(i+N,i)*g_(i+N,i)`$, calculate new parameters, estimate and repeat.


In addition to the above recursion, the following recursive operation takes place when the algorithm
has been give over N points, where N the length of the rectangular window.

$`e_(i+N,i) = d_i - (Φ_i)^T*w_(i+N,i)`$ , calculating new error.

$`g_(i+N,i+1) = P_(i+N,i)Φ_i / ( 1 - (Φ_i)^T*P_(i+N,i)Φ_i )`$, calculating new gain vector.

$`P_(i+N,i+1) = (P_(i+N,i) - g(n)(Φ_i)^T*P_(i+N,i))`$, calculate new covariance matrix.

$`w_(i+N,i+1) = w_(i+N,i) - e_(i+N,i)*g_(i+N,i+1) `$, calculate new parameters, estimate and repeat.


1) A general exponential estimator where the input are the regressors and the data
2) An exponential polynomial estimator where the regressors are calculated from the iterations 
of the class
3) A rectangular window estimator where only the N points are taken into account



