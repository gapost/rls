# Command line filtering utilities

## 1. ``rlsfilter``

This is a command line filtering utility that reads values of $[\phi(t),y(t)]$ from ``stdin`` and uses RLS to obtain the parameters of the model $y(t;\theta) = \theta^T \cdot \phi(t)$.

A forgetting factor (ff) or sliding-block algorithm can be selected.

The output is written to ``stdout`` and contains the current signal prediction $\hat{y}(t)$, the value of the cost function (sum of squares) and the values of the estimated parameters $\theta_1, \theta_2, ...$. Numbers are separated by tabs.

### Usage

Typing ``rlsfilter -h`` gives the following help message:

```
Usage:
  rlsfilter [options] 

Recursive least squares (RLS) fitting of a linear model 
y(t|θ)=θ^Τ * φ(t) to a time series using a forgetting factor or
sliding-block algorithm.

Options:
  -nX  : number of fitting parameters, X=1,2,3,... (default = 2) 
  -ffX : forgetting factor, 0<X<=1, (default = 0.98) 
  -wX  : block width, X=0 means no block (default)
  -sqrt : update square root of covariance
  -h   : display this help message

rlsfilter reads vector φ(t) and the signal y(t), t=1,2,..., from stdin
as a stream of ascii-coded real values until eof or a non-numeric input
is encountered.
The output is written to stdout. Each output line contains the following tab
separated data: estimated signal, sum of squared residuals,
estimated parameters

Use redirection for file input/output, e.g.:
  polyrlsfilter < in.dat > out.dat
```
## 2. ``polyrlsfilter``

This is a command line filtering utility that reads a signal $y(t)$ from ``stdin`` and uses the RLS algorithm to fit a polynomial $f(t) = a_0 + a_1\, t + a_2\, t^2 + ...$ to the signal. 

A forgetting factor (ff) or sliding-block RLS algorithm can be selected.

### Usage

Typing ``polyrlsfilter -h`` gives the following help message:

```
Usage:
  polyrlsfilter [options] 

Recursive least squares (RLS) fitting of a polynomial 
a_0 + a_1*t + ... to a time series using a forgetting factor or
sliding-block algorithm.

Options:
  -nX  : number of fitting parameters, X=2,3,... (default = 2) 
  -ffX : forgetting factor, 0<X<=1, (default = 0.98) 
  -wX  : block width, X=0 means no block (default)
  -sqrt : update square root of covariance
  -h   : display this help message

polyrlsfilter takes the input signal, y_t, t=1,2,..., from stdin
as a stream of ascii-coded real values until eof or a non-numeric input
is encountered.
The output is written to stdout. Each output line contains the following tab
separated data: estimated signal, estimated rate, sum of squared residuals,
estimated parameters

Use redirection for file input/output, e.g.:
  polyrlsfilter < in.dat > out.dat
```
## Building

The utilities can be built using ``cmake``. See the main [README](../README.md) for details. 

## Testing

The MATLAB/OVTAVE scripts ``rlstest.m`` and ``polyrlstest.m`` generate test data, save it to a file and use the cli utilities (rlsfilter and polyrlsfilter, respectively) to filter the data with RLS, load the output and generate some plots. The filter options can be changed in the script.

It is assumed the the scripts are run from within the ``filter/`` directory and the executables reside in ``.build/filter/``.