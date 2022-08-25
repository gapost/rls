# rlsfilter

This folder contains the source and cmake build file for the `rlsfilter` program.

It is a simple command line utility that filters a signal $`y_i`$ entered from ``stdin`` and writes the filter output to ``stdout``.

The filter uses the RLS algorithm to fit a polynomial $`y(t) = a_0 + a_1\, t + a_2\, t^2 + ...`$ to the signal. A forgetting factor (ff) or block RLS algorithm can be selected by the user.

Each line of output from `rlsfilter` contains the current signal prediction $`y(t_i)`$ and the values of estimated parameters $`a_0, a_1, ...$`. The numbers are separated by tabs.

## Usage

Typing ``rlsfilter -h`` gives the following help message:

    Usage:
      rlstst [options] 

    Options:
      -nX  : fitting polynomial order, X=1,2,3,... (default = 2) 
      -ffX : forgetting factor, 0<X<=1, (default = 0.98) 
      -wX  : block width, X=0 means no block (default)

    rlstst takes input from stdin and writes to stdout. Use redirection for
    file input/output, e.g.:
      rlstst < in.dat > out.dat

## Building

Use ``cmake`` to build the program. 

This can be done, e.g. on a Linux system, by opening a command prompt at the ``filter` directory and giving the following commands:

    > mkdir .build
    > cd .build
    > cmake ..
    > make

This will build the ``rlsfilter`` executable in the ``.build`` sub-directory.

## Testing

The MATLAB/OVTAVE script ``filtertst.m`` generates a test signal, passes it through ``rlsfilter`` and plots the output. The filter options can be changed in the script.

It is assumed the the script is run from within the ``filter/`` directory and the executable resides at ``filter/.build/rlsfilter``.