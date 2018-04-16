Fast Evaluation of Bessel Functions
=====================================

An algorithm for evaluating the Bessel functions J_nu and Y_nu of the first and seconds kinds of 
nonnegative real orders and positive real arguments.  It also calculates either the values of a 
nonoscillatory phase function and its derivative (when in the oscillatory regime) or the values of 
the logarithms of the functions J_nu and -Y_nu  (when in the nonoscillatory regime).  For the
most part, it operates via a table of numerically precomputed expansions, although asymptoptic 
expansions  are used in the cases of extremely large and extremely small  arguments and series 
expansions are used to evaluate Bessel functions of small orders and arguments.  On a typical 
laptop machine available circa  2017 the average time required for a single evaluation was 
approximately  5 x 10^(-7) seconds.

The algorithm is described in some detail in the preprint

    James Bremer, "An algorithm for the numerical evaluation of Bessel function real orders and
     arguments."  arXiv:1705.07820.

The file bessel_eval.f90 contains an implementation of the algorithm.  It uses
precomputed expansions contained in the files bessel_data.bin and bessel_data16.bin.
The former contains versions of the precomputed table which is accurate to
roughly 16 digits of precision while the latter contains a version of the 
precomputed table which achieves roughly 10^(-25) digits of accuracy (when computations 
are perfomed using extended precision arithemtic, of course).

To compile and test on a linux/GNU system:

`make`
