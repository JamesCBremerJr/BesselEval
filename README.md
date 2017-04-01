Fast Evalaution of Bessel Functions
=====================================


An algorithm for evaluating the Bessel functions J_nu and Y_nu of the first and seconds kinds of nonnegative real orders and positive real arguments.  It also calculates either the values of a nonoscillatory phase function and its derivative (when in the oscillatory regime) or the values of the logarithms of the functions J_nu and -Y_nu  (when in the nonoscillatory regime).  It operates via numerically precomputed expansions rather than asymptoptic expansions and other devices obtained through analyses of Bessel's differential equation.  On a typical laptop machine available circa 2017 the average time required  for a  single evaluation was approximately  3 x 10^(-7) seconds.

The implementation of this algorithm consists of the following three files:

*  bessel_eval.f90
*  bessel_data.f90
*  bessel_test.f90

Data used by the bessel_test.f90 to perform testing of the code is contained in the text
files:

*  bessel1.txt
*  bessel2.txt
*  bessel3.txt
*  bessel4.txt

The main documentation is in the file bessel_eval.f90, which contains the only user-callable subroutine, bessel_eval.  The test code for this subroutine can be compiled and executed with the following
command:

    gfortran -O3 -o bessel_test.out bessel_data.f90 bessel_eval.f90 bessel_test.f90
    ./bessel_test.out

