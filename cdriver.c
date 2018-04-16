/* Test bessel_eval calling from C, double prec. Barnett 4/16/18

Usage: (also see makefile)

gfortran cdriver.c bessel_eval.f90 -o cdriver -O3 -march=native
time ./cdriver

Note: -march=native flag doubles the speed on i7 with GCC 5.4.

Todo: fix calling names to the f90 module, seem too baroque and nonportable.
     Link against the .mod ?
*/

#include <time.h>
#include <stdio.h>
#include <math.h>

extern void __besseleval_MOD_bessel_eval_init(double *dsize);
extern void __besseleval_MOD_bessel_eval(double *nu,double *t,double *alpha,double *alphader,double *vallogj,double *vallogy,double *valj,double *valy);

int main()
{
  int i;
  double nu, t, dt, dsize;
  double alpha,alphader,vallogj,vallogy,valj,valy;  // for outputs
  
  clock_t start=clock();
  __besseleval_MOD_bessel_eval_init(&dsize);
  printf("dsize=%g MB, init time=%g sec\n",dsize, (clock()-start)/(double)CLOCKS_PER_SEC);

  nu = 7.0;
  t = 200;
  __besseleval_MOD_bessel_eval(&nu,&t,&alpha,&alphader,&vallogj,&vallogy,&valj,&valy);
  printf("J_%g(%g) = %.16g,\terror = %.3g\n",nu,t,valj,fabs(valj-0.055762660213175076655));

  nu = 100.0;    // >=300 gives "!" output for some values, guess underflow?
  start = clock();
  printf("making 1e7 calls at nu=%g...\n",nu);
  dt = 1e-3;
  for (i=1;i<=10000000;++i) {  // i=0 crashes it!
    t = dt*i;
    __besseleval_MOD_bessel_eval(&nu,&t,&alpha,&alphader,&vallogj,&vallogy,&valj,&valy);
  }
  printf("done in %g sec.\n",(clock()-start)/(double)CLOCKS_PER_SEC);
  return 0;
}

/* timing comparison on i7-7700HQ:

MATLAB R2016b:
>> x=1e-3*(1:1e7); tic; besselj(100,x); toc          % uses all 8 threads.
Elapsed time is 5.137431 seconds.

cf 2.6 sec for bessel_eval, 1 thread.
*/
