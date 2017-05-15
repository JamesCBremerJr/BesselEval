
program bessel_test3

use besseleval

implicit double precision (a-h,o-z)


pi = acos(-1.0d0)
call bessel_eval_init(dize)


print *,"Enter dnu"
read *,dnu

print *,"Enter t"
read *,t

call bessel_eval(dnu,t,alpha,alphader,vallogj,vallogy,valj,valy)

print *,""
print *,"alpha     = ",alpha
print *,"alphader  = ",alphader
print *,"vallogj   = ",vallogj
print *,"vallogy   = ",vallogy
print *,"valj      = ",valj
print *,"valy      = ",valy
print *,""


! if (alpha .ne. 0) then
! valj = sqrt(2/(pi*t)) * cos(alpha)/sqrt(alphader)
! valy = sqrt(2/(pi*t)) * sin(alpha)/sqrt(alphader)
! print *,valj
! print *,valy
! endif


end program
