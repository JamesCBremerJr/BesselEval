!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Copyright 2017 by James Bremer    
!
!  This program is free software: you can redistribute it and/or modify it under the terms 
!  of the GNU General Public License as published by the Free Software Foundation, either 
!  version 3 of the License, or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
!  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
!  See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  The code in this file tests the bessel_evalroutine through comparison
!  with highly accurate reference values which were computed in various ways
!  and are stored in the text files bessel?.txt.  
!
!  In addition to outputing results to the console, this program generates
!  LaTeX source code for some of the tables which appear in the companion paper.
!  This source code is written to text files named table?.tex.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine report_on_expansions()
use bessel_data
implicit double precision (a-h,o-z)

print *,"-------------------------------------------------------------------------------------"
print *,"Report on the characteristics of the expansions"
print *,"-------------------------------------------------------------------------------------"
print *,""

dnu1           = dnu11
dnu2           = dnu21
dmemory        = dmemory1
dtime          = dtime1

print "(A40,' = ',D10.3,' - ',D10.3)","Expansion dnu range",dnu1,dnu2
print "(A40,' = ',F8.2)","Memory used by expansion (in MB)",dmemory
print "(A40,' = ',F8.2)","Time to build expansion (in seconds)",dtime
print *,""

dnu1           = dnu12
dnu2           = dnu22
dmemory        = dmemory2
dtime          = dtime2

print "(A40,' = ',D10.3,' - ',D10.3)","Expansion dnu range",dnu1,dnu2
print "(A40,' = ',F8.2)","Memory used by expansion (in MB)",dmemory
print "(A40,' = ',F8.2)","Time to build expansion (in seconds)",dtime
print *,""


dnu1           = dnu13
dnu2           = dnu23
dmemory        = dmemory3
dtime          = dtime3

print "(A40,' = ',D10.3,' - ',D10.3)","Expansion dnu range",dnu1,dnu2
print "(A40,' = ',F8.2)","Memory used by expansion (in MB)",dmemory
print "(A40,' = ',F8.2)","Time to build expansion (in seconds)",dtime
print *,""

dnu1           = dnu14
dnu2           = dnu24
dmemory        = dmemory4
dtime          = dtime4

print "(A40,' = ',D10.3,' - ',D10.3)","Expansion dnu range",dnu1,dnu2
print "(A40,' = ',F8.2)","Memory used by expansion (in MB)",dmemory
print "(A40,' = ',F8.2)","Time to build expansion (in seconds)",dtime
print *,""

print "(A40,' = ',D10.3)","Total time = ",dtime1+dtime2+dtime3+dtime4
print "(A40,' = ',D10.3)","Total memory = ",dmemory1+dmemory2+dmemory3+dmemory4
print *,""

end subroutine



subroutine test_oscillatory(filename1,filename2)
implicit double precision (a-h,o-z)
character(len=*) :: filename1, filename2

double precision, allocatable :: dnus(:),ts(:),valsa(:),valsa0(:)
double complex, allocatable :: vals(:),vals0(:)
double complex              :: ima

ima = (0.0d0,1.0d0)
pi  = acos(-1.0d0)

print *,"-------------------------------------------------------------------------------------"
print *,"Testing in the oscillatory region ( ",filename1," )"
print *,"-------------------------------------------------------------------------------------"
print *,""

iw  = 20
iw2 = 21
open(iw, FILE=filename1)
open(iw2, FILE=filename2)

write(iw2,*) "\begin{tabular}{lcc}"
write(iw2,*) "\toprule"
write(iw2,*) "Range of $\nu$ & Maximum relative         & Average evaluation\\"
write(iw2,*) "               & error in $\alpha_\nu'$   & time (in seconds) \\"
write(iw2,*) "\midrule"


read(iw,*) nchunks

do ichunk=1,nchunks

read(iw,*) n
read(iw,*) dnu1
read(iw,*) dnu2

allocate(dnus(n),ts(n),valsa(n),valsa0(n),vals(n),vals0(n))

do i=1,n
read (iw,*) dnus(i)
read (iw,*) ts(i)
read (iw,*) valj
read (iw,*) valy
read (iw,*) valsa0(i)
vals0(i) = valj + ima*valy
end do

call elapsed(t1)
do i=1,n
dnu = dnus(i)
t   = ts(i)
call bessel_eval(dnu,t,alpha,alphader,beta1,beta2,valj,valy)
valsa(i) = alphader
vals(i)  = valj + ima * valy
end do

call elapsed(t2)

avg     = (t2-t1)/n
errrela = maxval( abs(valsa - valsa0) / abs(valsa0) )
errrel  = maxval( abs(vals - vals0)   / abs(vals0)  )

print "(A45,' = ',D25.16,' - ',D25.16)","Range of dnus",dnu1,dnu2
print "(A45,' =    ',I8.8)","Number of dnus",n
print "(A45,' = ',D25.16)","Average evaluation time",avg
print "(A45,' = ',D25.16)","Max relative error alpha'",errrela
print "(A45,' = ',D25.16)","Max relative error in J_nu + i Y_\nu ",errrel


call write_table_integer_range(iw2,floor(dnu1),floor(dnu2))
call write_table_next(iw2)
call write_table_double(iw2,errrela)
call write_table_next(iw2)
call write_table_double(iw2,avg)
call write_table_nextline(iw2)


print *,""
deallocate(dnus,ts,valsa,valsa0,vals,vals0)

end do


write (iw2,*) "\bottomrule"
write (iw2,*) "\end{tabular}"


close(iw)
close(iw2)

end subroutine


subroutine test_nonoscillatory(filename1,filename2)
implicit double precision (a-h,o-z)
character(len=*) :: filename1,filename2

double precision, allocatable :: dnus(:),ts(:)
double complex, allocatable   :: vals(:),vals0(:)
double precision, allocatable :: errs1(:),errs2(:),dnus1(:),dnus2(:)

double complex :: ima
ima = (0.0d0,1.0d0)

print *,"-------------------------------------------------------------------------------------"
print *,"Testing in the nonoscillatory regime ( ",filename1," )"
print *,"-------------------------------------------------------------------------------------"
print *,""

iw  = 20
open(iw, FILE=filename1)

iw2  = 21
open(iw2, FILE=filename2)


write(iw2,*) "\begin{tabular}{lcc}"
write(iw2,*) "\toprule"
write(iw2,*) "Range of $\nu$ & Maximum relative  & Average evaluation\\"
write(iw2,*) "               & error in          & time (in seconds) \\"
write(iw2,*) "               &$\log(J_\nu) +i \log(-Y_\nu)$\\"
write(iw2,*) "\midrule"

read(iw,*) nchunks

allocate(errs1(nchunks),errs2(nchunks),dnus1(nchunks),dnus2(nchunks))

do ichunk=1,nchunks

read(iw,*) n
read(iw,*) dnu1
read(iw,*) dnu2

dnus1(1) = dnu1
dnus2(1) = dnu2

allocate(dnus(n),ts(n),vals(n),vals0(n))

do i=1,n
read (iw,*) dnus(i)
read (iw,*) ts(i)
read (iw,*) vallogj0
read (iw,*) vallogy0

vals0(i) = vallogj0 + ima * vallogy0

end do


call elapsed(t1)
do i=1,n
dnu = dnus(i)
t   = ts(i)
call bessel_eval(dnu,t,alpha,alphader,vallogj,vallogy,valj,valy)
vals(i)     = vallogj + ima * vallogy
end do


call elapsed(t2)


avg    = (t2-t1)/n
err1   = maxval(abs(vals-vals0)/abs(vals0))

print "(A45,' = ',D25.16,' - ',D25.16)","Range of dnus",dnu1,dnu2
print "(A45,' =    ',I8.8)","Number of dnus",n
print "(A45,' = ',D25.16)","Average evaluation time",avg
print "(A45,' = ',D25.16)","Max error in  log(J_dnu) + i log(-Y_nu) ",err1


call write_table_integer_range(iw2,floor(dnu1),floor(dnu2))
call write_table_next(iw2)
call write_table_double(iw2,err1)
call write_table_next(iw2)
call write_table_double(iw2,avg)
call write_table_nextline(iw2)

print *,""

deallocate(dnus,ts,vals,vals0)
end do

write (iw2,*) "\bottomrule"
write (iw2,*) "\end{tabular}"

close(iw)
close(iw2)

end subroutine



subroutine test_vals(filename1,filename2)
implicit double precision (a-h,o-z)
character(len=*) :: filename1,filename2
double precision, allocatable :: ts(:)
double complex, allocatable   ::  vals(:),vals0(:)
double complex :: ima
ima = (0.0d0,1.0d0)

print *,"-------------------------------------------------------------------------------------"
print *,"Testing the values of J_nu and Y_nu ( ",filename1," )"
print *,"-------------------------------------------------------------------------------------"
print *,""

iw  = 20
open(iw,  FILE=filename1)



iw2  = 21
open(iw2, FILE=filename2)


write(iw2,*) "\begin{tabular}{lcc}"
write(iw2,*) "\toprule"


write(iw2,"(A)",advance="no") "& \multicolumn{2}{c}{\large \tt bessel\_eval} \\"
write(iw2,"(A)") "\addlinespace[.25em]"

write(iw2,"(A)",advance="no") "$\nu$ & Maximum relative & Average evaluation \\ "

write(iw2,"(A)",advance="no") "& error & time (in seconds) \\"


write(iw2,*) "\midrule"

read(iw,*) nchunks
read(iw,*) n

allocate(ts(n),vals0(n),vals(n))

do ichunk=1,nchunks

read(iw,*) nu
read(iw,*) dtime
dnu = nu


do i=1,n
read (iw,*) t
read (iw,*) valj
read (iw,*) valy

ts(i)    = t
vals0(i) = valj + ima * valy
end do

call elapsed(t1)
do i=1,n
t   = ts(i)
call bessel_eval(dnu,t,alpha,alphader,beta1,beta2,valj,valy)
vals(i)    = valj + ima * valy
end do
call elapsed(t2)

t_eval    = (t2-t1)/n
err1      = maxval(abs( vals - vals0) / abs(vals0) )

call write_integer_sep(iw2,nu)
call write_table_next(iw2)
call write_table_double(iw2,err1)
call write_table_next(iw2)
call write_table_double(iw2,t_eval)
call write_table_next(iw2)
call write_table_nextline(iw2)


print "(A45,' =    ',I10.10)","nu",nu
print "(A45,' =    ',I10.10)","Number of points",n
print "(A45,' = ',D25.16)","Average time eval",t_eval
print "(A45,' = ',D25.16)","Max error in J_\dnu + i Y_\dnu",err1

print *,""

end do

deallocate(ts,vals,vals0)


write (iw2,*) "\bottomrule"
write (iw2,*) "\end{tabular}"

close(iw)
close(iw2)

end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Support routines
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine elapsed(t)
implicit double precision (a-h,o-z)
integer*8 i,irate
!
!  Return the time elapsed since some preset moment using the 
!  system_clock routine.
!
call system_clock(i,irate)
dd = i
dd = dd/irate
t  = dd
return
end subroutine


subroutine write_table_double(iw,xx)
implicit double precision (a-h,o-z)

if (xx .eq. 0) then
i1 = 0
i2 = 0
i3 = 0
nexp = 0

write (iw,'(I1,".",I1,I1,"\e{+",I2.2,"}")',advance="no") i1,i2,i3,nexp

return
endif

nexp = 0
x    = xx

if (x .gt. 1) then

do while (x .ge. 10.0d0) 
x    = x / 10.0d0
nexp = nexp + 1
end do


i1 = floor(x)
i2 = floor(x*10-10*i1)
i3 = floor(x*100-100*i1-10*i2)

write (iw,'(I1,".",I1,I1,"\e{+",I2.2,"} ")',advance="no") i1,i2,i3,nexp

else

do while (x .le. 1.0d0) 
x    = x * 10.0d0
nexp = nexp + 1
end do


i1 = floor(x)
i2 = floor(x*10-10*i1)
i3 = floor(x*100-100*i1-10*i2)

write (iw,'(I1,".",I1,I1,"\e{-",I2.2,"}")',advance="no") i1,i2,i3,nexp

endif

end subroutine



subroutine write_table_integer_range(iw,nu1,nu2)
implicit double precision (a-h,o-z)

call write_integer_sep(iw,nu1)
write(iw,"(A)",advance="no") " - "
call write_integer_sep(iw,nu2)

end subroutine



subroutine write_integer_sep(iw,n)
implicit double precision (a-h,o-z)

integer :: idigits(1000)

m      = n
ncount = 0
isep   = 0

if (n .eq. 0) then
write(iw,"(A)",advance="no") "0"
return
endif

do while (m .ge. 1)

if (isep .eq. 3) then
ncount          = ncount+1
idigits(ncount) = -1
isep            = 0
endif

ncount          = ncount+1
idigits(ncount) = mod(m,10)
m               = m /10
isep           = isep+1

end do

do i=1,ncount
id = idigits(ncount-i+1)
if (id .ge. 0) then
write(iw,"(I1)",advance="no") id
else
write(iw,"(A)",advance="no") "\sep,"
endif

end do

end subroutine


subroutine write_table_nextline(iw)
implicit double precision (a-h,o-z)
write (iw,*) " \\"
end subroutine


subroutine write_table_next(iw)
implicit double precision (a-h,o-z)
write (iw,"(' & ')",advance="no") 
end subroutine





program bessel_test
implicit double precision (a-h,o-z)

double precision :: vals(10000)
double complex :: val,val0,ima,val2


call report_on_expansions()
call test_oscillatory("bessel1.txt","table1.tex")
call test_nonoscillatory("bessel2.txt","table2.tex")
call test_nonoscillatory("bessel3.txt","table3.tex")
call test_vals("bessel4.txt","table4.tex")

end program
