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
!  The code in this file tests compares the performance of the bessel_eval routine
!  with that the code of Amos described in 
!
!    Amos, D. E. Algorithm 644: a portable package for Bessel functions of a 
!     complex argument and nonnegative order. ACM Transactions on Mathematica 
!     Software 3 (1986), 265–273.
!
!  Amos' code is not included in this package; the user must obtain it on his
!  or her own in order to run these experiments.
!
!  In addition to outputing results to the console, this program generates
!  LaTeX source code for the fourth table appearing in the numerical
!  experiments section of the preprint
!
!     James Bremer, "An algorithm for the rapid numerical evaluation of Bessel 
!       functions of real orders and arguments."  arXiv: 
!
!  as well as a Mathematica script for generating the graphs appearing in
!  the only figure appearing that section.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine test_vals(filename1,filename2)
use besseleval
implicit double precision (a-h,o-z)
character(len=*) :: filename1,filename2
double precision, allocatable :: ts(:)

! double precision, allocatable :: errs(:) , errs_amos(:)
! double precision, allocatable :: times(:), times_amos(:)

double complex, allocatable   ::  vals(:),vals0(:),vals2(:)
double complex :: ima
ima = (0.0d0,1.0d0)

print *,"-------------------------------------------------------------------------------------"
print *,"Testing the values of J_nu and Y_nu ( ",filename1," )"
print *,"-------------------------------------------------------------------------------------"
print *,""




iw2  = 21
open(iw2, FILE=filename2)


write(iw2,*) "\begin{tabular}{l@{\hspace{2em}}cc@{\hspace{2em}}cc}"
write(iw2,*) "\toprule"


write(iw2,"(A)",advance="no") "& \multicolumn{2}{c}{\large \tt bessel\_eval} "
write(iw2,"(A)")              "& \multicolumn{2}{c}{\large \tt Amos' code} \\"

write(iw2,"(A)") "\addlinespace[.25em]"

write(iw2,"(A)",advance="no") "$n$ & Maximum relative & Average evaluation  "
write(iw2,"(A)")              " & Maximum relative & Average evaluation \\ "

write(iw2,"(A)",advance="no") "      & error in $H_n$ & time (in seconds) "
write(iw2,"(A)")              " & error in $H_n$ & time (in seconds) \\"

write(iw2,*) "\midrule"



iw  = 20
open(iw,  FILE=filename1)
read(iw,*) nchunks
read(iw,*) n

allocate(ts(n),vals0(n),vals(n),vals2(n))

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


!
!  Try to avoid cache effects by clearing val
!

do i=1,n
vals(i) = 0.0d0
end do



!
!  Run the bessel_eval code
!

call elapsed(t1)
do i=1,n
t   = ts(i)
call bessel_eval(dnu,t,alpha,alphader,beta1,beta2,valj,valy)
vals(i)    = valj + ima * valy
end do
call elapsed(t2)



t_eval1   = (t2-t1)/n
err1      = maxval(abs( vals - vals0) / abs(vals0) )

!
!  Try to avoid cache effects by clearing val2
!

do i=1,n
vals2(i) = 0.0d0
end do


!
!  Now Amos' code
!

call elapsed(t1)
do i=1,n
t   = ts(i)
! call zbesj(t, 0.0d0,dnu,1,1,valj, CYI, NZ, IERR)
! call zbesy(t,0.0d0, dnu,1, 1, valy, dwork1, NZ, dwork2,dwork3, IERR)
call zbesh(t,0.0d0,dnu,1,1,1,valj,valy,nz,ierr)

vals2(i)    = valj + ima * valy
end do
call elapsed(t2)

t_eval2   = (t2-t1)/n
err2      = maxval(abs( vals2 - vals0) / abs(vals0) )



call write_integer_sep(iw2,nu)
call write_table_next(iw2)
call write_table_double(iw2,err1)
call write_table_next(iw2)
call write_table_double(iw2,t_eval1)
call write_table_next(iw2)

if (abs(err2) .gt. 1.0d-1) then
call write_table_string(iw2,"-")
call write_table_next(iw2)
call write_table_string(iw2,"-")
else
call write_table_double(iw2,err2)
call write_table_next(iw2)
call write_table_double(iw2,t_eval2)
endif

call write_table_nextline(iw2)


print "(A45,' =    ',I10.10)","nu",nu
print "(A45,' =    ',I10.10)","Number of points",n
print "(A45,' = ',D25.16)","Average time (bessel_eval) ",t_eval1
print "(A45,' = ',D25.16)","Max error in H_nu (bessel_eval) ",err1

if (abs(err2) .gt. 1.0d-1) then
print "(A45,' = ',A)","Average time (Amos) ","N/A"
print "(A45,' = ',A)","Max error in H_nu (Amos) ","N/A"
else
print "(A45,' = ',D25.16)","Average time (Amos) ",t_eval2
print "(A45,' = ',D25.16)","Max error in H_nu (Amos) ",err2
endif

print *,""


end do

deallocate(ts,vals,vals0,vals2)


write (iw2,*) "\bottomrule"
write (iw2,*) "\end{tabular}"

close(iw)
close(iw2)

end subroutine



subroutine test_figures(filename1,filename2)
use besseleval
implicit double precision (a-h,o-z)
character(len=*) :: filename1,filename2
double precision, allocatable :: ts(:)

double complex, allocatable   ::  vals(:),vals0(:),vals2(:)
double complex :: ima
ima = (0.0d0,1.0d0)

pi   = acos(-1.0d0)
eps0 = epsilon(0.0d0)

print *,"-------------------------------------------------------------------------------------"
print *,"Testing the values of J_nu and Y_nu ( ",filename1," )"
print *,"-------------------------------------------------------------------------------------"
print *,""

iw   = 20
open(iw, FILE=filename1)

iw2  = 21
open(iw2, FILE=filename2)

if (eps0 .lt. 1.0d-17) then
write(iw2,*) "eps0       = 2^(-112);"
else
write(iw2,*) "eps0       = 2^(-52);"
endif

write(iw2,*) "f[t_]      = HankelH1[dnu,t];";
write(iw2,*) "g[dnu_,t_] = Log10[Abs[t * f'[t]/f[t]]*eps0];"
write(iw2,*) ""

read(iw,*) M
read(iw,*) n

allocate(ts(n),vals0(n),vals(n),vals2(n))

!
!  Read in the data
!

do j=1,M

!
!  Read in the data stored in the text file
!

read(iw,*) dnu
read(iw,*) a
read(iw,*) b

do i=1,n
read (iw,*) t
read (iw,*) valr
read (iw,*) vali
vals0(i) = valr+ima*vali
ts(i)    = t
end do

!
!  Use bessel_eval
!

write (iw2,*) "MM = 100;"
write (iw2,*) "nu = ",dnu," ;"
write (iw2,*) "a = ",a," ;"
write (iw2,*) "b = ",b," ;"

! write (iw2,*) "a  = nu/2;"
! write (iw2,*) "b  = 1000*nu;"
write(iw2,*) "X1 = Table[{Exp[Log[a] + (Log[b]-Log[a])*(i-1)/(MM-1)],"
write(iw2,*) "  N[g[nu, Exp[Log[a] + (Log[b] - Log[a])*(i - 1)/(MM - 1)]]]},{i,1,MM}];"
write (iw2,*) "Y1 = ListLinePlot[X1, PlotStyle -> {Black, Thickness[.005]}, "
write (iw2,*) "GridLines -> Automatic, Frame -> True, BaseStyle -> {FontWeight -> Bold, FontSize -> 10},"
if (eps0 .lt. 1.0d-17) then
write (iw2,*) "PlotRange->{{a,b},{-36,-10}}];"
else
write (iw2,*) "PlotRange->{{a,b},{-16,-6}}];"
endif



write(iw2,*) "X2 = {" 

do i=1,n

t = ts(i)
call bessel_eval(dnu,t,alpha,alphader,beta1,beta2,valj,valy)
vals(i) = valj+ima*valy

err = log10(abs(vals(i)-vals0(i)) / abs(vals0(i)))

!print *,dnu,t,err

if (i .lt. n) then
write (iw2,"('{',F45.15,' , ',F45.15,' },')") t,err
else
write (iw2,"('{',F45.15,' , ',F45.15,' }};')") t,err
endif

end do




!
!  Amos
!

write(iw2,*) "X3 = {" 

do i=1,n
t = ts(i)
call zbesh(t,0.0d0,dnu,1,1,1,valj,valy,nz,ierr)
! call zbesj(t, 0.0d0,dnu,1,1,valj, CYI, NZ, IERR)
! call zbesy(t,0.0d0, dnu,1, 1, valy, dwork1, NZ, dwork2,dwork3, IERR)
vals2(i)    = valj + ima * valy

err = log10(abs(vals2(i)-vals0(i)) / abs(vals0(i)))

if (i .lt. n) then
write (iw2,"('{',F45.15,' , ',F45.15,' },')") t,err
else
write (iw2,"('{',F45.15,' , ',F45.15,' }};')") t,err
endif

end do



write (iw2,*) "Y2 = ListPlot[X2,PlotStyle->{PointSize[.005]},GridLines->Automatic,Frame->True,"
write (iw2,*) "  BaseStyle->{FontWeight->Bold,FontSize->10},"
if (eps0 .lt. 1.0d-17) then
write (iw2,*) "PlotRange->{{a,b},{-36,-10}}];"
else
write (iw2,*) "PlotRange->{{a,b},{-16,-6}}];"
endif
write (iw2,*) "Z = Show[Y1,Y2];"

if (eps0 .lt. 1.0d-17) then
write (iw2,"(A,I1.1,A)") 'Export["BesselPlot',j,'_16.pdf",Z];'
else
write (iw2,"(A,I1.1,A)") 'Export["BesselPlot',j,'.pdf",Z];'
endif


write (iw2,*) "Y3 = ListPlot[X3,PlotStyle->{PointSize[.005]},GridLines->Automatic,Frame->True,"
write (iw2,*) "  BaseStyle->{FontWeight->Bold,FontSize->10},"
if (eps0 .lt. 1.0d-17) then
write (iw2,*) "PlotRange->{{a,b},{-36,-10}}];"
else
write (iw2,*) "PlotRange->{{a,b},{-16,-6}}];"
endif

write (iw2,*) "Z = Show[Y1,Y3];"

if (eps0 .lt. 1.0d-17) then
write (iw2,"(A,I1.1,A)") 'Export["AmosPlot',j,'_16.pdf",Z];'
else
write (iw2,"(A,I1.1,A)") 'Export["AmosPlot',j,'.pdf",Z];'
endif

end do


close(iw2)
close(iw)

! write(iw2,*) "\begin{tabular}{l@{\hspace{2em}}cc@{\hspace{2em}}cc}"
! write(iw2,*) "\toprule"


! write(iw2,"(A)",advance="no") "& \multicolumn{2}{c}{\large \tt bessel\_eval} "
! write(iw2,"(A)")              "& \multicolumn{2}{c}{\large \tt Amos' code} \\"

! write(iw2,"(A)") "\addlinespace[.25em]"

! write(iw2,"(A)",advance="no") "$\nu$ & Maximum relative & Average evaluation  "
! write(iw2,"(A)")              " & Maximum relative & Average evaluation \\ "

! write(iw2,"(A)",advance="no") "      & error in $H_n$ & time (in seconds) "
! write(iw2,"(A)")              " & error in $H_n$ & time (in seconds) \\"

! write(iw2,*) "\midrule"



! iw  = 20
! open(iw,  FILE=filename1)
! read(iw,*) nchunks
! read(iw,*) n

! allocate(ts(n),vals0(n),vals(n),vals2(n))

! do ichunk=1,nchunks

! read(iw,*) nu
! read(iw,*) dtime
! dnu = nu


! do i=1,n
! read (iw,*) t
! read (iw,*) valj
! read (iw,*) valy

! ts(i)    = t
! vals0(i) = valj + ima * valy
! end do


! !
! !  Run the bessel_eval code
! !

! call elapsed(t1)
! do i=1,n
! t   = ts(i)
! call bessel_eval(dnu,t,alpha,alphader,beta1,beta2,valj,valy)
! vals(i)    = valj + ima * valy
! end do
! call elapsed(t2)

! t_eval1   = (t2-t1)/n
! err1      = maxval(abs( vals - vals0) / abs(vals0) )

! !
! !  Now Amos' code
! !

! call elapsed(t1)
! do i=1,n
! t   = ts(i)
! call zbesh(t,0.0d0,dnu,1,1,1,valj,valy,nz,ierr)
! vals2(i)    = valj + ima * valy
! end do
! call elapsed(t2)

! t_eval2   = (t2-t1)/n
! err2      = maxval(abs( vals2 - vals0) / abs(vals0) )


! call write_integer_sep(iw2,nu)
! call write_table_next(iw2)
! call write_table_double(iw2,err1)
! call write_table_next(iw2)
! call write_table_double(iw2,t_eval1)
! call write_table_next(iw2)

! if (abs(err2) .gt. 1.0d-1) then
! call write_table_string(iw2,"-")
! call write_table_next(iw2)
! call write_table_string(iw2,"-")
! else
! call write_table_double(iw2,err2)
! call write_table_next(iw2)
! call write_table_double(iw2,t_eval2)
! endif

! call write_table_nextline(iw2)


! print "(A45,' =    ',I10.10)","nu",nu
! print "(A45,' =    ',I10.10)","Number of points",n
! print "(A45,' = ',D25.16)","Average time (bessel_eval) ",t_eval1
! print "(A45,' = ',D25.16)","Max error in H_nu (bessel_eval) ",err1

! if (abs(err2) .gt. 1.0d-1) then
! print "(A45,' = ',A)","Average time (Amos) ","N/A"
! print "(A45,' = ',A)","Max error in H_nu (Amos) ","N/A"
! else
! print "(A45,' = ',D25.16)","Average time (Amos) ",t_eval2
! print "(A45,' = ',D25.16)","Max error in H_nu (Amos) ",err2
! endif

! print *,""

! end do

! deallocate(ts,vals,vals0,vals2)


! write (iw2,*) "\bottomrule"
! write (iw2,*) "\end{tabular}"

! close(iw)
! close(iw2)

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



subroutine write_table_string(iw,str)
implicit double precision (a-h,o-z)
character (len = *) :: str
write (iw,'(A)',advance="no") str
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
use besseleval
implicit double precision (a-h,o-z)

eps0 = epsilon(0.0d0)

call elapsed(t1)
call bessel_eval_init(dsize)
call elapsed(t2)

write (*,"(A,F4.2,A)") "time required to read bessel data = ",t2-t1," seconds"
write (*,"(A,F4.2,A)") "memory occupied by bessel data    = ",dsize," MB"
print *,""

!
!  conduct the fourth experiment
!


!
!  conduct the fifth experiment
! 

if (eps0 .gt. 1.0d-17) then
call test_figures("bessel4.txt","generate_figs.m")
endif

if (eps0 .lt. 1.0d-17) then
call test_vals("bessel5.txt","table5_16.tex")
else
call test_vals("bessel5.txt","table5.tex")
endif

end program
