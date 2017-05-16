! gaussfilter for time series
subroutine gaussfilter(sigin,sigout,wn,alpha,dt,n)
integer n,n21,i,xi,nptsmax,npts
parameter (nptsmax=131072)
real dt,df,alpha,fac
real freq,fact,filt
real frequp,freqlw,wn
real sigin(nptsmax),sigout(nptsmax)
complex temp(nptsmax)
npts=1
do while(npts.lt.n)
      npts=npts*2
enddo
n21=npts/2+1
df=1.0/npts/dt
temp=cmplx(0.0,0.0)
do i=1,n
      temp(i)=cmplx(sigin(i),0.0)
enddo
call four(temp,npts,-1,dt,df)
fac = sqrt(3.1415927/alpha)
frequp = (1.0+fac)*wn
freqlw = (1.0-fac)*wn
if(freqlw .le. 0.0)freqlw = df
do i=1,n21
       xi = i - 1
       freq = xi *df
       if(freq.ge.freqlw .and. freq.le.frequp)then
               fact = -alpha*((freq-wn)/(wn))**2
               filt = exp(fact)
               temp(i) = filt*temp(i)
       else
               temp(i) = cmplx(0.0,0.0)
       endif
       temp(npts+2-i)=cmplx(0.0,0.0)
enddo
call four(temp,npts,+1,dt,df)
do i=1,n
      sigout(i)=real(temp(i))
enddo
end subroutine
