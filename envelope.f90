subroutine envelope(sigin,sigout,n,dt)
integer n,nmax,npts,i,halfpts
parameter(nmax=2000000)
real sigin(n),sigout(n)
real dt,df,w
real twopi
complex temp(nmax)
twopi = acos(-1.0)
npts=1
do while (npts.lt.n)
       npts=npts*2
enddo
halfpts=npts/2
df=1.0/npts/dt
temp=cmplx(0.0,0.0)
do i=1,n
      temp(i)=cmplx(sigin(i),0.0)
enddo
call four(temp,npts,-1,dt,df) ! do fft
do i=1,halfpts
       temp(i)=2.0*temp(i)*dt
       temp(i+halfpts)=cmplx(0.0,0.0)
enddo
call four(temp,npts,1,dt,df)
do i=1,n
       sigout(i)=cabs(temp(i))/dt
enddo
end subroutine 
