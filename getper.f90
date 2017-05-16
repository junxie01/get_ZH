subroutine getperr(wn,NX,pmin,pmax,nper,permin,permax) 
!       automatically create periods from the [pmin,pmax] limits
!       pmin    R*4 minimum period
!       pmax    R*4 maximum period
!       nper    I*4 number of periods generated
!       permin  R*4 - minimum period to be used in trace
!       permax  R*4 - maximum period to be used in trace
!----
integer MAXFRQ
parameter (MAXFRQ=100)
real wn(NX)
real pmin, pmax
real permin, permax
integer perrmin,perrmax
integer NX, nper
perrmin=int(max(pmin,permin))
perrmax=int(min(pmax,permax))
nper=1
wn(nper)=1.0/perrmin
p=perrmin
do while(p.lt.perrmax)
   if(p.lt.30)then
      p=p+2
   else
      p=p+4
   endif
   wn(nper+1)=1.0/p
   nper=nper+1
enddo
return
end subroutine
