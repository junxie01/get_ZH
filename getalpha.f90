!***************************************************************
!
!    variable Gaussian Width
!
!***************************************************************

module getalpha
contains

subroutine getalpha_hp(distance,period,alpha,id) ! chenhaopeng
real distance,period,alpha
integer :: id
   alpha = distance/8000.0*((6.8146+51.2927/period)**2)

return
end subroutine

subroutine getalpha_cps(distance,period,alpha,id)
real distance,period,alpha
integer i,id
parameter(npts = 6)
real x(npts), a(npts), width(npts),dist(npts)
data x /25., 35., 45., 55., 60., 300./
data a /0.025,0.011,0.009,0.008,0.007125,0.003/
data width /3, 6.5, 12.5, 25, 50, 75./
data dist /100, 250, 500, 750, 1000, 3000./
            
if (id.eq.2) then     
      if(period .le. x(1)) then
           alpha = 0.5 * distance * a(1)
           return
      else if(period .gt. x(npts))then
           alpha = 0.5 * distance * a(npts)
           return
      end if
!      
      i = 1
10    i = i + 1
      if(period .gt. x(i)) go to 10
      
      dx = period - x(i-1)
      dadx = (a(i) - a(i-1)) / (x(i) - x(i-1))
      alpha = a(i-1) + dx * dadx
      alpha = alpha * 0.5 * distance
else if(id.eq.3)then
      i=1
11    i=i+1
      if (distance.gt.dist(i).and.i.lt.6)goto 11
      alpha=width(i)
      if(distance.gt.3000)alpha=100
endif
return
end subroutine

subroutine getalpha_ftan(distance,period,alpha,id)
real distance,alpha
integer id
alpha=2*20.0*sqrt(distance/1000.)
end subroutine getalpha_ftan

end module getalpha
