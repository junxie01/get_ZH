!subroutine get_amp_v(x,xfrq,pharr,a,u,nsamp,dist,ve,vb,t0,dt,per)
subroutine get_amp_v(amp,phas,pharr,a,u,nsamp,dist,ve,vb,t0,dt,per)
parameter(MAX4=10)
integer nsamp,i,kk,j
complex amp(nsamp),phas(nsamp)
real dxdt(nsamp),per
real x(nsamp),xfrq(nsamp),pharr(nsamp)
real aa(MAX4),uu(MAX4)
real ufr(MAX4),p(MAX4)
real time,phase,t0
real ampl,vel,freq
integer iustrt,iuend
kk=0
do i=1,nsamp/2
       j=(i-1)*2
       x(j+1)=real(amp(i))
       x(j+2)=aimag(amp(i))
       xfrq(j+1)=real(phas(i))
       xfrq(j+2)=aimag(phas(i))
enddo
do i=1,nsamp-1
       dxdt(i)=x(i+1)-x(i)
enddo
!write(*,*)'dist=',dist
iustrt=max(int((dist/ve-t0)/dt)+1,1)
iuend =min(int((dist/vb-t0)/dt)+1,nsamp)
!write(*,*)iustrt,iuend
do i=iustrt+1,iuend-1
      if(sign(1.0,dxdt(i-1)).ne.sign(1.0,dxdt(i)))then
            pv = (dxdt(i-1) - 0.0)/(dxdt(i-1) - dxdt(i))
            time = t0 + (i-2 +pv) * dt
            if(time.gt.0.0)then
                   vel = dist/time
                   freq = pv*xfrq(i-1) +(1-pv)*xfrq(i)
                   ampl = pv*x(i-1) + (1-pv)*x(i)
                   phase = pv*pharr(i-1) +(1-pv)*pharr(i)
                   kk = kk + 1
                   call srt(aa,uu,ufr,p,vel,ampl,freq,phase,kk,MAX4)
            endif
      endif
enddo
if(kk.gt.MAX4)kk=MAX4
!do i=1,kk
!       write(*,*)per,uu(i),aa(i),freq
!enddo
!write(*,*)per,uu(1),aa(1)
a=aa(1)
u=uu(1)
per=ufr(1)
return
end subroutine

subroutine srt(a,u,ufr,p,vel,ampl,frq2,phase,k,MAX4)

!-----
!       subroutine to maintain a list of the MAX4
!       largest values of ampl in a(i)
!       also saving the corresponding vel in v(i)
!
!       This is being done to avoid have two very large
!       arrays in overhead and because we do not need
!       to do a full sort
!-----
parameter (MAXVEL=10)
real*4 a(MAXVEL), u(MAXVEL), ufr(MAXVEL), vel, ampl
real*4 p(MAXVEL)
integer*4 MAX4, k
integer*4 MAX41
integer*4 key(11)
real*4 tmp(11)
real*4 tamp(11), tvel(11), tfrq(11), tp(11)
MAX41 = MAX4 + 1
if(k .eq. 1)then
      do i=1,MAX4
            a(i) = 0.0
            u(i) = 0.0
            ufr(i) = 0.0
            p(i) = 0.0
      enddo
endif
if(k.lt.MAX4)then
      kup = k
else
      kup = MAX4
endif
!-----
!       assume amplitudes arranged in decreasing order
!-----
!-----
!       we now know that the value will replace one of the amplitudes,
!       at least the lowest
!-----
do i=1,MAX4
      tmp(i) = a(i)
      tamp(i) = a(i)
      tvel(i) = u(i)
      tfrq(i) = ufr(i)
      tp(i) = p(i)
enddo
tmp(MAX41) = ampl
tamp(MAX41) = ampl
tvel(MAX41) = vel
tfrq(MAX41) = frq2
tp(MAX41) = phase
call sort(tmp,key,MAX41)
do i= 1,MAX4
      kk = key(MAX41 +1 - i)
      a(i) = tamp(kk)
      u(i) = tvel(kk)
      ufr(i) = tfrq(kk)
      p(i) = tp(kk)
enddo
return
end
