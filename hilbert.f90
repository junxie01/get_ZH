      subroutine hilbert(data,n,dt,sig1) 
      integer limit
      parameter(limit=131104)
      dimension data(*)
      complex sig(limit),sgn
      real sig1(n),df
      integer n,n12,i,nn
      real t
      sig=cmplx(0.0,0.0)
      do i=1,n
         sig(i)=cmplx(data(i),0.0)
      enddo
      nn=1
      do while(nn<n) 
               nn=nn*2
      enddo
      df=1.0/(nn*dt)
!      call four1(sig,nn,-1)
      call zfour(sig,nn,-1,dt,df)
!      sig=sig*dt
      n12=nn/2 
!      sig(1)    =cmplx(0.0,0.0)
!      sig(n12+1)=cmplx(0.0,0.0)
!      sig(n)=cmplx(0.0,0.0)
      do i=1,nn
!             if (i==1)sgn=cmplx(0.0,0.0)
              if (i==1)then
!                 sig(i)=sig(i)/2.0
                   sgn=cmplx(1/2.0,0)
              else if (i==2)then
                   sgn=cmplx(0,0)
!             if (i==n12+1)sgn=cmplx(1,0.0)
              else if (i==n12+1)then
                   sig(i)=cmplx(real(sig(i)),0.0)
                   sgn=cmplx(1,0)
              else if (i==n12+2)then
                   sgn=cmplx(0,0)
              else if (i>2.and.i<n12+1)then
                   sgn=cmplx(0,-2)
              else 
                   sgn=cmplx(0,0)
              endif
              sig(i)=sig(i)*sgn
      enddo
!      call four1(sig,nn,1)
      call zfour(sig,nn,1,dt,df)
      !sig1=real(sig)
      do i=1,n
             !sig1(i)=real(sig(i))*df
             sig1(i)=real(sig(i))
      enddo
      end subroutine
