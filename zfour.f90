        subroutine zfour(zarr,nn,isign,dt,df) 
!-----
!     THE input is a complex array
!     which has numbers stored in memory as
!     R1, I1, R2, I2, ..., Rnn, Inn
!     where nn must be a power of 2 R and I are the real and imaginary
!     parts of the complex number
!
!     For isign -1 this is a complex time series
!     For isign +1 this is a complex frequency series with
!        index 1 (in fortran corresponding to f=0
!              2                              f=df
!            nn/2 + 1                         f = 1/2dt = Nyquist
!            nn - 1                           f = -2df
!             nn                              f = -df

!-----
!     the cooley-tookey fast fourier transform in usasi basic fortran
!     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
!     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
!     dimensional complex array (i.e., the real and imaginary parts of
!     datc are located immediately adjacent in storage, such as fortran
!     places them) whose length nn is a power of two.  isign
!     is +1 or -1, giving the sign of the transform.  transform values
!     are returned in array datc, replacing the input datc.  the time is
!     proportional to n*log2(n), rather than the usual n**2
!     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
!     b is the number of bits in the floating point fraction.
!
!     the program computes df from dt, dt from df and checks to see
!     if they are consistent. In addition, the transforms are multiplied
!     by dt or df to make the results dimensionally correct
!
!     This is a slightly modified version of the original Brenner routine
!     The changes were to add physical dimensions to the transform
!     and to make it all complex
!-----
        complex zarr(*) 
        integer nn, isign
        real dt, df

        complex ztemp
!-----
!       ensure that the dt and df are defined and
!       consistent
!-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
!-----
!       now begin the transform
!-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
!-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
!-----
!       use trig relations to compute the next sin/cos
!       without actually calling sin() or cos()
!-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
!-----
!       transform is done
!-----
   10   continue 
!-----
!     give the arrays the proper physical dimensions
!-----
        if(isign.lt.0)then
!-----
!             time to frequency domain
!-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
!-----
!             frequency to time domain
!-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end

