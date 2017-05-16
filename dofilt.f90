subroutine dofilt(datas,data,datad,pharr,dt,df,n,n21,alpha,idva,ampmx,wn,pumin,pumax,dist,t0,nsamp)
!-----
!       perform narrow band pass filter
!
!       datas   C   Complex array of saved spectrum
!       data    C   Complex array of filtered signal
!       datad   C   Complex array of derivative of filtered signal
!       pharr   R
!       dt  R   sampling interval, seconds
!       df  R   frequency sampling interval
!       n   I   number of points in time series (power of 2)
!       n21 I   n/2 + 1
!       alpha   R   Gaussian filter parameter
!       idva    I   0 displacment time series
!               1 velocity time series
!               2 acceleration time series
!               we want to have displacement spectrum
!       ampmx   R   Maximum amplitude for this filter frequency
!       wn  R   Filter frequency
!       pumin  R*4      minimum group velocity for plots
!       pumax  R*4      maximum group velocity for plots
!       dist    R   epicentral distance in km
!       t0  R   travel time of first sample
!       nsamp   I number of samples in original time series
!-----
        integer NPTSER, NPTFRQ
!        parameter (NPTSER=132000, NPTFRQ=65000)
        complex datas(n),data(n),datad(n)
        real pharr(n)
        complex s, s2

        integer nsamp
        real pumin, pumax
!-----
!       define cutoff for filter
!-----
        fac = sqrt(3.1415927/alpha)
        frequp = (1.0+fac)*wn
        freqlw = (1.0-fac)*wn
        if(freqlw .le. 0.0)freqlw = df
        do 1002 i=1,n21
                xi = i - 1
                freq = xi *df
                if(freq.ge.freqlw .and. freq.le.frequp)then
                    fact = -alpha*((freq-wn)/(wn))**2
                    filt = exp(fact)
                    s = cmplx(0.0, 6.2831853*freq)
                    s2 = s * s
!-----
!       convert from velocity to displacement, or from
!       acceleration to displacement here
!       since we are bandpass filtgering we never divide by zero??
!-----
                    data(i) = filt*datas(i)
                    if(idva.eq.1)then
                        data(i) = data(i)/s
                    else if(idva.eq.2)then
                        data(i) = data(i)/s2
                    endif
!------
!       set up array to compute time derivative which is
!       necessary for computing the instantaneous period/frequency
!-----
                    datad(i) = data(i) * s
                else
                    data(i) = cmplx(0.0,0.0)
                    datad(i) = cmplx(0.0,0.0)
                endif
!-----
!           zero out negative frequencies
!-----
                if(i.gt.1)then
                    data(n+2-i)=cmplx(0.0,0.0)
                    datad(n+2-i)=cmplx(0.0,0.0)
                endif
 1002   continue
        call four(data,n,+1,dt,df)
        call four(datad,n,+1,dt,df)
        
!-----
!           save envelope for contour plot
!           wn is the filter frequency
!-----
        
        afac = sqrt(alpha/3.1415927)/wn
        amx = 0.0
!-----
!       compute the spectral amplitude and instantaneous
!       frequency and store into first half of the
!       complex data and datad arrays
!       Note we want the maximum amplitude within a given group velocity window
!-----
        tuend = dist/pumin - t0
        iuend = tuend/dt + 1
        if(iuend.gt.nsamp)iuend = nsamp
        if(iuend.lt.1)iuend = 1
        tustrt = dist/pumax - t0
        iustrt = tustrt/dt + 1
        if(iustrt.gt.nsamp)iustrt = nsamp
        if(iustrt.lt.1)iuend = 1
        do 5006 KK = 1,n
                xamp = cabs(data(KK))*afac
                if(xamp.gt.0)then
                    xfrq = aimag(datad(KK)/data(KK))
                    xpha = atan2(aimag(data(KK)),real(data(KK)))
                else
                    xfrq = 1.0e+5
                    xpha = 0.0
                endif
                if(KK.ge.iustrt .and. KK.le. iuend)then
                    if(xamp.gt.amx)amx=xamp
                endif
                data(KK)=cmplx(xamp,0.0)
                datad(KK)=cmplx(xfrq,0.0)
                pharr(KK) = xpha
 5006   continue
!-----
!       Now systematically pack the complex array
!       in a manner that is the same as the equivalence
!       Basically
!       R 0 R I R I R I
!       where we read in R 0 or R I and then want to
!       store as
!       R R R R
!       For each  complex number that we readm, we will with
!       two successive real values
!-----
        KJ = 1
        do 5007 KK=1,n
                xfrq = real(datad(KK))
                xamp = real(data (KK))
                if(mod(KK,2).eq.1)then
                    data (KJ) = cmplx(xamp,0.0)
                    datad(KJ) = cmplx(xfrq,0.0)
                else
                    data (KJ) = cmplx(real(data (KJ)),xamp)
                    datad(KJ) = cmplx(real(datad(KJ)),xfrq)
                    KJ = KJ + 1
                endif
                
 5007   continue
        !call wsac0('mft_z.sac',dum,real(data),nerr)
        if(amx.gt.ampmx)ampmx = amx
        return 
        end
