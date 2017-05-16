        subroutine four(data,nn,isign,dt,df) 
        integer NPTSER
        parameter (NPTSER=132000)
        dimension data(NPTSER*2) 
        n = 2 * nn 
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
        j = 1 
        do 5 i=1,n,2 
!       if(i-j)1,2,2 
        if(i .lt. j) then
            go to 1
        else
            go to 2
        endif
    1 tempr = data(j) 
        tempi = data(j+1) 
        data(j) = data(i) 
        data(j+1)=data(i+1) 
        data(i) = tempr 
        data(i+1) = tempi 
    2 m = n/2 
!    3 if(j-m) 5,5,4 
    3 continue
        if(j.le.m) then
            go to 5
        else 
            go to 4
        endif
    4 j = j-m 
        m = m/2 
!       if(m-2)5,3,3 
        if(m.lt.2)then
            go to 5
        else
            go to 3
        endif
    5 j=j+m 
        mmax = 2 
!    6 if(mmax-n) 7,10,10 
    6 continue
        if(mmax .lt. n)then
            go to 7
        else if(mmax .ge. n)then
            go to 10
        endif
    7 istep= 2 *mmax 
        theta = 6.283185307/float(isign*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 m=1,mmax,2 
        do 8 i=m,n,istep 
        j=i+mmax 
        tempr=wr*data(j)-wi*data(j+1) 
        tempi=wr*data(j+1)+wi*data(j) 
        data(j)=data(i)-tempr 
        data(j+1)=data(i+1)-tempi 
        data(i)=data(i)+tempr 
    8 data(i+1) = data(i+1)+tempi 
        tempr = wr 
        wr = wr*wstpr-wi*wstpi + wr 
    9 wi = wi*wstpr+tempr*wstpi + wi 
        mmax = istep 
        go to 6 
   10 continue 
        if(isign.lt.0) go to 1002 
!     frequency to time domain 
        do 1001 iiii = 1,n 
 1001 data(iiii) = data(iiii) * df 
        return 
 1002 continue 
!     time to frequency domain 
        do 1003 iiii = 1,n 
 1003 data(iiii) = data(iiii) * dt 
        return 
        end 

