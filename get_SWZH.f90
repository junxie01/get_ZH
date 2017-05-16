!      measure the amplitude ratio (at each period) between z and r component of Rayleigh wave
!      2016.07.22 junxie
!         2016.07.22 add four alpha algorithms:
!                    1: from chen haopeng
!                    2: cps1
!                    3: cps2
!                    4: ftan
program get_EQ_ZH
use sacio
use getalpha
implicit none
type(sac_head):: sacheadz,sacheadr
integer nptsmax,npmax
integer nper,idva,ip
integer halfpts,npts
integer nb,ne,np,nr,nz,id
integer n1,n2,nn,nno,imark
integer nargs,nsamp,nerr,i,j
character(8) kstnm
character(16) kevnm
character(180) args,sacfile_z,sacfile_r,output
parameter (nptsmax=4000000,npmax=100)
real inpr,inpz               ! instaneous frequency(period)
real snr_r,snr_z
real t1,t2,sd,ccc
real hbr(nptsmax)
real pharr(nptsmax)
real pharz(nptsmax),baz
real pb,pe,pmin,pmax,dp
real amp_z,amp_r,u_r,u_z
real maxnoisez,maxnoiser,snr
real wmax,wmin,ampmax,amplmax
real gfz(nptsmax),gfr(nptsmax)
real a1(nptsmax),a2(nptsmax),zh
real per(npmax),vb,ve,wn(npmax)
real sigz(nptsmax),sigr(nptsmax)
real stla,stlo,evla,evlo,rmax,zmax
real en_gfz(nptsmax),en_gfr(nptsmax)
real begin_tr,begin_tz,dt,width,period,df,dist,evdp
complex hb_r(nptsmax)
complex sig_z(nptsmax),sig_r(nptsmax)
complex sig_zs(nptsmax),sig_rs(nptsmax)
complex sig_zd(nptsmax),sig_rd(nptsmax)

  nno=100                    ! number of trailer noise points
  nargs=iargc()
  if(nargs.ne.6)then
     write(*,*)'Usage: get_EQZH sac_file_z_comp sac_file_r_comp pmin pmax id output'
     stop
  endif
  call getarg(1,sacfile_z)
  call getarg(2,sacfile_r)
  call getarg(3,args)
  read(args,'(bn,f20.0)')pmin
  call getarg(4,args)
  read(args,'(bn,f20.0)')pmax
  call getarg(5,args)
  read(args,'(bn,i20)')id
  call getarg(6,output)


! read z component
  call read_sachead(sacfile_z,sacheadz, nerr)
  call read_sac(sacfile_z,sigz,sacheadz,nerr)
  if(nerr.eq.-1) then
     write(*,*)'STOP: Error in reading file: ',sacfile_z
     stop
  endif
  dist=sacheadz%dist
  kstnm=sacheadz%kstnm
  kevnm=sacheadz%kevnm
  evlo=sacheadz%evlo
  evla=sacheadz%evla
  stla=sacheadz%stla
  stlo=sacheadz%stlo
  evdp=sacheadz%evdp
  baz=sacheadz%baz
  begin_tz=sacheadz%b
! read r component
  call read_sachead(sacfile_r,sacheadr, nerr)
  call read_sac(sacfile_r,sigr,sacheadr,nerr)
  if(nerr.eq.-1) then
     write(*,*)'STOP: Error reading in file: ',sacfile_r
     stop
  endif
  if(sacheadr%delta/=sacheadz%delta)then
     write(*,*)'STOP: The sampling rate is different, please check! '
     stop
  endif
  if(sacheadr%npts/=sacheadz%npts)then
     write(*,*)'STOP: The number of points is different, please check! '
     stop
  endif
  nsamp=sacheadr%npts
  dt=sacheadr%delta
  begin_tr=sacheadr%b
  call hilbert(sigr,nsamp,dt,hbr) ! do hilbert transform to r component

  npts=1
  do while(npts.lt.nsamp)
     npts=npts*2
  enddo
  halfpts=npts/2+1
  df=1.0/(npts*dt)
  sig_z=cmplx(0,0)
  sig_r=cmplx(0,0)
  hb_r=cmplx(0,0)
  sig_z(1:nsamp)=cmplx(sigz(1:nsamp),0)
  sig_r(1:nsamp)=cmplx(sigr(1:nsamp),0)
  call four(sig_z,npts,-1,dt,df)
  call four(sig_r,npts,-1,dt,df)

  idva=0;vb=1.0;ve=5.0
  n1=int((dist/ve-begin_tz)/dt)
  n2=int((dist/vb-begin_tz)/dt)
  if(n1.lt.1) n1=1
  if(n2.gt.nsamp) n2=nsamp
  dp=1.0/dt/2.0;dp=1.0/dp ! the minimum period
! shortest wave length should be one third of distance
  pe=dist/3.0/ve;pb=dt;pb=max(pb,dp)

  call getperr(wn,npmax,pmin,pmax,np,pb,pe)  ! generate period chains
  open(10,file=output)
  write(10,'("# ",1a8,1x,2f10.4,1x,1a16,1x,5f13.4)')adjustl(kstnm),stlo,stla,trim(kevnm),evlo,evla,evdp,dist, baz
  do ip=1,np                                          ! loop over period
       select case(id)
         case (1) ! chen hao peng
         call getalpha_hp( dist,1.0/wn(ip),width, 1)  ! get the alpha parameter
         case (2) ! cps 1
         call getalpha_cps( dist,1.0/wn(ip),width,2)  ! get the alpha parameter
         case (3) ! cps 2
         call getalpha_cps( dist,1.0/wn(ip),width,3)  ! get the alpha parameter
         case (4) ! ftan
         call getalpha_ftan(dist,1.0/wn(ip),width,4)  ! get the alpha parameter
       end select

       call gaussfilter(sigz,gfz,wn(ip),width,dt,nsamp) ! do gaussian filter to z component
       call gaussfilter(hbr, gfr,wn(ip),width,dt,nsamp) ! do gaussian filter to r component
       call envelope(gfz,en_gfz,nsamp,dt)               ! get envelope of z
       call envelope(gfr,en_gfr,nsamp,dt)               ! get envelope of r
       imark=0;zmax=0;rmax=0
       maxnoisez=sqrt(sum(en_gfz(nsamp-nno+1:nsamp)**2)/nno)  !calculate rms of trailer noise with length of nno
       maxnoiser=sqrt(sum(en_gfr(nsamp-nno+1:nsamp)**2)/nno)
       do i=n1,n2                         ! find the max amplitude of the envolope waveform
          if(zmax.lt.en_gfz(i))then       ! zmax is the peak amplitude of the envolope wavform of Z component after gaussfilter
             zmax=en_gfz(i)
             imark=i
          endif
          rmax=amax1(rmax,en_gfr(i))      ! maxmium amplitude of r component at the correponding time
       enddo
       snr_z=zmax/maxnoisez
       snr_r=rmax/maxnoiser               ! calculate the SNR

       nb=imark-int(4.0/wn(ip)/dt)        ! define the surface wave window with half length of 4T (guess it is too long)
       ne=imark+int(4.0/wn(ip)/dt)
       if(nb.lt.1)nb=1
       if(ne.gt.nsamp)ne=nsamp
       nn=ne-nb
       sd=0;a1=0;a2=0;ccc=0
       zh=sum(en_gfz(nb:ne-1)/en_gfr(nb:ne-1))/nn ! calculate average Z/H ratio in a window
       a1(1:nn)=gfz(nb:ne-1)
       a2(1:nn)=gfr(nb:ne-1)
       call correlate(a1,a2,nn,ccc)      ! calculate the cross correlation coefficient
       do i=1,nn
          sd=sd+(zh-en_gfz(i+nb-1)/en_gfr(i+nb-1))**2/nn    
       enddo
       sd=sqrt(sd)                       ! standard error
       call dofilt(sig_z,sig_zs,sig_zd,pharz,dt,df,npts,halfpts,width,idva,ampmax,wn(ip),vb,ve,dist,begin_tz,nsamp) ! do gaussian filter for z and r component in frequency domain
       call dofilt(sig_r,sig_rs,sig_rd,pharr,dt,df,npts,halfpts,width,idva,ampmax,wn(ip),vb,ve,dist,begin_tr,nsamp) ! do gaussian filter for z and r component in frequency domain
       call get_amp_v(sig_zs,sig_zd,pharz,amp_z,u_z,npts,dist,ve,vb,begin_tz,dt,inpz) ! get maximum amplitude
       call get_amp_v(sig_rs,sig_rd,pharr,amp_r,u_r,npts,dist,ve,vb,begin_tr,dt,inpr)
       if(ccc.ne.ccc)ccc=0
       write(10,'(1f7.2,6f10.4,2f7.3,2f9.2, 1f7.3, 1f7.2, 1f10.2, 1f7.2)') 1.0/wn(ip),zmax/rmax,zmax/en_gfr(imark),amp_z/amp_r,&
                                                   snr_z/snr_r,zh,sd,u_z, u_r, snr_z,snr_r,ccc, evdp, dist, baz
enddo
stop
end program
