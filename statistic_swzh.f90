!
!V1. Jun Xie, Jiajun Chong @ WHIGG  2016-07-20
!
!Format of SWZH files for statistical analysis: (1-15)
!1.0/wn(ip),zmax/rmax,zmax/en_gfr(imark),amp_z/amp_r,snr_z/snr_r,zh,sd,u_z, u_r, snr_z,snr_r,ccc, evdp, dist, baz
!
! u_z, u_r: 8,9  (group velocity)
! snr_z,snr_r: 10, 11
! ccc: 12
! evdp, dist, baz: 13, 14, 15
!
program main ! do statistical analysis for Z/H ratio
implicit none
integer nmax,npmax
real smv
parameter(nmax=1000,npmax=100,smv=1.0e-4)
logical ext
real tmpf(15)
real snr_min,cc_min
real pb,pe,pmin,pmax
real dif_phase
real std(npmax),mean_zh(npmax)
real std1(npmax),mean_zh1(npmax)
real period(nmax,npmax),periods(npmax)
real zh(nmax,npmax),wn(npmax),wnin(npmax)
real snrr(nmax,npmax),snrz(nmax,npmax),cc(nmax,npmax)
integer nper,num(npmax),num1(npmax),i,j,izh,ii
integer ip,nzh,flag(nmax,npmax),flag1(nmax,npmax)
character(80) par,tmp
character(80) filelist,output,zh_file(nmax)
if(iargc().ne.6)then
     write(*,*)'Usage:statistic_swzh swzh_filelist pmin pmax snr_min cc_min output_file'
     call exit(-1)
endif
call getarg(1,filelist)
call getarg(2,par)
read(par,'(bn,f20.0)')pmin
call getarg(3,par)
read(par,'(bn,f20.0)')pmax
call getarg(4,par)
read(par,'(bn,f20.0)')snr_min
call getarg(5,par)
read(par,'(bn,f20.0)')cc_min
call getarg(6,output)

inquire(file=filelist,exist=ext)
if(.not.ext)then
     write(*,*)'Please put individual swzh measurement file names into ',trim(filelist)
     call exit(-1)
endif

pe=200;pb=1

call getperr(periods,npmax,pmin,pmax,nper,pb,pe)  ! generate period chains
periods(1:nper)=1.0/periods(1:nper)
nzh=0;std=0;flag=0;flag1=0

open(10,file=filelist)
11 read (10,*,err=12,end=12)tmp ! read filelist
   inquire(file=tmp,exist=ext) 
   if(.not.ext)continue
   nzh=nzh+1                        ! if zh_file exists
   zh_file(nzh)=tmp
   goto 11
12 continue

close(11)

mean_zh=0;mean_zh1=0;zh=0;num=0

do izh=1,nzh         ! loop over z/h ratio file
   !write(*,*)"read file:",trim(zh_file(izh))
   open(20,file=zh_file(izh))        ! read Z/H ratio
   read(20,*)tmp
21 read(20,*,err=22,end=22)(tmpf(i),i=1,15)

   do ii=1,nper                      ! loop over target periods
      if(abs(tmpf(1)-periods(ii)).le.smv)then ! locate the period
         zh(izh,ii)=tmpf(4)
         dif_phase = (tmpf(14)/tmpf(8)-tmpf(14)/tmpf(9))/tmpf(1)*180.0/3.14159265
         dif_phase = abs(dif_phase)
         if(tmpf(12).eq.tmpf(12))then !check CCC
            if(tmpf(10).ge.snr_min.and.tmpf(11).ge.snr_min.and.tmpf(12).ge.cc_min.and.tmpf(4).lt.3.and.&
tmpf(4).gt.0.2 .and. dif_phase .lt. 22.5)then
               num(ii)=num(ii)+1
               mean_zh(ii)=mean_zh(ii)+tmpf(4)
            !write(*,*)tmpf(1),tmpf(4)
               flag(izh,ii)=1
            endif
         endif
         goto 21
      endif
   enddo   ! end loop over period
   goto 21
22 continue
   close(20)
enddo

open(30,file=output)

do ip=1,nper      ! loop over period
   mean_zh(ip)=mean_zh(ip)/num(ip) ! raw mean
   std(ip)=0

   do izh=1,nzh
      if(flag(izh,ip).eq.1)then
         std(ip)=std(ip)+(zh(izh,ip)-mean_zh(ip))**2
      endif
   enddo

   !write(*,*)periods(ip),num(ip)
   std(ip)=sqrt(std(ip)/(num(ip)-1))
   num1(ip)=0
   mean_zh1(ip)=0
   std1(ip)=0

!recount for zh_mean+-std
   do izh=1,nzh
      if(flag(izh,ip).eq.1.and.zh(izh,ip).gt.mean_zh(ip)-std(ip).and.zh(izh,ip).lt.mean_zh(ip)+std(ip))then
         mean_zh1(ip)=mean_zh1(ip)+zh(izh,ip)
         num1(ip)=num1(ip)+1
         flag1(izh,ip)=1
      endif
   enddo
   mean_zh1(ip)=mean_zh1(ip)/num1(ip)
   do izh=1,nzh
      if(flag1(izh,ip).eq.1)then
         std1(ip)=std1(ip)+(zh(izh,ip)-mean_zh1(ip))**2
      endif
   enddo
   std1(ip)=sqrt(std(ip)/(num1(ip)-1))
!     mean_zh: mean of all zh, (std, err_average); mean_zh1: mean of mean_zh+-2*std_zh;
   write(30,'(4f7.3,i5,3f7.3,i5)')periods(ip),mean_zh(ip),std(ip),std(ip)/sqrt(real(num(ip))),num(ip),mean_zh1(ip),std1(ip),&
   std1(ip)/sqrt(real(num1(ip))),num1(ip)
enddo
close(30)
end program
