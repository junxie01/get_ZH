      function correlate(str1,str2,n,cc)
      dimension str1(*),str2(*)
      integer n,i
      real ua/0/,ub/0/,uab/0/,ava/0/,avb/0/
      real cc
      do i=1,n
         ava=ava+str1(i)/n
         avb=avb+str2(i)/n
      enddo
      do i=1,n
         ua=ua+(str1(i)-ava)**2
         ub=ub+(str2(i)-avb)**2
         uab=uab+(str1(i)-ava)*(str2(i)-avb)
      enddo
      cc=uab/sqrt(ua*ub) 
      return
      end function
