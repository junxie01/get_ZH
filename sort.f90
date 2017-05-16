       subroutine sort(x,key,n)
!-----
!     Starting with x(1) ,,, x(n)
!     return   the xarray sorted in increasing order
!     also return the pointers key to the initial array. 
!     For example given x = [ 3, 1, 2 ]
!     the returned values are
!                       x = [ 1, 2, 3 ]        
!                     key = [ 2, 3, 1 ]
!-----
!      Reference: http://en.wikipedia.org/wiki/Bubble_sort
!-----
       integer n
       real x(n)
       integer key(n)

       do i=1,n
           key(i) = i
       enddo
       do i = n, 1, -1
           do j = 1 , i -1
               if(x(j) .gt. x(j+1))then
                   tmp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = tmp
                   ktmp = key(j)
                   key(j) = key(j+1)
                   key(j+1) = ktmp
                endif
           enddo
       enddo
       return
       end
