      subroutine zonfind(x,nx,zxget,i)
!
      use prec_mod
      
      implicit none
      
      integer, intent(in) :: nx
      real( gpu ) :: x(nx),zxget
      real( gpu ) :: dx
      integer i, nxm, i1, i2, ij, ii
!
!  find index i such that x(i).le.zxget.le.x(i+1)
!
!  x(1...nx) is strict increasing and x(1).le.zxget.le.x(nx)
!  (this is assumed to already have been checked -- no check here!)
!
      nxm=nx-1
      if((i.lt.1).or.(i.gt.nxm)) then
         i1=1
         i2=nx-1
         go to 10
      endif
!
      if(x(i).gt.zxget) then
!  look down
         dx=x(i+1)-x(i)
         if((x(i)-zxget).gt.4*dx) then
            i1=1
            i2=i-1
            go to 10
         else
            i2=i-1
            do ij=i2,1,-1
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            i=1
            return
         endif
      else if(x(i+1).lt.zxget) then
!  look up
         dx=x(i+1)-x(i)
         if((zxget-x(i+1)).gt.4*dx) then
            i1=i+1
            i2=nxm
            go to 10
         else
            i2=i+1
            do ij=i2,nxm
               if((x(ij).le.zxget).and.(zxget.le.x(ij+1))) then
                  i=ij
                  return
               endif
            enddo
            ij=nxm
            return
         endif
      else
!  already there...
         return
      endif
!
!---------------------------
!  binary search
!
 10   continue
!
      if(i1.eq.i2) then
! found by proc. of elimination
         i=i1
         return
      endif
!
      ii=(i1+i2)/2
!
      if(zxget.lt.x(ii)) then
         i2=ii-1
      else if(zxget.gt.x(ii+1)) then
         i1=ii+1
      else
! found
         i=ii
         return
      endif
!
      go to 10
!
      end
 
