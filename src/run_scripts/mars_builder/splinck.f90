subroutine splinck( x, nx, ilinx, tol, iflag )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Check if a grid is strictly ascending and if it is evenly spaced
!     to w/in tol.
!
! *** Calling variables:
!
!     x		GPU             i/p     : x-axis                            [-]
!     nx        I               i/p     : x-axis dimension                  [-]
!     ilinx     I               o/p     : x-axis equal spacing flag         [-]
!     tol       GPU             i/p     : Spacing tolerance                 [-]
!     iflag     I               o/p     : Error flag                        [-]
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,                intent(in)    :: nx
      real( gpu ),            intent(in)    :: x(nx) ! grid to check
      integer,                intent(out)   :: ilinx  ! =1 if evenly spaced
                                                      ! =2 O.W.
      real( gpu ),            intent(in)    :: tol    ! spacing check tolerance
      integer,                intent(out)   :: iflag  ! =0 if OK

      real( gpu )                           :: dxavg, zeps, zdiffx, zdiff
      integer                               :: ix

!     iflag = 1 is returned if x(1...nx) is NOT STRICTLY ASCENDING.

!     Initialize error flag:

      iflag = 0
      
      ilinx = 1
      
      if( nx <= 1 )return

      dxavg = ( x(nx) - x(1) )/real( nx-1, gpu )
      zeps  = abs( tol*dxavg )

      do ix = 2, nx
      
         zdiffx = ( x(ix) - x( ix-1 ) )
         
         if( zdiffx <= 0.0_gpu )iflag = 2
         
         zdiff  = zdiffx-dxavg
         
         if( abs(zdiff) > zeps )then
             ilinx = 2
         endif
         
      enddo

end subroutine splinck

