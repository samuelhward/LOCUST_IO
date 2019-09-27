subroutine ibc_ck( ibc, slbl, xlbl, imin, imax, iflag )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Check that spline routine ibc flag is in range.
!
! *** Calling variables:
!
!     ibc       I               i/p     : Flag value                        [-]
!     slbl      char            i/p     : Subroutine name                   [-]
!     xlbl      char            i/p     : Axis label                        [-]
!     imin      I               i/p     : Min allowed value                 [-]
!     imax      I               i/p     : Max      -"-                      [-]
!     iflag     I               o/p     : Error flag                        [-]
!
! *** Error flags:
!
!     iflag = 0         : Return status OK
!             1         : ERROR : Overwritten to 1 if error is detected
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,                intent(in)    :: ibc   ! flag value
      character*(*),          intent(in)    :: slbl  ! subroutine name
      character*(*),          intent(in)    :: xlbl  ! axis label
      integer,                intent(in)    :: imin  ! min allowed value
      integer,                intent(in)    :: imax  ! max allowed value
      integer,                intent(inout) :: iflag ! set =1 if error detected

      if( ( ibc < imin ) .or. ( ibc > imax ) ) then
         iflag = 1
         write(io(1),100) slbl, xlbl, ibc, imin, imax
      endif

      return
      
!     Format statements:
      
100   format( ':', A, ' -- ibc', A, ' = ', I9,                                 &
              ' out of range ', I2, ' to ', I2 ) 
      
end subroutine ibc_ck

