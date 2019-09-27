subroutine evspline_(xget,x,nx,ilinx,f,ict,fval,ier)
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Wrapper to PSPLINES routine evspline allowing for "extrapolation" out
!     of range.
!
! *** Calling variables:
!
!     xget      GPU             i/p     : Target of the interpolation       [-]
!     x         GPU             i/p     : Ordered x grid                    [-]
!     nx        I               i/p     : Grid size                         [-]
!     ilinx     GPU             i/p     : ilinx=1 => assume x evenly spaced [-]
!     f         GPU             i/p     : Function data                     [-]
!     ict       I               i/p     : Code specifying output desired    [-]
!     fval      GPU             o/p     : Output data                       [-]
!     ier       I               o/p     : Error code                        [-]
!
!
! *** Author:
!
!     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
!
! *** Error flags:               
!
!     iflag = 0         : Return status OK.
!             1         : ERROR : 
!
! *** Creation Date:
!
!     18/05/2016
!
!------------------------------------------------------------------------------
      use prec_mod

      implicit none

      integer,     intent(in)  :: nx
      real( gpu ), intent(in)  :: xget
      real( gpu ), intent(in)  :: x(nx)
      integer,     intent(in)  :: ilinx
      real( gpu ), intent(in)  :: f(2,nx)
      integer,     intent(in)  :: ict(3)
      real( gpu ), intent(out) :: fval(*)
      integer,     intent(out) :: ier

!     Local variables:

      real( gpu ) :: a_, b_, c_, d_, dx

      ier = 0

      if( xget >= x(1) .and. xget <= x(nx) )then

          call evspline(xget,x,nx,ilinx,f,ict,fval,ier)
          return

      elseif( xget < x(1) )then      

          dx =  x(2)-x(1)
          a_ = (f(2,2)-f(2,1))/(6.0_gpu*dx)
          b_ = (f(2,1)-6.0_gpu*a_*x(1))/2.0_gpu
          c_ =((f(1,2)-a_*x(2)**3-b_*x(2)**2) - (f(1,1)-a_*x(1)**3-b_*x(1)**2))/dx
          d_ = (f(1,1)-a_*x(1)**3-b_*x(1)**2) -c_*x(1)

      else

          dx =  x(nx)-x(nx-1)
          a_ = (f(2,nx  )-f(2,nx-1))/(6.0_gpu*dx)
          b_ = (f(2,nx-1)-6.0_gpu*a_*x(nx-1))/2.0_gpu
          c_ =((f(1,nx  )-a_*x(nx)**3-b_*x(nx)**2) - (f(1,nx-1)-a_*x(nx-1)**3-b_*x(nx-1)**2))/dx
          d_ = (f(1,nx-1)-a_*x(nx-1)**3-b_*x(nx-1)**2) -c_*x(nx-1)

      endif
      
      if( ict(1)==1 )then
          fval(1) = a_*xget**3 + b_*xget**2 + c_*xget + d_
      endif
      if( ict(2)==1 )then
          fval(2) = 3.0_gpu*a_*xget**2 + 2.0_gpu*b_*xget + c_
      endif
      if( ict(3)==1 )then
          fval(3) = 6.0_gpu*a_*xget + 2.0_gpu*b_
      endif
      
end subroutine evspline_

