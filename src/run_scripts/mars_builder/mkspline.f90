subroutine mkspline( x, nx, fspl, ibcxmin, bcxmin, ibcxmax,                   &
                     bcxmax, ilinx, iflag )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Make a 2-coefficient 1d spline.
!
! *** Details:
!
!     Only 2 coefficients, the data and its 2nd derivative, are needed to
!     fully specify a spline.  See e.g. Numerical Recipies in Fortran-77
!     (2nd edition) chapter 3, section on cubic splines.
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,                intent(in)    :: nx     ! no. of data points
      real( gpu ),            intent(in)    :: x(nx)  ! x axis data (ascending)

!     f(1,*): data in; f(2,*): coeffs out :

      real( gpu ),            intent(inout) :: fspl( 2, nx )                   
!
!     f(1,j) = f(x(j))  on input (unchanged on output)
!     f(2,j) = f''(x(j)) (of interpolating spline) (on output).
!
!     Boundary conditions:

!     b.c. flag @ x=xmin=x(1)  :

      integer,                intent(in)    :: ibcxmin 
      real( gpu ),            intent(in)    :: bcxmin ! b.c. data @xmin

!     b.c. flag @ x=xmax=x(nx) :

      integer,                intent(in)    :: ibcxmax                   
      real( gpu ),            intent(in)    :: bcxmax ! b.c. data @xmax

!     ibcxmin=-1 -- periodic boundary condition
!                   (bcxmin,ibcxmax,bcxmax are ignored)
!
!                   the output spline s satisfies
!                   s'(x(1))=s'(x(nx)) ..and.. s''(x(1))=s''(x(nx))
!
!     If non-periodic boundary conditions are used, then the xmin and xmax
!     boundary conditions can be specified independently:
!
!     ibcxmin (ibcxmax) = 0 -- this specifies a "not a knot" boundary
!                   condition, see "cubsplb.for".  This is a common way
!                   for inferring a "good" spline boundary condition
!                   automatically from data in the vicinity of the
!                   boundary.  ... bcxmin (bcxmax) are ignored.
!
!     ibcxmin (ibcxmax) = 1 -- boundary condition is to have s'(x(1))
!                   ( s'(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!     ibcxmin (ibcxmax) = 2 -- boundary condition is to have s''(x(1))
!                   ( s''(x(nx)) ) match the passed value bcxmin (bcxmax).
!
!     ibcxmin (ibcxmax) = 3 -- boundary condition is to have s'(x(1))=0.0
!                   ( s'(x(nx))=0.0 )
!
!     ibcxmin (ibcxmax) = 4 -- boundary condition is to have s''(x(1))=0.0
!                   ( s''(x(nx))=0.0 )
!
!     ibcxmin (ibcxmax) = 5 -- boundary condition is to have s'(x(1))
!                   ( s'(x(nx)) ) match the 1st divided difference
!                   e.g. at x(1):  d(1)/h(1), where
!                                  d(j)=f(1,j+1)-f(1,j)
!                                  h(j)=x(j+1)-x(j)
!
!     ibcxmin (ibcxmax) = 6 -- BC is to have s''(x(1)) ( s''(x(nx)) )
!                   match the 2nd divided difference
!                   e.g. at x(1):
!                   e(1) = [d(2)/h(2) - d(1)/h(1)]/(0.5*(h(1)+h(2)))
!
!     ibcxmin (ibcxmax) = 7 -- BC is to have s'''(x(1)) ( s'''(x(nx)) )
!                   match the 3rd divided difference
!                   e.g. at x(1): [e(2)-e(1)]/(0.33333*(h(1)+h(2)+h(3)))
!

!     =1: hint, x axis is ~evenly spaced:

      integer,                intent(out)   :: ilinx
                     
!     let dx[avg] = (x(nx)-x(1))/(nx-1)
!     let dx[j] = x(j+1)-x(j), for all j satisfying 1.le.j.lt.nx
!
!     if for all such j, abs(dx[j]-dx[avg]).le.(1.0e-3*dx[avg]) then
!     ilinx=1 is returned, indicating the data is (at least nearly)
!     evenly spaced.  Even spacing is useful, for speed of zone lookup,
!     when evaluating a spline.
!
!     if the even spacing condition is not satisfied, ilinx=2 is returned.

      integer,                intent(out)   :: iflag ! exit code, 0=OK

!     An error code is returned if the x axis is not strict ascending,
!     or if nx.lt.4, or if an invalid boundary condition specification was
!     input.
!
!
!     N.B. this routine calls traditional 4 coefficient spline software, and
!     translates the result to 2 coefficient form.
!
!     This could be done more efficiently but we decided out of conservatism
!     to use the traditional software.

!     Local variables:

      real( gpu ),            allocatable   :: fspl4(:,:) ! 4-spline
      real( gpu ),            allocatable   :: wk   (:)   ! cspline workspace
      integer                               :: i, inwk

      allocate( fspl4( 4, nx ), wk(nx) )

!     Make the traditional call:

      do i = 1, nx
         fspl4( 1, i ) = fspl( 1, i )
         fspl ( 2, i ) = 0.0_gpu ! for now
      enddo

      inwk = nx

!     Boundary conditions imposed by cspline...

      call cspline( x, nx, fspl4, ibcxmin, bcxmin, ibcxmax, bcxmax,           &
                    wk, inwk, ilinx, iflag )

      if( iflag == 0 )then

!         Copy the output -- careful of end point.

          do i = 1, nx-1
             fspl( 2, i ) = 2.0_gpu*fspl4( 3, i )
          enddo

          fspl( 2, nx ) = 2.0_gpu*fspl4( 3, nx-1 ) +                          &
              ( x(nx)-x(nx-1) )*6.0_gpu*fspl4( 4, nx-1 )
      endif

      deallocate( fspl4, wk )

end subroutine mkspline

