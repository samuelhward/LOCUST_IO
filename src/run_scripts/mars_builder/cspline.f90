subroutine cspline( x, nx, fspl, ibcxmin, bcxmin, ibcxmax, bcxmax,            &
                    wk, iwk, ilinx, iflag )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     cspline -- dmc 15 Feb 1999
!
!     A standard interface to the 1d spline setup routine
!     modified dmc 3 Mar 2000 -- to use Wayne Houlberg's v_spline code.
!     New BC options added.
!
! *** Details:

!     N.B. wk(...) array is not used unless ibcxmin=-1 (periodic spline
!     evaluation)
!
!     This routine computes spline coefficients for a 1d spline --
!     evaluation of the spline can be done by cspeval.for subroutines
!     or directly by inline code.
!
!     The input x axis x(1...nx) must be strictly ascending, i.e.
!     x(i+1).gt.x(i) is required for i=1 to nx-1.  This is checked and
!     iflag=1 is set and the routine exits if the test is not satisfied.
!
!     On output, ilinx=1 is set if, to a reasonably close tolerance,
!     all grid spacings x(i+1)-x(i) are equal.  This allows a speedier
!     grid lookup algorithm on evaluation of the spline.  If on output
!     ilinx=2, this means the spline x axis is not evenly spaced.
!
!     The input data for the spline are given in f[j] = fspl(1,j).  The
!     output data are the spline coefficients fspl(2,j),fspl(3,j), and
!     fspl(4,j), j=1 to nx.  The result is a spline s(x) satisfying the
!     boundary conditions and with the properties
!
!        s(x(j)) = fspl(1,j)
!        s'(x) is continuous even at the grid points x(j)
!        s''(x) is continuous even at the grid points x(j)
!
!     The formula for evaluation of s(x) is:
!
!        let dx = x-x(i), where x(i).le.x.le.x(i+1).  Then,
!        s(x)=fspl(1,i) + dx*(fspl(2,i) +dx*(fspl(3,i) + dx*fspl(4,i)))
!
!     ==>boundary conditions.  Complete specification of a 1d spline
!     requires specification of boundary conditions at x(1) and x(nx).
!
!     This routine provides 4 options:
!
!    -1 ***** PERIODIC BC
!       ibcxmin=-1  --  periodic boundary condition.  This means the
!       boundary conditions s'(x(1))=s'(x(nx)) and s''(x(1))=s''(x(nx))
!       are imposed.  Note that s(x(1))=s(x(nx)) (i.e. fspl(1,1)=fspl(1,nx))
!       is not required -- that is determined by the fspl array input data.
!       The periodic boundary condition is to be preferred for periodic
!       data.  When splining periodic data f(x) with period P, the relation
!       x(nx)=x(1)+n*P, n = the number of periods (usually 1), should hold.
!       (ibcxmax, bcxmin, bcxmax are ignored).
!
!     If a periodic boundary condition is set, this covers both boundaries.
!     for the other types of boundary conditions, the type of condition
!     chosen for the x(1) boundary need not be the same as the type chosen
!     for the x(nx) boundary.
!
!     0 ***** NOT A KNOT BC
!       ibcxmin=0 | ibcxmax=0 -- this specifies a "not a knot" boundary
!       condition -- see cubsplb.for.  This is a common way for inferring
!       a "good" spline boundary condition automatically from data in the
!       vicinity of the boundary.  (bcxmin | bcxmax are ignored).
!
!     1 ***** BC:  SPECIFIED SLOPE
!       ibcxmin=1 | ibcxmax=1 -- boundary condition is to have s'(x(1)) |
!       s'(x(nx)) match the passed value (bcxmin | bcxmax).
!
!     2 ***** BC:  SPECIFIED 2nd DERIVATIVE
!       ibcxmin=2 | ibcxmax=2 -- boundary condition is to have s''(x(1)) |
!       s''(x(nx)) match the passed value (bcxmin | bcxmax).
!
!     3 ***** BC:  SPECIFIED SLOPE = 0.0
!       ibcxmin=3 | ibcxmax=3 -- boundary condition is to have s'(x(1)) |
!       s'(x(nx)) equal to ZERO.
!
!     4 ***** BC:  SPECIFIED 2nd DERIVATIVE = 0.0
!       ibcxmin=4 | ibcxmax=4 -- boundary condition is to have s''(x(1)) |
!       s''(x(nx)) equal to ZERO.
!
!     5 ***** BC:  1st DIVIDED DIFFERENCE
!       ibcxmin=5 | ibcxmax=5 -- boundary condition is to have s'(x(1)) |
!       s'(x(nx)) equal to the slope from the 1st|last 2 points
!
!     6 ***** BC:  2nd DIVIDED DIFFERENCE
!       ibcxmin=6 | ibcxmax=6 -- boundary condition is to have s''(x(1)) |
!       s''(x(nx)) equal to the 2nd derivative from the 1st|last 3 points
!
!     7 ***** BC:  3rd DIVIDED DIFFERENCE
!       ibcxmin=7 | ibcxmax=7 -- boundary condition is to have s'''(x(1)) |
!       s'''(x(nx)) equal to the 3rd derivative from the 1st|last 4 points
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,                intent(in)    :: nx
      real( gpu ),            intent(in)    :: x(nx)      ! x axis
      real( gpu ),            intent(inout) :: fspl(4,nx) ! Spline data
      integer,                intent(in)    :: ibcxmin    ! x(1)  BC flag 
      real( gpu ),            intent(in)    :: bcxmin     ! x(1)  BC data 
      integer,                intent(in)    :: ibcxmax    ! x(nx) BC flag 
      real( gpu ),            intent(in)    :: bcxmax     ! x(nx) BC data 
      integer,                intent(in)    :: iwk
      real( gpu ),            intent(in)    :: wk(iwk)    ! workspace size>=nx
      integer,                intent(out)   :: ilinx      ! even spacing flag
      integer,                intent(out)   :: iflag      ! =0 means OK

!     Local variables:

      integer                               :: inum, i

!     Initialize error flag:

      iflag = 0

!     Error checks:

      if( nx < 2 )then
          write(io(1),'('':cspline     : ERROR : at least 2 x points&
                         & required.'')')
          iflag = 1
      endif
      
      call ibc_ck( ibcxmin, 'cspline     : ERROR :', 'xmin', -1, 7, iflag )
      if( ibcxmin >= 0 )call ibc_ck( ibcxmax, 'cspline     : ERROR :',        &
                                     'xmax', 0, 7, iflag )

!     x axis check:

      call splinck( x, nx, ilinx, 1.0e-03_gpu, iflag )

      if( iflag /= 0 )then
      
          iflag = 2
      
          write(io(1),'('':cspline     : ERROR : x axis not strict&
                         & ascending'')')
      endif

      if( ibcxmin == -1 )then
          inum = nx
          if( iwk < inum )then
              write(io(1),100) inum, iwk, nx
              iflag = 3
          endif
      endif

      if( iflag /= 0 )return

!     OK -- evaluate spline:

      if(     ibcxmin == 1 )then
          fspl(2,1) = bcxmin
      elseif( ibcxmin == 2 )then
          fspl(3,1) = bcxmin
      endif

      if(     ibcxmax == 1 )then
          fspl(2,nx) = bcxmax
      elseif( ibcxmax == 2 )then
          fspl(3,nx) = bcxmax
      endif

      call v_spline( ibcxmin, ibcxmax, nx, x, fspl, wk )

      do i = 1, nx
         fspl(3,i) = 0.5_gpu*fspl(3,i)
         fspl(4,i) =         fspl(4,i)/6.0_gpu
      enddo

      return

!     Format statements:

100   format(' :cspline:     : ERROR : workspace too small.  need:  ',i6,     &
             ' got:  ',i6/'  (need = nx, nx=',i6)

end subroutine cspline

