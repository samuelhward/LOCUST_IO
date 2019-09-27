subroutine v_spline( k_bc1, k_bcn, n, x, f, wk )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Evaluate the coefficients for a 1d cubic interpolating spline
!
! *** Details:
!
!     References:
!                 Forsythe, Malcolm, Moler, Computer Methods for Mathematical
!                 Computations, Prentice-Hall, 1977, p.76.
!                 Engeln-Muellges, Uhlig, Numerical Algorithms with Fortran,
!                 Springer, 1996, p.251.
!
! *** Calling variables:
!
!     Input:
!
!     k_bc1 : option for BC at x1 = x(1)
!
!           =-1 periodic, ignore k_bcn
!           = 0 not-a-knot
!           = 1 s'(x1) = input value of f(2,1)
!           = 2 s''(x1) = input value of f(3,1)
!           = 3 s'(x1) = 0.0
!           = 4 s''(x1) = 0.0
!           = 5 match first derivative to first 2 points
!           = 6 match second derivative to first 3 points
!           = 7 match third derivative to first 4 points
!           = else use not-a-knot
!
!     k_bcn : option for boundary condition at xn = x(n)
!
!           = 0 not-a-knot
!           = 1 s'(xn) = input value of f(2,n)
!           = 2 s''(xn) = input value of f(3,n)
!           = 3 s'(xn) = 0.0
!           = 4 s''(xn) = 0.0
!           = 5 match first derivative to last 2 points
!           = 6 match second derivative to lasst 3 points
!           = 7 match third derivative to last 4 points
!           = else use knot-a-knot
!
!     n     : number of data points or knots - ( n >= 2 )
!     x(n)  : abscissas of the knots in strictly increasing order
!     f(1,i): ordinates of the knots
!     f(2,1): input value of s'(x1) for k_bc1=1
!     f(2,n): input value of s'(xn) for k_bcn=1
!     f(3,1): input value of s''(x1) for k_bc1=2
!     f(3,n): input value of s''(xn) for k_bcn=2
!     wk(n) : scratch work area for periodic BC
!
!     Output:
!     f(2,i): s'(x(i))
!     f(3,i): s''(x(i))
!     f(4,i): s'''(x(i))
!
! *** Comments:
!
!     s(x)=f(1,i)+f(2,i)*(x-x(i))+f(3,i)*(x-x(i))**2/2!
!          +f(4,i)*(x-x(i))**3/3! for x(i).le.x.le.x(i+1)
!
!     W_SPLINE can be used to evaluate the spline and its derivatives
!     The cubic spline is twice differentiable (C2)
!
! *** Modifications -- dmc 18 Feb 2010:
!
!     Deal with n.lt.4 -- the general tridiagonal spline method 
!     does not have the right formulation for n==3 "not a knot" or periodic
!     boundary conditions, nor for n==2 with any boundary conditions.
!
!     Apply boundary conditions even for n==2, when the "spline" is really
!     just a single cubic polynomial.
!     In this case, several boundary condition (BC) options are mapped to 
!     BC option 5, 1st divided difference.  If 5 is used on both sides
!     of an n==2 "spline" you get a linear piece which is what the old
!     code always gave, regardless of BC option settings.  The following
!     BC controls are mapped to 5:  
!        periodic (-1)
!        not a knot (0) (for n=2 no grid point exists for knot location).
!        option (5) is preserved
!        options 6 and 7 -- mapped to (5); higher divided differences 
!        need n>2; in fact 7 needs n>3; for n=3 option 6 is substituted.
!
!     The boundary condition screening is done at the start of the code;
!     passed controls k_bc1 and k_bcn are mapped to i_bc1 and i_bcn.
!
!     ALSO: for n=3, "not a knot" from both left and right needs special
!           interpretation, since the 2 boundary conditions overlap.  The
!           chosen interpretation is a parabolic fit to the 3 given data 
!           points. And so f''' = 0 and f'' = constant.  If "not a knot" is
!           used on one side only, the solution is a single cubic piece and
!           special code is also needed.
!     ALSO: for n=3, "periodic" boundary condition needs special code; this
!           is added.
!
!  *** Bug Fixes -- dmc 24 Feb 2004:
!
!      (a) fixed logic for not-a-knot:
!            !    Set f(3,1) for not-a-knot
!                      if(k_bc1.le.0.or.k_bc1.gt.7)then ...
!          instead of
!            !    Set f(3,1) for not-a-knot
!                      if(k_bc1.le.0.or.k_bc1.gt.5)then ...
!          and similarly for logic after cmt
!            !    Set f(3,n) for not-a-knot
!          as required since k_bc*=6 and k_bc*=7 are NOT not-a-knot BCs.
!
!      (b) the BCs to fix 2nd derivative at end points did not work if that
!          2nd derivative were non-zero.  The reason is that in those cases
!          the off-diagonal matrix elements nearest the corners are not
!          symmetric; i.e. elem(1,2).ne.elem(2,1) and 
!          elem(n-1,n).ne.elem(n,n-1) where I use "elem" to refer to
!          the tridiagonal matrix elements.  The correct values for the
!          elements is:   elem(1,2)=0, elem(2,1)=x(2)-x(1)
!                         elem(n,n-1)=0, elem(n-1,n)=x(n)-x(n-1)
!          the old code in effect had these as all zeroes.  Since this
!          meant the wrong set of linear equations was solved, the
!          resulting spline had a discontinuity in its 1st derivative
!          at x(2) and x(n-1).  Fixed by introducing elem21 and elemnn1
!          to represent the non-symmetric lower-diagonal values.  Since
!          elem21 & elemnn1 are both on the lower diagonals, logic to 
!          use them occurs in the non-periodic forward elimination loop
!          only.  DMC 24 Feb 2004.
!
! *** Author:
!
!     W.A.Houlberg, D.McCune
!
! *** Creation date:
!
!     March 2000
!
!------------------------------------------------------------------------------

      use prec_mod

      implicit none

      integer,                intent(in)    :: k_bc1
      integer,                intent(in)    :: k_bcn
      integer,                intent(in)    :: n

      real( gpu ),            intent(in)    :: x (  *)
      real( gpu ),            intent(inout) :: wk(  *)
      real( gpu ),            intent(inout) :: f (4,*)

!     Local variables

      integer                               :: i, ib, imax, imin
      real( gpu )                           :: d1, dn, b1, bn, q, t, hn
      real( gpu )                           :: elem21, elemnn1
      integer                               :: i_bc1 ! screened BC cntrls
      integer                               :: i_bcn ! screened BC cntrls

      integer                               :: iord1, iord2 
      real( gpu )                           :: h, f0, fh 
      real( gpu )                           :: h1,h2, h3, dels
      real( gpu )                           :: f1,f2, f3, aa, bb 
      integer                               :: i3knots, i3perio

!     Screen the BC options (DMC Feb. 2010...)

      i_bc1 = k_bc1
      i_bcn = k_bcn

!     Outside [-1:7] -> not-a-knot:

      if( ( i_bc1 < -1 ) .or. ( i_bc1 > 7 ) )i_bc1 = 0 

!     Outside [ 0:7] -> not-a-knot:

      if( ( i_bcn <  0 ) .or. ( i_bcn > 7 ) )i_bcn = 0   

      if( i_bc1 == -1 )i_bcn = -1 ! periodic BC

      i3knots = 0
      i3perio = 0
      
      if( n == 3 )then
          i_bc1 = min( 6, i_bc1 )
          i_bcn = min( 6, i_bcn )
          if( i_bc1 ==  0 )i3knots = i3knots + 1
          if( i_bcn ==  0 )i3knots = i3knots + 1
          if( i_bc1 == -1 )i3perio = 1
      endif

      if( n == 2 )then
          if( i_bc1 == -1 )then
              i_bc1 = 5
              i_bcn = 5
          endif

          if( ( i_bc1 == 0 ) .or. ( i_bc1 > 5 ) )i_bc1 = 5
          if( ( i_bcn == 0 ) .or. ( i_bcn > 5 ) )i_bcn = 5

          if( ( i_bc1 == 1 ) .or. ( i_bc1 == 3 ) .or. ( i_bc1 == 5 ) )then
                iord1 =  1  ! 1st derivative match on LHS
          else
                iord1 =  2  ! 2nd derivative match on LHS
          endif

          if( ( i_bcn == 1 ) .or. ( i_bcn == 3 ) .or. ( i_bcn == 5 ) )then
                iord2 =  1  ! 1st derivative match on RHS
          else
                iord2 =  2  ! 2nd derivative match on RHS
          endif
      endif

!     Set default range:

      imin = 1
      imax = n

!     Set first and second BC values:

      d1   = 0.0_gpu
      b1   = 0.0_gpu
      dn   = 0.0_gpu
      bn   = 0.0_gpu
      
      if(     i_bc1 == 1 )then
        d1 = f(2,1)
      elseif( i_bc1 == 2 )then
        b1 = f(3,1)
      elseif( i_bc1 == 5 )then
        d1 =           ( f(1,2)-f(1,1) )/( x(2)-x(1) )
      elseif( i_bc1 == 6 )then
        b1 = 2.0_gpu*( ( f(1,3)-f(1,2) )/( x(3)-x(2) )                        &
                      -( f(1,2)-f(1,1) )/( x(2)-x(1) ) )/( x(3)-x(1) )
      endif
      
      if(     i_bcn == 1 )then
        dn = f(2,n)
      elseif( i_bcn == 2 )then
        bn = f(3,n)
      elseif( i_bcn == 5 )then
        dn = ( f(1,n)-f(1,n-1) )/( x(n)-x(n-1) )
      elseif( i_bcn == 6 )then
        bn = 2.0_gpu*( ( f(1,n  )-f(1,n-1) )/( x(n  )-x(n-1) )                &
                      -( f(1,n-1)-f(1,n-2) )/( x(n-1)-x(n-2) ) )/             &
                                             ( x(n  )-x(n-2) )
      endif
      
!     Clear f(2:4,n):

      f(2,n) = 0.0_gpu
      f(3,n) = 0.0_gpu
      f(4,n) = 0.0_gpu
      
      if( n == 2 )then

         if( ( i_bc1 == 5 ) .and. ( i_bcn == 5 ) )then
         
!           Coefficients for n=2 (this was the original code):

            f(2,1) = ( f(1,2)-f(1,1) )/( x(2)-x(1) )
            f(3,1) = 0.0_gpu
            f(4,1) = 0.0_gpu
            f(2,2) = f(2,1)
            f(3,2) = 0.0_gpu
            f(4,2) = 0.0_gpu
           
         elseif( ( iord1 == 1 ) .and. ( iord2 == 1 ) )then
        
!           LHS: match d1 for 1st deriv; RHS: match dn for 1st deriv.

            f(2,1) = d1
            f(2,2) = dn
            h      = ( x(2)-x(1) )
            f0     = f(1,1)
            fh     = f(1,2)

            ! setting xx = x-x(1),
            ! f = c1*xx**3 + c2*xx**2 + d1*xx + f0
            !   -->  c1*h**3   + c2*h**2 = fh - f0 - d1*h
            !   and  3*c1*h**2 + 2*c2*h  = dn - d1
            ! f' = 3*c1*xx*2 + 2*c2*xx + d1
            ! f'' = 6*c1*xx + 2*c2
            ! f''' = 6*c1

            ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2:

!           2*c2:

            f(3,1) = ( 3.0_gpu*(fh-f0)/(h*h) - (2.0_gpu*d1 + dn)/h)*2.0_gpu 

!           6*c1:

            f(4,1) = (-2.0_gpu*(fh-f0)/(h*h*h) + (d1 + dn)/(h*h))*6.0_gpu  

            f(4,2) = f(4,1)
            f(3,2) = f(4,1)*h + f(3,1)

         elseif( ( iord1 == 1 ) .and. ( iord2 == 2 ) )then
        
!           LHS: match d1 for 1st deriv; RHS: match bn for 2nd deriv.

            f(2,1) = d1
            f(3,2) = bn
            h      = ( x(2)-x(1) )
            f0     = f(1,1)
            fh     = f(1,2)

            ! setting xx = x-x(1),
            ! f = c1*xx**3 + c2*xx**2 + d1*xx + f0
            !   -->  c1*h**3   + c2*h**2 = fh - f0 - d1*h
            !   and  6*c1*h    + 2*c2    = bn
            ! f' = 3*c1*xx*2 + 2*c2*xx + d1
            ! f'' = 6*c1*xx + 2*c2
            ! f''' = 6*c1

            ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

!           2*c2:

            f(3,1) = (-bn/4.0_gpu + 3.0_gpu*(fh-f0)/                          &
                     (2.0_gpu*h*h) - 3.0_gpu*d1/(2.0_gpu*h))*2.0_gpu       

!           6*c1:
 
            f(4,1) = (bn/(4.0_gpu*h) - (fh-f0)/(2.0_gpu*h*h*h) +              &
                      d1/(2.0_gpu*h*h))*6.0_gpu

            f(4,2) = f(4,1)
            f(2,2) = f(4,1)*h*h/2 + f(3,1)*h + d1

         elseif( ( iord1 == 2 ) .and. ( iord2 == 1 ) )then

            ! LHS: match b1 for 2nd deriv; RHS: match dn for 1st deriv.
           
            f(3,1) = b1
            f(2,2) = dn
            h      = ( x(2)-x(1) )
            f0     = f(1,1)
            fh     = f(1,2)

            ! setting xx = x-x(1), 
            ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
            !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
            !   and  3*c1*h**2 + c3   = dn - b1*h
            ! f' = 3*c1*xx*2 + b1*xx + c3
            ! f'' = 6*c1*xx + b1
            ! f''' = 6*c1

            ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

!           c3:

            f(2,1) = 3.0_gpu*(fh-f0)/(2.0_gpu*h) - b1*h/4.0_gpu - dn/2.0_gpu

!           6*c1:

            f(4,1) = (dn/(2.0_gpu*h*h) - (fh-f0)/(2.0_gpu*h*h*h) -            &
                      b1/(4.0_gpu*h))*6.0_gpu

            f(4,2) = f(4,1)
            f(3,2) = f(4,1)*h + f(3,1)

         elseif( ( iord1 == 2 ) .and. ( iord2 == 2 ) )then
           
!           LHS: match b1 for 2nd deriv; RHS: match bn for 2nd deriv.

            f(3,1) = b1
            f(3,2) = bn
            h      = ( x(2)-x(1) )
            f0     = f(1,1)
            fh     = f(1,2)

            ! setting xx = x-x(1), 
            ! f = c1*xx**3 + (b1/2)*xx**2 + c3*xx + f0
            !   -->  c1*h**3   + c3*h = fh - f0 - b1*h**2/2
            !   and  6*c1*h           = bn - b1
            ! f' = 3*c1*xx*2 + b1*xx + c3
            ! f'' = 6*c1*xx + b1
            ! f''' = 6*c1

            ! solve 2x2 system for c1 -> f(4,1)/6 and c2 -> f(3,1)/2

!           c3:

            f(2,1) = (fh-f0)/h -b1*h/3.0_gpu -bn*h/6.0_gpu

!           6*c1:

            f(4,1) = (bn-b1)/h

            f(4,2) = f(4,1)
            f(2,2) = f(4,1)*h*h/2 + b1*h + f(2,1)

         endif

      elseif( i3perio == 1 )then
      
!        Special case - nx=3 periodic spline:

         h1     = x(2)-x(1)
         h2     = x(3)-x(2)
         h      = h1 + h2

         dels   = ( f(1,3)-f(1,2) )/h2 - ( f(1,2)-f(1,1) )/h1

         f(2,1) = ( f(1,2)-f(1,1) )/h1 + ( h1*dels )/h
         f(3,1) =  -6.0_gpu*dels/h
         f(4,1) =  12.0_gpu*dels/(h1*h)

         f(2,2) = ( f(1,3)-f(1,2) )/h2 - ( h2*dels )/h
         f(3,2) =   6.0_gpu*dels/h
         f(4,2) = -12.0_gpu*dels/(h2*h)

         f(2,3) = f(2,1)
         f(3,3) = f(3,1)
         f(4,3) = f(4,1)

      elseif( i3knots == 2 )then
      
!        Special case - nx=3, not-a-knot on both sides:

         h1     = x(2)-x(1)
         h2     = x(3)-x(2)
         h      = h1 + h2

!        This is just a quadratic fit through 3 pts:

         f1     = f(1,1)-f(1,2)
         f2     = f(1,3)-f(1,2)

!        Quadratic around origin at (x(2),f(1,2))
!                  aa*h1**2 - bb*h1 = f1
!                  aa*h2**2 + bb*h2 = f2

         aa     = ( f2*h1 + f1*h2 )/( h1*h2*h )
         bb     = ( f2*h1*h1 - f1*h2*h2 )/( h1*h2*h )

         f(4,1:3) = 0.0_gpu    ! f''' = 0 (quadratic polynomial)
         f(3,1:3) = 2.0_gpu*aa ! f''  = const

         f(2,  1) = bb - 2.0_gpu*aa*h1
         f(2,  2) = bb
         f(2,  3) = bb + 2.0_gpu*aa*h2

      elseif( i3knots == 1 )then

!        Special cases - nx=3, not-a-knot on single side:

         if( ( i_bc1 == 1 ) .or. ( i_bc1 == 3 ) .or. ( i_bc1 == 5 ) )then

!           f' LHS condition; not-a-knot RHS
!           d1 = f' LHS BC

            h2    = x(2)-x(1)
            h3    = x(3)-x(1)

            f2    = f(1,2)-f(1,1)
            f3    = f(1,3)-f(1,1)

!           Find cubic aa*xx**3 + bb*xx**2 + d1*xx
!           satisfying aa*h2**3 + bb*h2**2 + d1*h2 = f2
!           and aa*h3**3 + bb*h3**2 + d1*h3 = f3:

            aa    = d1/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
            bb    =-d1*(h3*h3-h2*h2)/(h2*h3*(h3-h2)) +                        &
                        f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))

            f(2,1) =         d1
            f(3,1) = 2.0_gpu*bb
            f(4,1) = 6.0_gpu*aa

            f(2,2) = 3.0_gpu*aa*h2*h2 + 2.0_gpu*bb*h2 + d1
            f(3,2) = 6.0_gpu*aa*h2    + 2.0_gpu*bb
            f(4,2) = 6.0_gpu*aa

            f(2,3) = 3.0_gpu*aa*h3*h3 + 2.0_gpu*bb*h3 + d1
            f(3,3) = 6.0_gpu*aa*h3    + 2.0_gpu*bb
            f(4,3) = 6.0_gpu*aa

         else if( ( i_bc1 == 2 ) .or. ( i_bc1 == 4 ) .or. ( i_bc1 == 6 ) )then

!           f'' LHS condition; not-a-knot RHS
!           b1 = f'' LHS BC

            h2     = x(2)-x(1)
            h3     = x(3)-x(1)

            f2     = f(1,2)-f(1,1)
            f3     = f(1,3)-f(1,1)

!           Find cubic aa*xx**3 + (b1/2)*xx**2 + bb*xx
!           satisfying aa*h2**3 + bb*h2 = f2 -(b1/2)*h2**2
!           and aa*h3**3 + bb*h3 = f3 -(b1/2)*h3**2:

            aa     =-(b1/2)*(h3-h2)/(h3*h3-h2*h2) -                           &
                      f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
            bb     =-(b1/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2) +                     &
                      f2*h3*h3/(h2*(h3*h3-h2*h2)) -                           &
                      f3*h2*h2/(h3*(h3*h3-h2*h2))

            f(2,1) =         bb
            f(3,1) =         b1
            f(4,1) = 6.0_gpu*aa

            f(2,2) = 3.0_gpu*aa*h2*h2 + b1*h2 + bb
            f(3,2) = 6.0_gpu*aa*h2    + b1
            f(4,2) = 6.0_gpu*aa

            f(2,3) = 3.0_gpu*aa*h3*h3 + b1*h3 + bb
            f(3,3) = 6.0_gpu*aa*h3    + b1
            f(4,3) = 6.0_gpu*aa

         elseif( ( i_bcn == 1 ) .or. ( i_bcn == 3 ) .or. ( i_bcn == 5 ) )then

!           f' RHS condition; not-a-knot LHS
!           dn = f' RHS BC

            h2     = x(2)-x(3)
            h3     = x(1)-x(3)

            f2     = f(1,2)-f(1,3)
            f3     = f(1,1)-f(1,3)

!           Find cubic aa*xx**3 + bb*xx**2 + dn*xx
!           satisfying aa*h2**3 + bb*h2**2 + dn*h2 = f2
!           and aa*h3**3 + bb*h3**2 + dn*h3 = f3:

            aa     = dn/(h2*h3) + f3/(h3*h3*(h3-h2)) - f2/(h2*h2*(h3-h2))
            bb     =-dn*(h3*h3-h2*h2)/(h2*h3*(h3-h2)) +                       &
                         f2*h3/(h2*h2*(h3-h2)) - f3*h2/(h3*h3*(h3-h2))

            f(2,3) =         dn
            f(3,3) = 2.0_gpu*bb
            f(4,3) = 6.0_gpu*aa

            f(2,2) = 3.0_gpu*aa*h2*h2 + 2.0_gpu*bb*h2 + dn
            f(3,2) = 6.0_gpu*aa*h2    + 2.0_gpu*bb
            f(4,2) = 6.0_gpu*aa

            f(2,1) = 3.0_gpu*aa*h3*h3 + 2.0_gpu*bb*h3 + dn
            f(3,1) = 6.0_gpu*aa*h3    + 2.0_gpu*bb
            f(4,1) = 6.0_gpu*aa

         elseif( ( i_bcn == 2 ) .or. ( i_bcn == 4 ) .or. ( i_bcn == 6 ) )then

!           f'' RHS condition; not-a-knot LHS
!           bn = f'' RHS BC

            h2     = x(2)-x(3)
            h3     = x(1)-x(3)

            f2     = f(1,2)-f(1,3)
            f3     = f(1,1)-f(1,3)

!           Find cubic aa*xx**3 + (bn/2)*xx**2 + bb*xx
!           satisfying aa*h2**3 + bb*h2 = f2 -(bn/2)*h2**2
!           and aa*h3**3 + bb*h3 = f3 -(bn/2)*h3**2:

            aa     =-(bn/2)*(h3-h2)/(h3*h3-h2*h2) -                           &
                      f2/(h2*(h3*h3-h2*h2)) + f3/(h3*(h3*h3-h2*h2))
            bb     =-(bn/2)*h2*h3*(h3-h2)/(h3*h3-h2*h2) +                     &
                      f2*h3*h3/(h2*(h3*h3-h2*h2)) -                           &
                      f3*h2*h2/(h3*(h3*h3-h2*h2))

            f(2,3) =         bb
            f(3,3) =         bn
            f(4,3) = 6.0_gpu*aa

            f(2,2) = 3.0_gpu*aa*h2*h2 + bn*h2 + bb
            f(3,2) = 6.0_gpu*aa*h2    + bn
            f(4,2) = 6.0_gpu*aa

            f(2,1) = 3.0_gpu*aa*h3*h3 + bn*h3 + bb
            f(3,1) = 6.0_gpu*aa*h3    + bn
            f(4,1) = 6.0_gpu*aa

         endif
         
      elseif( n > 2 )then

!        Set up tridiagonal system for A*y=B where y(i) are the second
!        derivatives at the knots.
!
!        f(2,i) are the diagonal elements of A
!        f(4,i) are the off-diagonal elements of A
!        f(3,i) are the B elements/3, and will become c/3 upon solution

         f(4,1) = x(2)-x(1)
         f(3,2) = (f(1,2)-f(1,1))/f(4,1)
         do i = 2, n-1
            f(4,i  ) = x(i+1)-x(i)
            f(2,i  ) = 2.0_gpu* (f(4,i-1)+f(4,i) )
            f(3,i+1) = ( f(1,i+1)-f(1,i) )/f(4,i)
            f(3,i  ) =   f(3,i+1)-f(3,i)
         enddo

!        Save now (dmc):

         elem21  = f(4,  1)
         elemnn1 = f(4,n-1)

!        BC's:

!        Left:

         if( i_bc1 == -1 )then
             f(2,1)    = 2.0_gpu*( f(4,1)+f(4,n-1) )
             f(3,1)    = ( f(1,2)-f(1,1) )/f(4,1)-( f(1,n)-f(1,n-1) )/f(4,n-1)
             wk(1    ) = f(4,n-1)
             wk(2:n-3) = 0.0_gpu
             wk(  n-2) = f(4,n-2)
             wk(  n-1) = f(4,n-1)
        elseif( ( i_bc1 == 1 ) .or. ( i_bc1 == 3 ) .or. ( i_bc1 == 5 ) )then
             f(2,1)    = 2.0_gpu*f(4,1)
             f(3,1)    = ( f(1,2)-f(1,1) )/f(4,1) - d1
        elseif( ( i_bc1 == 2 ) .or. ( i_bc1 == 4 ) .or. ( i_bc1 == 6 ) )then
             f(2,1)    = 2.0_gpu*f(4,1)
             f(3,1)    = f(4,1)*b1/3.0_gpu
             f(4,1)    = 0.0_gpu  ! upper diagonal only (dmc: cf elem21)
        elseif( ( i_bc1 == 7 ) )then
             f(2,1)    =-f(4,1)
             f(3,1)    = f(3,3)/(x(4)-x(2))-f(3,2)/(x(3)-x(1))
             f(3,1)    = f(3,1)*f(4,1)**2/(x(4)-x(1))
        else                             ! not a knot:
             imin      = 2
             f(2,2)    = f(4,1) + 2.0_gpu*f(4,2)
             f(3,2)    = f(3,2)*f(4,2)/( f(4,1)+f(4,2) )
        endif
        

!       Right:
        if(     ( i_bcn == 1 ) .or. ( i_bcn == 3 ) .or. ( i_bcn == 5 ) )then
             f(2,n)    = 2.0_gpu*f(4,n-1)
             f(3,n)    =-( f(1,n)-f(1,n-1) )/f(4,n-1) + dn
        elseif( ( i_bcn == 2 ) .or. ( i_bcn == 4 ) .or. ( i_bcn == 6 ) )then
             f(2,n)    = 2.0_gpu*f(4,n-1)
             f(3,n)    = f(4,n-1)*bn/3.0_gpu
!            f(4,n-1)  = 0.0_gpu  ! dmc: preserve f(4,n-1) for back subst.
             elemnn1   = 0.0_gpu  ! lower diagonal only (dmc)
        elseif( ( i_bcn == 7 ) )then
             f(2,n)    =-f(4,n-1)
             f(3,n)    = f(3,n-1)/( x(n)-x(n-2) )-f(3,n-2)/( x(n-1)-x(n-3) )
             f(3,n)    =-f(3,n)*f(4,n-1)**2/( x(n)-x(n-3) )
        elseif( ( i_bc1 /=-1 ) )then  ! not a knot:
             imax      = n-1
              f(2,n-1) = 2.0_gpu*f(4,n-2)+f(4,n-1)
              f(3,n-1) = f(3,n-1)*f(4,n-2)/( f(4,n-1)+f(4,n-2) )
        endif

        if( i_bc1 == -1 )then

!           Solve system of equations for second derivatives at the knots.
!           Periodic BC. Forward elimination:

            do i = 2, n-2
               t        = f(4,i-1)/f(2,i-1)
               f(2,i)   = f(2,i)-t*f(4,i-1)
               f(3,i)   = f(3,i)-t*f(3,i-1)
               wk(i)    = wk(i)-t*wk(i-1)
               q        = wk(n-1)/f(2,i-1)
               wk(n-1)  =-q*f(4,i-1)
               f(2,n-1) = f(2,n-1)-q*wk(i-1)
               f(3,n-1) = f(3,n-1)-q*f(3,i-1)
            enddo
          
!           Correct the n-1 element:

            wk(n-1)  = wk(n-1) + f(4,n-2)
            
!           Complete the forward elimination:
!           wk(n-1) and wk(n-2) are the off-diag elements of the lower corner.

            t        = wk(n-1)/f(2,n-2)
            f(2,n-1) = f(2,n-1)-t*wk(n-2)
            f(3,n-1) = f(3,n-1)-t*f(3,n-2)
            
!           Back substitution:

            f(3,n-1) = f(3,n-1)/f(2,n-1)
            f(3,n-2) =(f(3,n-2)-wk(n-2)*f(3,n-1))/f(2,n-2)

            do ib = 3, n-1
               i      = n-ib
               f(3,i) =(f(3,i)-f(4,i)*f(3,i+1)-wk(i)*f(3,n-1))/f(2,i)
            enddo
            
            f(3,n)   = f(3,1)
            
        else

!          Non-periodic BC, Forward elimination.
!          For Not-A-Knot BC the off-diagonal end elements are not equal.

           do i = imin+1, imax

              if( (i == n-1 )   .and. ( imax == n-1 ) )then
                   t = ( f(4,i-1)-f(4,i) )/f(2,i-1)
              else
                   if(     i == 2)then
                      t = elem21 /f(2,i-1)
                   elseif( i == n )then
                      t = elemnn1/f(2,i-1)
                   else
                      t = f(4,i-1)/f(2,i-1)
                   endif
              endif
             
              if( ( i == imin+1 ) .and. ( imin == 2 ) )then
                    f(2,i) = f(2,i)-t*( f(4,i-1)-f(4,i-2) )
              else
                    f(2,i) = f(2,i)-t*f(4,i-1)
              endif

              f(3,i) = f(3,i)-t*f(3,i-1)

           enddo
          
!          Back substitution:

           f(3,imax) = f(3,imax)/f(2,imax)
          
           do ib =1, imax-imin
              i = imax-ib
              if( ( i == 2 ) .and. ( imin == 2 ) )then
                    f(3,i) = ( f(3,i)-(f(4,i)-f(4,i-1) )*f(3,i+1))/f(2,i)
              else
                    f(3,i) = ( f(3,i)- f(4,i)*f(3,i+1) )/f(2,i)
              endif
           enddo
          
!          Reset d array to step size:

           f(4,1)   = x(2)-x(1)
           f(4,n-1) = x(n)-x(n-1)

!          Set f(3,1) for not-a-knot:

           if( i_bc1 <= 0 .or. i_bc1 > 7 )then
             f(3,1) = (f(3,2)*( f(4,1)+f(4,2) )-f(3,3)*f(4,1))/f(4,2)
           endif
          
!          Set f(3,n) for not-a-knot:

           if( i_bcn <= 0 .or. i_bcn > 7 )then
             f(3,n) = f(3,n-1)+( f(3,n-1)-f(3,n-2) )*f(4,n-1)/f(4,n-2)
           endif

        endif
        
!       f(3,i) is now the sigma(i) of the text and f(4,i) is the step size.
!       Compute polynomial coefficients:

        do i = 1, n-1
          f(2,i) = ( f(1,i+1)-f(1,i) )/f(4,i)-f(4,i)*(f(3,i+1) + &
                     2.0_gpu*f(3,i))
          f(4,i) = ( f(3,i+1)-f(3,i) )/f(4,i)
          f(3,i) =   6.0_gpu*f(3,i)
          f(4,i) =   6.0_gpu*f(4,i)
        enddo

        if( i_bc1 == -1 )then
            f(2,n) = f(2,1)
            f(3,n) = f(3,1)
            f(4,n) = f(4,1)
        else
            hn     = x(n)-x(n-1)
            f(2,n) = f(2,n-1) + hn*(f(3,n-1)+0.5_gpu*hn*f(4,n-1))
            f(3,n) = f(3,n-1) + hn*f(4,n-1)
            f(4,n) = f(4,n-1)
            if(     ( i_bcn == 1 ) .or. ( i_bcn == 3 ) .or.                   &
                    ( i_bcn == 5 ) )then
                      f(2,n) = dn
            elseif( ( i_bcn == 2 ) .or. ( i_bcn == 4 ) .or.                   &
                    ( i_bcn == 6 ) )then
                      f(3,n) = bn
            endif
        endif
      endif
      
end subroutine v_spline

