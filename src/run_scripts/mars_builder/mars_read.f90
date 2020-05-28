#ifndef NC
# define NC 3
#endif
#ifndef TOKAMAK
# define TOKAMAK 1
#endif

subroutine B_F( B_Re, B_Im, SH, NT, M1, M2, n, s, chi, phi, B, itype, iflag )
!------------------------------------------------------------------------------
!
! *** Description
!
!     Calculate MARS-F contravariant field B1, B2 or B3.
!
! *** Calling variables:
!
!     B_Re      GPU             i/p     : Re{B} spline coefficients         [m]
!     B_Im      GPU             i/p     : Im{B} spline coefficients         [m]
!     SH        GPU             i/p     : s grid matching above             [-]
!     NT        I               i/p     : # radial points                   [-]
!     M1        I               i/p     : Poloidal harmonic cut-off (lower) [-]
!     M2        I               i/p     : Poloidal harmonic cut-off (upper) [-]
!     n         I               i/p     : Toroidal mode number of pertbn.   [-]
!     s         GPU             i/p     : Requested s                       [-]
!     chi       GPU             i/p     : Requested chi                     [-]
!     phi       GPU             i/p     : Requested phi                     [-]
!     B         GPU             o/p     : [Re{B},Im{B}] output           [T.m2]
!     itype     I               i/p     : 0=default, 1=more info            [-]
!     iflag     I               o/p     : Error flag                        [-]
!
! *** Author:
!
!     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
!
! *** Error flags:               
!
!     iflag = 0         : Return status OK.
!             1         : ERROR : evspline_ failure @ point 1
!             2         : ERROR :            -"-            2
!
! *** Creation Date:
!
!     12/05/2016
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,       intent(in)    :: NT, M1, M2, n
      real( gpu ),   intent(in)    :: s, chi, phi
      real( gpu ),   intent(in)    :: B_Re( 2, NT, M2-M1+1 )
      real( gpu ),   intent(in)    :: B_Im( 2, NT, M2-M1+1 )
      real( gpu ),   intent(in)    :: SH(NT)
      real( gpu ),   intent(out)   :: B(18)
      integer,       intent(in)    :: itype
      integer,       intent(out)   :: iflag

      integer                      :: M, ict(3), ierr
      real( gpu )                  :: sin_, cos_, f1(3), f2(3), dat(4), M_

!     Initialize effor flag:

      if( itype<=0 )then

!         Request value and 1st derivative:
      
          ict = [1,0,0]
      else

!         Request value and 1st derivative:
      
          ict = [1,1,1]
      endif

      iflag  = 0
      
!     Initialize fields:

      B(1:18) = 0.0_gpu

!     Sum over harmonics:

      do M=M1,M2
         sin_ = sin( real(M,gpu)*chi + real(n,gpu)*phi )
         cos_ = cos( real(M,gpu)*chi + real(n,gpu)*phi )

         call evspline_( s, SH, NT, 2, B_Re(1:2, 1:NT, M-M1+1), ict, f1, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':B_F   : ERROR : evspline_ failure @ point 1'
             call display ( bold, lun=io(1) )
             iflag = 1
             return
         endif
         
         call evspline_( s, SH, NT, 2, B_Im(1:2, 1:NT, M-M1+1), ict, f2, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':B_F   : ERROR : evspline_ failure @ point 2'
             call display ( bold, lun=io(1) )
             iflag = 2
             return
         endif

!        B:

         B(1) = B(1) + ( f1(1)*cos_ - f2(1)*sin_ )
         B(2) = B(2) + ( f1(1)*sin_ + f2(1)*cos_ )

         if( itype==-1 )then
             write(6,'(A16,I3,A26,2E16.8)') ':mars_read : m:', M,             &
                     ' B1(partial,contribution):', B(1), f1(1)
         endif

         if( itype==1 )then

             M_   = real(M**2,gpu)

!            dB/ds:

             B( 3) = B( 3) + ( f1(2)*cos_ - f2(2)*sin_ )
             B( 4) = B( 4) + ( f1(2)*sin_ + f2(2)*cos_ )

!            d2B/ds2:

             B( 5) = B( 5) + ( f1(3)*cos_ - f2(3)*sin_ )
             B( 6) = B( 6) + ( f1(3)*sin_ + f2(3)*cos_ )

!            dB/dchi:

             B( 7) = B( 7) - ( f1(1)*sin_ + f2(1)*cos_ )*real(M,gpu)
             B( 8) = B( 8) + ( f1(1)*cos_ - f2(1)*sin_ )*real(M,gpu)

!            d2B/dsdchi:

             B( 9) = B( 9) - ( f1(2)*sin_ + f2(2)*cos_ )*real(M,gpu)
             B(10) = B(10) + ( f1(2)*cos_ - f2(2)*sin_ )*real(M,gpu)

!            d3B/ds2dchi:

             B(11) = B(11) - ( f1(3)*sin_ + f2(3)*cos_ )*real(M,gpu)
             B(12) = B(12) + ( f1(3)*cos_ - f2(3)*sin_ )*real(M,gpu)

!            d2B/dchi2:

             B(13) = B(13) - ( f1(1)*cos_ - f2(1)*sin_ )*M_
             B(14) = B(14) - ( f1(1)*sin_ + f2(1)*cos_ )*M_

!            d3B/dsdchi2:

             B(15) = B(15) - ( f1(2)*cos_ - f2(2)*sin_ )*M_
             B(16) = B(16) - ( f1(2)*sin_ + f2(2)*cos_ )*M_

!            d4B/ds2dchi2:

             B(17) = B(17) - ( f1(3)*cos_ - f2(3)*sin_ )*M_
             B(18) = B(18) - ( f1(3)*sin_ + f2(3)*cos_ )*M_

         endif
      enddo

end subroutine B_F

subroutine rspos( R_Re, R_Im, Z_Re, Z_Im, SH, NT, M0, s, chi, R, Z, itype,    &
                  iflag )
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Calculate, R, dR/ds, dR/dchi, Z, dZ/ds and dZ/dchi
!
! *** Calling variables:
!
!     R_Re      GPU             i/p     : Re{R} spline coefficients         [m]
!     R_Im      GPU             i/p     : Im{R} spline coefficients         [m]
!     Z_Re      GPU             i/p     : Re{Z} spline coefficients         [m]
!     Z_Im      GPU             i/p     : Im{Z} spline coefficients         [m]
!     SH        GPU             i/p     : s grid matching above             [-]
!     NT        I               i/p     : # radial points                   [-]
!     M0        I               i/p     : Poloidal harmonic cut-off         [-]
!     s         GPU             i/p     : Requested s                       [-]
!     chi       GPU             i/p     : Requested chi                     [-]
!     R         GPU             o/p     : [R,dR/ds,dR/dchi] output          [m]
!     Z         GPU             o/p     : [Z,dZ/ds,dZ/dchi] output          [m]
!     itype     I               i/p     : 0=default, 1=more info            [-]
!     iflag     I               o/p     : Error flag                        [-]
!
! *** Author:
!
!     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
!
! *** Error flags:               
!
!     iflag = 0         : Return status OK.
!             1         : ERROR : evspline_ failure @ point 1
!             2         : ERROR :            -"-            2
!             3         : ERROR :            -"-            3
!             4         : ERROR :            -"-            4
!             5         : ERROR :            -"-            5
!             6         : ERROR :            -"-            6
!             7         : ERROR :            -"-            7
!             8         : ERROR :            -"-            8
!             9         : ERROR :            -"-            9
!             10        : ERROR :            -"-            10
!
! *** Creation Date:
!
!     12/05/2016
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      integer,       intent(in)   :: NT, M0
      real( gpu ),   intent(in)   :: s, chi
      real( gpu ),   intent(in)   :: R_Re( 2, NT, M0 ), R_Im( 2, NT, M0 )
      real( gpu ),   intent(in)   :: Z_Re( 2, NT, M0 ), Z_Im( 2, NT, M0 )
      real( gpu ),   intent(in)   :: SH  (NT)
      real( gpu ),   intent(out)  :: R(15), Z(15)
      integer,       intent(in)   :: itype
      integer,       intent(out)  :: iflag
      
      integer                     :: M, ict(3), ierr
      real( gpu )                 :: sin_, cos_, f1(3), f2(3), dat(4)
      real( gpu )                 :: M_1, M_2, M_3

!     Initialize effor flag:

      iflag = 0
      
!     Initialize dR/dchi & dZ/dchi:

      R(3:15)  = 0.0_gpu
      Z(3:15)  = 0.0_gpu

!     Request value and 1st derivative:
      
      if( itype==0 )then
          ict = [1,1,0]
      else
          ict = [1,1,1]
      endif

!     m=0 part:
      
      call evspline_( s, SH, NT, 2, R_Re( 1:2, 1:NT, 1 ), ict, f1, ierr )
      
      if( ierr /= 0 )then
          call display ( bold, lun=io(1) )
               write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 1'
          call display ( bold, lun=io(1) )
          iflag = 1
          return
      endif

      R(1:2) = f1(1:2)

      call evspline_( s, SH, NT, 2, Z_Re( 1:2, 1:NT, 1 ), ict, f1, ierr )

      if( ierr /= 0 )then
          call display ( bold, lun=io(1) )
               write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 2'
          call display ( bold, lun=io(1) )
          iflag = 2
          return
      endif

      Z(1:2) = f1(1:2)

!     Add in the higher m harmonics:

      do M=2, M0

         M_1  = real(M-1,gpu)

         sin_ = sin( M_1*chi )
         cos_ = cos( M_1*chi )

         call evspline_( s, SH, NT, 2, R_Re( 1:2, 1:NT, M ), ict, f1, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 3'
             call display ( bold, lun=io(1) )
             iflag = 3
             return
         endif
         
         call evspline_( s, SH, NT, 2, R_Im( 1:2, 1:NT, M ), ict, f2, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 4'
             call display ( bold, lun=io(1) )
             iflag = 4
             return
         endif

!        R, dR/ds & dR/dchi:

         R(1) = R(1) + 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )
         R(2) = R(2) + 2.0_gpu*( f1(2)*cos_ - f2(2)*sin_ )
         R(3) = R(3) - 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_1

         if( itype==1 )then

             M_2 = real((M-1)**2,gpu)
             M_3 = real((M-1)**3,gpu)

             R( 4) = R( 4) + 2.0_gpu*( f1(3)*cos_ - f2(3)*sin_ )
             R( 6) = R( 6) - 2.0_gpu*( f1(2)*sin_ + f2(2)*cos_ )*M_1
             R( 7) = R( 7) - 2.0_gpu*( f1(3)*sin_ + f2(3)*cos_ )*M_1
             R( 9) = R( 9) - 2.0_gpu*( f1(2)*cos_ - f2(2)*sin_ )*M_2
             R(10) = R(10) - 2.0_gpu*( f1(3)*cos_ - f2(3)*sin_ )*M_2
             R(12) = R(12) - 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )*M_2
             R(13) = R(13) + 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_3
             R(14) = R(14) + 2.0_gpu*( f1(2)*sin_ + f2(2)*cos_ )*M_3
             R(15) = R(15) + 2.0_gpu*( f1(3)*sin_ + f2(3)*cos_ )*M_3

!            3rd derivative needs an extra call:

             ict = [3,0,0]

             call evspline_(s, SH, NT, 2, R_Re( 1:2, 1:NT, M ), ict, f1, ierr)

             if( ierr /= 0 )then
                 call display ( bold, lun=io(1) )
                      write(io(1),*)                                          &
                      ':rspos : ERROR : evspline_ failure @ point 5'
                 call display ( bold, lun=io(1) )
                 iflag = 5
                 return
             endif
         
             call evspline_(s, SH, NT, 2, R_Im( 1:2, 1:NT, M ), ict, f2, ierr)

             if( ierr /= 0 )then
                 call display ( bold, lun=io(1) )
                      write(io(1),*)                                          &
                      ':rspos : ERROR : evspline_ failure @ point 6'
                 call display ( bold, lun=io(1) )
                 iflag = 6
                 return
             endif

             R( 5) = R( 5) + 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )
             R( 8) = R( 8) - 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_1
             R(11) = R(11) - 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )*M_2

             ict = [1,1,1]

         endif

         call evspline_( s, SH, NT, 2, Z_Re( 1:2, 1:NT, M ), ict, f1, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 7'
             call display ( bold, lun=io(1) )
             iflag = 7
             return
         endif

         call evspline_( s, SH, NT, 2, Z_Im( 1:2, 1:NT, M ), ict, f2, ierr )

         if( ierr /= 0 )then
             call display ( bold, lun=io(1) )
                  write(io(1),*) ':rspos : ERROR : evspline_ failure @ point 8'
             call display ( bold, lun=io(1) )
             iflag = 8
             return
         endif

!        Z, dZ/ds & dZ/dchi:

         Z(1) = Z(1) + 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )
         Z(2) = Z(2) + 2.0_gpu*( f1(2)*cos_ - f2(2)*sin_ )
         Z(3) = Z(3) - 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_1

         if( itype==1 )then

             Z( 4) = Z( 4) + 2.0_gpu*( f1(3)*cos_ - f2(3)*sin_ )
             Z( 6) = Z( 6) - 2.0_gpu*( f1(2)*sin_ + f2(2)*cos_ )*M_1
             Z( 7) = Z( 7) - 2.0_gpu*( f1(3)*sin_ + f2(3)*cos_ )*M_1
             Z( 9) = Z( 9) - 2.0_gpu*( f1(2)*cos_ - f2(2)*sin_ )*M_2
             Z(10) = Z(10) - 2.0_gpu*( f1(3)*cos_ - f2(3)*sin_ )*M_2
             Z(12) = Z(12) - 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )*M_2
             Z(13) = Z(13) + 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_3
             Z(14) = Z(14) + 2.0_gpu*( f1(2)*sin_ + f2(2)*cos_ )*M_3
             Z(15) = Z(15) + 2.0_gpu*( f1(3)*sin_ + f2(3)*cos_ )*M_3

!            3rd derivative needs an extra call:

             ict = [3,0,0]

             call evspline_(s, SH, NT, 2, Z_Re( 1:2, 1:NT, M ), ict, f1, ierr)

             if( ierr /= 0 )then
                 call display ( bold, lun=io(1) )
                      write(io(1),*) &
                      ':rspos : ERROR : evspline_ failure @ point 9'
                 call display ( bold, lun=io(1) )
                 iflag = 9
                 return
             endif
         
             call evspline_(s, SH, NT, 2, Z_Im( 1:2, 1:NT, M ), ict, f2, ierr)

             if( ierr /= 0 )then
                 call display ( bold, lun=io(1) )
                      write(io(1),*) &
                      ':rspos : ERROR : evspline_ failure @ point 10'
                 call display ( bold, lun=io(1) )
                 iflag = 10
                 return
             endif

             Z( 5) = Z( 5) + 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )
             Z( 8) = Z( 8) - 2.0_gpu*( f1(1)*sin_ + f2(1)*cos_ )*M_1
             Z(11) = Z(11) - 2.0_gpu*( f1(1)*cos_ - f2(1)*sin_ )*M_2

             ict = [1,1,1]

         endif

      enddo

end subroutine rspos

program mars_read
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Build [R,Z] BPLASMA data from Fourier Harmonic BPLASMA format (from
!
!     "MARS-F computed raw field data for mapping", Yueqiang Liu,
!     March 20, 2016.
!
! *** Notes
!
!     The MARS-F field is given by (e.g.):
!
!     BR = Re{BR}*cos( -|n|*PHI_M ) - Im{BR}*sin( -|n|*PHI_M ) [1]
!
!     where n is -ve ( e.g. n=-3 ) and PHI_M is in MARS-F coordinates
!     (i.e. PHI running clockwise looking down on the machine, i.e. the
!     [R,Z,PHI] right handed coordinate system). ITER coordinates are
!     for the [R,PHI,Z] right handed coordinate system, i.e. PHI_I runs in the
!     opposite direction counter clockwise.
!
!     The ITER Filaments model field is for +ve n, e.g. n=+3.
!
!     Eqn. [1] can be recast as follows:
!
!     BR = Re{BR}*cos( +|n|*PHI_I ) - Im{BR}*sin( +|n|*PHI_I ) [2]
!
!     because PHI_I = -PHI_M. In other words, -|n| in the MARS-F
!     coordinate system is the same field as +|n| in the ITER coordinate 
!     system - almost....
!
!     The returned field will still be in MARS-F coordinates, i.e.
!     BT will be pointing in the wrong direction and will have to be flipped.
!
!     Eqn. [2] is effectively then for a +|n| field in the ITER coordinate
!     system but with BT the wrong sign.
!
!     The phase shift [86,0,34] etc. is applied once in the ITER coordinate
!     system.  Additional phase shifts due to the offset coils depend
!     upon both the primary mode number n0 and the harmonic mode number
!     n.

! *** Compiler options:
!
!     -DOLD      : Load in MARS-F output with all coil sets included.
!     -DFXXYY    : Calculate and append fxx, fyy and fxxyy.
!     -DSCHIPLOT : Output some data for plotting.
!
! *** Calling variables:
!
! *** Author:
!
!     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323 : original
!     Samuel Ward, 81/147, ITER organization, samuel.ward@iter.org, +33 4 42 17 81 59 : LOCUST_IO workflow additions
!
! *** Error flags:               
!
!     iflag = 0         : Return status OK.
!             1         : ERROR : 
!
! *** Creation Date:
!
!     12/05/2016
!
!------------------------------------------------------------------------------

      use prec_mod
      use disp_mod

      implicit none

      character( len=1000 )       :: root = '/home/ITER/wards2/' !prepend this to all file writes

#if (TOKAMAK==1)
#if defined (MATCH)
      character( len=1006  ),                                                    &
                    dimension(3)  :: TAIL = ['_U_VAC','_M_VAC','_L_VAC']
#else
#if defined (VAC)
      character( len=1006  ),                                                    &
                    dimension(3)  :: TAIL = ['_U_VAC','_M_VAC','_L_VAC']
#elif defined (PLS)
      character( len=1006  ),                                                    &
                    dimension(3)  :: TAIL = ['_U_PLS','_M_PLS','_L_PLS']
#else
      character( len=1006  ),                                                    &
                    dimension(3)  :: TAIL = ['_U_VAC','_M_VAC','_L_VAC']
#endif
#endif
#else
      character( len=10011 ),                                                    &
                    dimension(2)  :: TAIL = ['_upper','_lower']
#endif

#if defined (MATCH)

  #if (TOKAMAK==1)

    #if   (MATCH==1)
      character( len=10031 )         :: mtch = 'data/&
                                            &15MA_ELM_upper.txt_cleaned'
    #elif (MATCH==2)
      character( len=10032 )         :: mtch = 'data/&
                                            &15MA_ELM_middle.txt_cleaned'
    #elif (MATCH==3)
      character( len=10031 )         :: mtch = 'data/&
                                            &15MA_ELM_lower.txt_cleaned'
    #else
      stop ! This will throw a compiler error.
    #endif

  #else

    #if   (MATCH==1)
      character( len=1001  )         :: mtch = '?'
    #elif (MATCH==2)
      character( len=1001  )         :: mtch = '?'
    #else
      stop ! This will throw a compiler error.
    #endif

  #endif

#endif

#if defined (UMP3D)
  #if (TOKAMAK==1)
      character( len=10031 )         :: mtch = 'data/&
                                            &15MA_ELM_coils.txt_cleaned'
  #else
      character( len=1001  )         :: mtch = '?'
  #endif
#endif

#ifndef OLD

  #if (TOKAMAK==1)

    #ifndef ORI
      #if (NC==3)
      character( len=10057 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_MOD_EQ/BPLASMA_MARSF_MOD_n3'
      #elif (NC==6)
      character( len=10057 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_MOD_EQ/BPLASMA_MARSF_MOD_n6'
      #elif (NC==12)
      character( len=10058 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_MOD_EQ/BPLASMA_MARSF_MOD_n12'
      #elif (NC==15)
      character( len=10058 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_MOD_EQ/BPLASMA_MARSF_MOD_n15'
      #else
      character( len=10057 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_MOD_EQ/BPLASMA_MARSF_MOD_n3'
      #endif
    #else
      #if (NC==3)
      character( len=10053 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_ORI_EQ/BPLASMA_MARSF_n3'
      #elif (NC==6)
      character( len=10053 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_ORI_EQ/BPLASMA_MARSF_n6'
      #elif (NC==12)
      character( len=10054 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_ORI_EQ/BPLASMA_MARSF_n12'
      #elif (NC==15)
      character( len=10054 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_ORI_EQ/BPLASMA_MARSF_n15'
      #else
      character( len=10053 )         :: file = 'data/ITER_15MA_10470/&
                                     &ELM_COIL_ORI_EQ/BPLASMA_MARSF_n3'
      #endif
    #endif

  #else

    #if (NC==2)
      character( len=10040 )         :: file = 'data/33143_2730/&
                                     &BPLASMA_newformat_vac'
    #else
       stop
    #endif

  #endif

#else

  #if (TOKAMAK==1)
      character( len=10040 )         :: file = 'data/&
                                     &BPLASMA_MARSF_VAC_n3_3COILS_OLD.txt'
  #else
      character( len=?  )         :: file = '?'
  #endif

#endif

      character( len=100100 )        :: STR
      character( len=1001 )          :: STR_C
      character( len=1002 )          :: STR_N
      character( len=1004 )          :: STR_R, STR_Z
      character( len=1001 ),                                                     &
                     dimension(3) :: TAG_I = ['U','M','L']
      integer,        parameter   :: lun   = 5
      integer                     :: ios
      integer                     :: idat(7)
      integer                     :: nmde, M0, M1, M2, NR, NV, NT, ICOIL
      real( gpu )                 :: ICURR, IN_IC
      integer                     :: M, ilinx, iflag, ict(3), imsk(3)
      integer                     :: ierr
      integer                     :: NT_in, NT_out
      integer                     :: N1_in, N1_out
      integer                     :: N2_in, N2_out
      integer                     :: N3_in, N3_out
      integer                     :: iR, iZ, iR_, iZ_, j, k, iter, jR, jZ
      integer( long )             :: i
#if defined (FXXYY)
      integer,        parameter   :: jR0=1
      integer,        parameter   :: jR1=7
      integer,        parameter   :: jZ0=1
      integer,        parameter   :: jZ1=7
#else
      integer,        parameter   :: jR0=4
      integer,        parameter   :: jR1=4
      integer,        parameter   :: jZ0=4
      integer,        parameter   :: jZ1=4 
#endif
      real( gpu ),    parameter               :: SMATCH= 0.5_gpu

#ifndef OLD
  #if (TOKAMAK==1 )
      real( gpu ),    parameter, dimension(3) :: IKATN = [45.0_gpu,           &
                                                          45.0_gpu,           &
                                                          45.0_gpu]
  #else

!     No scaling for AUG simulations:

      real( gpu ),    parameter, dimension(2) :: IKATN = [1.0_gpu,            &
                                                          1.0_gpu] 
  #endif

!     MARS-F data does not necessarily already contain the n*30 and n*26.7 phase shifts:

  #if (TOKAMAK==1)

    #if( NC==3 )
      real( gpu ),    parameter, dimension(3) :: PH0   = [ 0.0_gpu,           &
                                                           0.0_gpu,           &
                                                           0.0_gpu]
    #else
      real( gpu ),    parameter, dimension(3) :: PH0   = [ 30.0_gpu,          &
                                                           26.7_gpu,          &
                                                           30.0_gpu]
    #endif

      real( gpu ),    parameter, dimension(3) :: dPH   = [ 29.4_gpu,          &
                                                           20.9_gpu,          &
                                                           30.5_gpu]
  #else
      real( gpu ),    parameter, dimension(2) :: PH0   = [ 0.0_gpu,           &
                                                           0.0_gpu]

      real( gpu ),    parameter, dimension(3) :: dPH   = [ 0.0_gpu,          &
                                                           0.0_gpu]          &
  #endif

  #if (TOKAMAK==1)

!     ITER phase adjustment settings:

      real( gpu ),               dimension(3) :: PH1   = 0.0_gpu
    #if (NC==3)
      real( gpu ),               dimension(3) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0,                         &
                                                 1.0_gpu ]
    #elif (NC==6)
      real( gpu ),               dimension(3) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0,                         &
                                                 1.0_gpu ]
    #elif (NC==12)
      real( gpu ),               dimension(3) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0,                         &
                                                 1.0_gpu ]
    #elif (NC==15)
      real( gpu ),               dimension(3) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0,                         &
                                                 1.0_gpu ]
    #else
      real( gpu ),               dimension(3) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0_gpu,                     &
                                                 1.0_gpu ]
    #endif

  #else

      real( gpu ),               dimension(3) :: PH1   = 0.0_gpu
    #if (NC==2)
      real( gpu ),               dimension(2) :: IMUL  =                      &
                                               [ 1.0_gpu,                     &
                                                 1.0_gpu ]
    #endif

  #endif

#else
  #if (TOKAMAK==1)
      real( gpu ),    parameter, dimension(1) :: IKATN = [45.0_gpu]
  #else
      real( gpu ),    parameter, dimension(1) :: IKATN = [ 1.0_gpu]
  #endif
      real( gpu ),    parameter, dimension(1) :: PH0   = [ 0.0_gpu]
      real( gpu ),    parameter, dimension(1) :: PH1   = [ 0.0_gpu]
      real( gpu ),    parameter, dimension(1) :: dPh   = [ 0.0_gpu]
      real( gpu ),               dimension(1) :: IMUL  = [ 1.0_gpu]    
#endif

#if (TOKAMAK==1)
      integer,        parameter               :: N_coils = 9
#else
      integer,        parameter               :: N_coils = 8
#endif       

#ifndef RESSCAN
#if defined (FINE)
      real( gpu )                 :: dXR   = 0.002_gpu
      real( gpu )                 :: dXZ   = 0.002_gpu
#elif defined (PRODN)
      real( gpu )                 :: dXR   = 0.0050_gpu
      real( gpu )                 :: dXZ   = 0.0050_gpu
#elif defined (COARSE)
      real( gpu )                 :: dXR   = 0.005_gpu
      real( gpu )                 :: dXZ   = 0.005_gpu
#elif defined (BPLASMAOLD)
      real( gpu )                 :: dXR   = 0.0283502_gpu
      real( gpu )                 :: dXZ   = 0.0259352_gpu
#elif defined (RYAN)
      real( gpu )                 :: dXR   = 2.0_gpu/399.0_gpu
      real( gpu )                 :: dXZ   = 3.0_gpu/599.0_gpu
#elif defined (ASCOTAUG)
      real( gpu )                 :: dXR   = 0.0300_gpu
      real( gpu )                 :: dXZ   = 0.0235_gpu
#else
      real( gpu )                 :: dXR   = 0.02_gpu
      real( gpu )                 :: dXZ   = 0.02_gpu
#endif
#else
      real( gpu )                 :: dXR   = RESSCAN
      real( gpu )                 :: dXZ   = RESSCAN
#endif

      integer                     :: nR_f, nR_f_
      integer                     :: nZ_f, nZ_f_
      integer                     :: nT_f  = 0
      integer                     :: iRF
      integer                     :: iZF
      integer,        parameter   :: NMAX  = 20
      integer( long ),parameter   :: NS    = 10000000
      integer( long ),parameter   :: NL    = 3601
      integer                     :: irbox
      integer                     :: izbox
      integer                     :: nR_, nZ_, nT_, i0, j0
      real( gpu )                 :: Rmin_, Rmax_, Zmin_, Zmax_, Tmin_, Tmax_
      real( gpu ),    allocatable :: BMTCH(:,:,:,:), Re_F(:,:,:), Im_F(:,:,:)
      integer,        external    :: omp_get_thread_num
      real( gpu ),    parameter   :: DPTOL = 1.0e-14_gpu
      real( gpu ),    allocatable :: ddat(:), S(:), SM(:), SH(:), q(:)
      real( gpu ),    allocatable :: RF_Re (:,:), RF_Im (:,:)
      real( gpu ),    allocatable :: ZF_Re (:,:), ZF_Im (:,:)
      real( gpu ),    allocatable :: RFM_Re(:,:), RFM_Im(:,:)
      real( gpu ),    allocatable :: ZFM_Re(:,:), ZFM_Im(:,:)
      real( gpu ),    allocatable :: RFH_Re(:,:), RFH_Im(:,:)
      real( gpu ),    allocatable :: ZFH_Re(:,:), ZFH_Im(:,:)
      real( gpu ),    allocatable :: B1_Re (:,:), B1_Im (:,:)
      real( gpu ),    allocatable :: B2_Re (:,:), B2_Im (:,:)
      real( gpu ),    allocatable :: B3_Re (:,:), B3_Im (:,:)
      real( gpu ),    allocatable :: RFH_Re_fspl_in (:,:,:)
      real( gpu ),    allocatable :: RFH_Im_fspl_in (:,:,:)
      real( gpu ),    allocatable :: ZFH_Re_fspl_in (:,:,:)
      real( gpu ),    allocatable :: ZFH_Im_fspl_in (:,:,:)
      real( gpu ),    allocatable :: RFH_Re_fspl_out(:,:,:)
      real( gpu ),    allocatable :: RFH_Im_fspl_out(:,:,:)
      real( gpu ),    allocatable :: ZFH_Re_fspl_out(:,:,:)
      real( gpu ),    allocatable :: ZFH_Im_fspl_out(:,:,:)
      real( gpu ),    allocatable :: B1_Re_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B1_Im_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B2_Re_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B2_Im_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B3_Re_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B3_Im_fspl_in  (:,:,:)
      real( gpu ),    allocatable :: B1_Re_fspl_out (:,:,:)
      real( gpu ),    allocatable :: B1_Im_fspl_out (:,:,:)
      real( gpu ),    allocatable :: B2_Re_fspl_out (:,:,:)
      real( gpu ),    allocatable :: B2_Im_fspl_out (:,:,:)
      real( gpu ),    allocatable :: B3_Re_fspl_out (:,:,:)
      real( gpu ),    allocatable :: B3_Im_fspl_out (:,:,:)
      real( gpu ),    allocatable :: fspl  (:,:), cells(:,:,:)
      real( gpu ),    allocatable :: cells_h(:,:,:)
      real( gpu ),    allocatable :: B(:,:,:,:,:)
      real( gpu ),    allocatable :: BR0(:,:,:), BZ0(:,:,:), BT0(:,:,:)
      real( gpu )                 :: f(3), R(15), Z(15), tht_
      real( gpu )                 :: B1(18), B2(18), B3(18), X_(8)
      real( gpu )                 :: R_p, Z_p, d_p, det, ds, dchi, JAC, sin_, cos_
      real( gpu )                 :: R_mag_, Z_mag_, PHSE(3), DOTP(3), BF(4)
      real( gpu ),    allocatable :: RAN (:,:)
      real( gpu ),    allocatable :: RAN_(:), rho_2(:), tht(:)
      real( gpu ),    allocatable :: PHI_(:)
      real( gpu )                 :: DUM(1001,361,2)
      real( gpu )                 :: s_, chi_, cp, sp, Rmag, Zmag, dat(4)
      real( gpu )                 :: dR_, dZ_, R0F_, R1F_, Z0F_, Z1F_
      real( gpu )                 :: R0C_, R1C_, Z0C_, Z1C_
      real( gpu ),    allocatable :: B_(:,:), RAX(:), ZAX(:), FDER(:,:,:)
      real( gpu ),    allocatable :: fxx(:,:), fyy(:,:), fxxyy(:,:)
      real( gpu ),      parameter :: dROFF=1.0e-04_gpu
      real( gpu ),      parameter :: dZOFF=1.0e-04_gpu
      real( double ),   external  :: omp_get_wtime

      write(io(1),*) ':mars_read : Request toroidal mode |n| = ', NC
      write(io(1),*) ':mars_read : R resolution [m]          = ', dXR
      write(io(1),*) ':mars_read : Z    -"-     [m]          = ', dXZ
#if defined (RESSCAN)
      write(io(1),*) ':mars_read : nR                        = ', &
                       nint((R1F-R0F)/dXR) + 1
      write(io(1),*) ':mars_read : nZ                        = ', &
                       nint((Z1F-Z0F)/dXZ) + 1
#endif

#if (TOKAMAK==1)

  #if defined (UPHASE)
      write(io(1),*) ':mars_read : U Phase set to  : ', UPHASE
      PH1(1) = real(UPHASE,gpu)
  #else
      PH1(1) = 86.0_gpu
      write(io(1),*) ':mars_read : U Phase DEFAULT : ', PH1(1)
  #endif
  #if defined (MPHASE)
      write(io(1),*) ':mars_read : M Phase set to  : ', MPHASE
      PH1(2) = real(MPHASE,gpu)
  #else
      PH1(2) = 0.0_gpu
      write(io(1),*) ':mars_read : M Phase DEFAULT : ', PH1(2)
  #endif
  #if defined (LPHASE)
      write(io(1),*) ':mars_read : L Phase set to  : ', LPHASE
      PH1(3) = real(LPHASE,gpu)
  #else
      PH1(3) = 34.0_gpu
      write(io(1),*) ':mars_read : L Phase DEFAULT : ', PH1(3)
  #endif

#else
      PH1(1) = 0.0_gpu
      PH1(2) = 0.0_gpu
#endif

#if defined (COILROW)
      write(io(1),*) ':mars_read : Coil Row separation requested - blank rows...'
#if defined (MATCH)
      write(io(1),*) ':mars_read : MATCH option does not work with COILROW!'
      stop
#endif
#if defined (OLD)
      write(io(1),*) ':mars_read : OLD option does not work with COILROW!'
      stop
#endif
      if( COILROW < 1 .or. COILROW > size(IMUL) )then
          write(io(1),*) ':mars_read : COILROW does not match IMUL size!'
          stop
      endif

      do i=1,size(IMUL)
         if( i/=COILROW )IMUL(i) = 0.0_gpu
      enddo

      write(io(1),'(A20,3F11.8)') ':mars_read : IMUL =', IMUL
#endif

#if defined (MATCH)

!     Load data:

      write(io(1),*) ':mars_read : Load matching data....'

      open ( unit=lun, file=mtch, form='formatted' )

      read( lun, '(I8,2E16.8)' ) nR_, Rmin_, Rmax_
      read( lun, '(I8,2E16.8)' ) nZ_, Zmin_, Zmax_
      read( lun, '(I8,2E16.8)' ) nT_, Tmin_, Tmax_

      allocate( BMTCH(nR_,nZ_,nT_,3) )

      write(io(1),*) ':mars_read : Reading matching field....'
      write(io(1),*) ':mars_read : Tmin : ', Tmin_/360.0_gpu
      write(io(1),*) ':mars_read : Tmax : ', Tmax_/360.0_gpu
      write(io(1),*) ':mars_read : WARNING : Matching produces a bespoke&
                     & BPLASMA grid!'

!     NOTE: The resulting BPLASMA file is not corrected

      read( lun, '(6E16.8)' ) BMTCH

      close( lun )

      R0F_ = Rmin_
      R1F_ = Rmax_
      Z0F_ = Zmin_
      Z1F_ = Zmax_
      nR_f = nR_
      nZ_f = nZ_

#elif defined (UMP3D)

!     Load data:
 
      write(io(1),*) ':mars_read : Load matching data....'

      open ( unit=lun, file=mtch, form='formatted' )

      read( lun, '(I8,2E16.8)' ) nR_, Rmin_, Rmax_
      read( lun, '(I8,2E16.8)' ) nZ_, Zmin_, Zmax_
      read( lun, '(I8,2E16.8)' ) nT_, Tmin_, Tmax_

      close( lun )

      R0F_ = Rmin_
      R1F_ = Rmax_
      Z0F_ = Zmin_
      Z1F_ = Zmax_

      nR_f  =nint((R1F_-R0F_)/dXR) + 1
      nZ_f  =nint((Z1F_-Z0F_)/dXZ) + 1

      nT_f = 37 ! + 1 ! N.B. cleaned version is missing the mirrored data at 
                       ! 2pi.  However, for some unknown reason, the tricubic
                       ! spline generator from PSPLINES intriduced a series
                       ! of glitches between nodes in the toroidal direction
                       ! at the first "peak" in PHI - removing one point
                       ! in the PHI direction fixes the problem.
                       !
                       ! Moral of the story : ALWAYS CHECK THAT THE FIELD
                       ! MAP IS CLEAN!
      write(io(1),*) ':mars_read : OVERRIDE RESOLUTION : ', &
          (R1F_-R0F_)/real(nR_f-1,gpu), (Z1F_-Z0F_)/real(nZ_f-1,gpu)

#else

      R0F_  = R0F
      R1F_  = R1F
      Z0F_  = Z0F
      Z1F_  = Z1F

      nR_f  = nint((R1F_-R0F_)/dXR) + 1
      nZ_f  = nint((Z1F_-Z0F_)/dXZ) + 1

#endif

      iRF   = max(nR_f/100,1)
      iZF   = max(nZ_f/100,1)
      irbox = max(nint(0.02_gpu*(R1F_-R0F_)/dXR),1)
      izbox = max(nint(0.02_gpu*(Z1F_-Z0F_)/dXZ),1)

      write(io(1),*) ':mars_read : # R nodes           : ', nR_f
      write(io(1),*) ':mars_read : # Z nodes           : ', nZ_f
      write(io(1),*) ':mars_read : Radial   patch size : ', iRF*2+1
      write(io(1),*) ':mars_read : Vertical     -"-    : ', iZF*2+1
      write(io(1),*) ':mars_read : Rmin,Rmax           : ', R0F_,R1F_
      write(io(1),*) ':mars_read : Zmin,Zmax           : ', Z0F_,Z1F_

!     Load data:

      write(io(1),*) ':mars_read : Load MARS-F data....'

#ifndef OLD

      write(io(1),*) ':mars_read : Assimilate 3 coil sets and scale.'

#ifndef MATCH
#if( TOKAMAK==1 )
      do j=1,3
#else
      do j=1,2
#endif

#else
      do j=MATCH,MATCH
#endif

      allocate( ddat(2) )

      write(io(1),*) ':mars_read : FILE : ', TRIM(file)//TRIM(TAIL(j))

      open ( unit=lun, file=TRIM(file)//TRIM(TAIL(j)), form='formatted' )
#else

      write(io(1),*) ':mars_read : 3 coil sets pre-combined and pre-scaled.'

      do j=1,1

      allocate( ddat(2) )

      open ( unit=lun, file=file, form='formatted' )
#endif

      read ( lun, * ) idat, ddat
      
#if (TOKAMAK==2)
      print*, '********** FUDGE **********'
      ddat(1) = 1.0d3
      ddat(2) = 1.0d0
#endif

      nmde  = idat(1)

      if( nmde >= 0 )then
          write(io(1),*) ':mars_read : ERROR : MARS-F n should be -ve!'
          stop
      endif

      M0    = idat(2)
      M1    = idat(3)
      M2    = idat(4)
      NR    = idat(5)
      NV    = idat(6)
      ICOIL = idat(7)
      ICURR = ddat(1)
      IN_IC = ddat(2)

#if (TOKAMAK==1)
      if (IN_IC==1.0)then ! MARSF data does not contain coil scaling - calculate here:
          
          write(io(1),*) ':mars_read : WARNING : IN_IC=1 in MARS-F file --> no harmonic scaling due to finite coil width!'
          write(io(1),*) ':mars_read : WARNING : recalculating IN_IC to re-enable scaling!'
          IN_IC=N_coils*sin(NC*dPH(j)*pi/(2.0*180.0_gpu))/(NC*pi)

      !checks:
        !ITER IN_IC=
                    !n   Upper   Middle  Lower
                    !3   0.6645  0.4968  0.6840
                    !4   0.6126  0.4774  0.6264
                    !5   0.5494  0.4530  0.5565
                    !6   0.4772  0.4243  0.4773
                    !12  0.0150  0.1946 -0.0125
                    !15 -0.1240  0.0754 -0.1436
                    !21 -0.1065 -0.0867 -0.0872
      endif
#endif

      NT    = NR+NV
      
      write(io(1),*) ':mars_read : n     :',nmde
      write(io(1),*) ':mars_read : M0    :',M0
      write(io(1),*) ':mars_read : M1    :',M1
      write(io(1),*) ':mars_read : M2    :',M2
      write(io(1),*) ':mars_read : NR    :',NR
      write(io(1),*) ':mars_read : NV    :',NV
      write(io(1),*) ':mars_read : ICOIL :',ICOIL
      write(io(1),*) ':mars_read : ICURR :',ICURR*1.0e-3_gpu,'kAt'
      write(io(1),*) ':mars_read : -REQ- :',IKATN(j),        'kAt'
      write(io(1),*) ':mars_read : In_Ic :',IN_IC
      write(io(1),*) ':mars_read : IMUL  :',IMUL(j)

!     Tuning factor:

      IN_IC = IN_IC*IMUL(j)

#ifndef MATCH
      if( j==1 )then
#endif

      allocate( S(NT), SM(NT), SH(NT*2), q(NT) )

#ifndef MATCH
      endif
#endif

      deallocate( ddat )
      allocate  ( ddat(3) )

      do i=1,NT
         read ( lun, * ) ddat

#ifndef MATCH
         if( j==1 )then
#endif
         S (i)     = ddat(1)
         SM(i)     = ddat(2)
         q (i)     = ddat(3)
         SH(2*i-1) = ddat(1)
         SH(2*i  ) = ddat(2)

#ifndef MATCH
         endif
#endif
      enddo

      deallocate( ddat    )
      allocate  ( ddat(8) )

#ifndef MATCH
      if(j==1)then
#endif
            
      allocate( RF_Re (NT,  M0), RF_Im (NT,  M0),                             &
                ZF_Re (NT,  M0), ZF_Im (NT,  M0) )
      allocate( RFM_Re(NT,  M0), RFM_Im(NT,  M0),                             &
                ZFM_Re(NT,  M0), ZFM_Im(NT,  M0) )
      allocate( RFH_Re(NT*2,M0), RFH_Im(NT*2,M0),                             &
                ZFH_Re(NT*2,M0), ZFH_Im(NT*2,M0) )
      
      allocate( B1_Re(NT,M2-M1+1), B1_Im(NT,M2-M1+1) )
      allocate( B2_Re(NT,M2-M1+1), B2_Im(NT,M2-M1+1) )
      allocate( B3_Re(NT,M2-M1+1), B3_Im(NT,M2-M1+1) )
      
      B1_Re = 0.0_gpu
      B1_Im = 0.0_gpu
      B2_Re = 0.0_gpu
      B2_Im = 0.0_gpu
      B3_Re = 0.0_gpu
      B3_Im = 0.0_gpu

#ifndef MATCH
      endif
#endif

      do M=1,M0
         do i=1,NT
            read ( lun, * ) ddat
#ifndef MATCH
            if( j==1 )then
#endif

            RF_Re (i,    M) = ddat(1)
            RF_Im (i,    M) = ddat(2)
            ZF_Re (i,    M) = ddat(3)
            ZF_Im (i,    M) = ddat(4)
            RFM_Re(i,    M) = ddat(5)
            RFM_Im(i,    M) = ddat(6)
            ZFM_Re(i,    M) = ddat(7)
            ZFM_Im(i,    M) = ddat(8)
            RFH_Re(2*i-1,M) = ddat(1)
            RFH_Im(2*i-1,M) = ddat(2)
            ZFH_Re(2*i-1,M) = ddat(3)
            ZFH_Im(2*i-1,M) = ddat(4)
            RFH_Re(2*i  ,M) = ddat(5)
            RFH_Im(2*i  ,M) = ddat(6)
            ZFH_Re(2*i  ,M) = ddat(7)
            ZFH_Im(2*i  ,M) = ddat(8)

#ifndef MATCH
            endif
#endif
         enddo
      enddo
      
      deallocate( ddat    )
      allocate  ( ddat(6) )

!     Phase offset. See notes in header. The Min_n3 phase offsets applied to
!     the data in the ITER coordinate sustem for a +|n| mode require a
!     [+86,0,+34]*|n| phase rotation. -nmde = +ve.

#if (TOKAMAK==1)

if( NC == 3 )then
      cp = +cos( -3.0_gpu*PH1(j)*pi/180.0_gpu )
      sp =  sin( -3.0_gpu*PH1(j)*pi/180.0_gpu )
elseif( NC==6 )then
      cp = +cos( (-9.0_gpu*PH0(j) + 3.0_gpu*PH1(j) )*pi/180.0_gpu )
      sp =  sin( (-9.0_gpu*PH0(j) + 3.0_gpu*PH1(j) )*pi/180.0_gpu )
elseif( NC==12 )then
      cp = +cos( (-9.0_gpu*PH0(j) - 3.0_gpu*PH1(j) )*pi/180.0_gpu )
      sp =  sin( (-9.0_gpu*PH0(j) - 3.0_gpu*PH1(j) )*pi/180.0_gpu )
elseif( NC==15 )then
      cp = +cos( (-18.0_gpu*PH0(j) + 3.0_gpu*PH1(j) )*pi/180.0_gpu )
      sp =  sin( (-18.0_gpu*PH0(j) + 3.0_gpu*PH1(j) )*pi/180.0_gpu )
endif

#else
if( NC == 2 )then
      cp = +cos( -2.0_gpu*PH1(j)*pi/180.0_gpu )
      sp =  sin( -2.0_gpu*PH1(j)*pi/180.0_gpu )
else
      stop
endif
#endif

      write(io(1),*) ':mars_read : PHASE SHIFT:',                             &
                       -nmde*( PH0(j)-PH1(j) ),'deg'

      do M=1,M2-M1+1
         do i=1,NT
            read( lun, * ) ddat

            B1_Re (i,M) = B1_Re (i,M) + (ddat(1)*cp - ddat(2)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
            B1_Im (i,M) = B1_Im (i,M) + (ddat(2)*cp + ddat(1)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
            B2_Re (i,M) = B2_Re (i,M) + (ddat(3)*cp - ddat(4)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
            B2_Im (i,M) = B2_Im (i,M) + (ddat(4)*cp + ddat(3)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
            B3_Re (i,M) = B3_Re (i,M) + (ddat(5)*cp - ddat(6)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
            B3_Im (i,M) = B3_Im (i,M) + (ddat(6)*cp + ddat(5)*sp)*IKATN(j)*   &
                          IN_IC/(ICURR*1.0e-03_gpu)
         enddo
      enddo
      
      close( unit=lun )

      deallocate( ddat )

      enddo

!     The MARS-F contravariant field vectors B1, B2 and B3 now correspond to
!     the +|n| field in ITER coordinates.

#if defined (MATCH)

      write(io(1),*) ':mars_read : Fourier decompose matching field:'

      allocate( PHI_(nT_), Re_F(nR_,nZ_,3), Im_F(nR_,nZ_,3) )

      PHI_ = (/(real(i,gpu),i=0,nT_-1)/)*2.0_gpu*pi/real(nT_,gpu)

      write(io(1),*) ':mars_read : Range : ', PHI_(1), PHI_(nT_)/(2.0_gpu*pi)

!     Fourier decompose for +n (note the -ve sign infront of nmde):

      do i=1,nR_
         do j=1,nZ_
            do k=1,3
               Re_F(i,j,k) = sum(BMTCH(i,j,1:nT_,k)*                          &
                             cos(PHI_*(-nmde)))*(PHI_(2)-PHI_(1))/pi
               Im_F(i,j,k) = sum(BMTCH(i,j,1:nT_,k)*                          &
                             sin(PHI_*(-nmde)))*(PHI_(2)-PHI_(1))/pi
            enddo
         enddo
      enddo

#endif

!     Determine inner and outer zones:

      do i=1,2*NT
         if( SH(i)>1.0_gpu )then
             NT_in  = i-1
             NT_out = 2*NT-NT_in+1 ! Include separatrix
             exit
         endif
      enddo

      write(io(1),*) ':mars_read : [R,Z]   Inner zone range ',NT_in
      write(io(1),*) ':mars_read : [R,Z]   Outer zone range ',NT_out

!     Build spline coefficients for spatial mapping:

      allocate( RFH_Re_fspl_in ( 2, NT_in,  M0 ),                             &
                RFH_Im_fspl_in ( 2, NT_in,  M0 ) )
      allocate( ZFH_Re_fspl_in ( 2, NT_in,  M0 ),                             &
                ZFH_Im_fspl_in ( 2, NT_in,  M0 ) )
      allocate( RFH_Re_fspl_out( 2, NT_out, M0 ),                             &
                RFH_Im_fspl_out( 2, NT_out, M0 ) )
      allocate( ZFH_Re_fspl_out( 2, NT_out, M0 ),                             &
                ZFH_Im_fspl_out( 2, NT_out, M0 ) )

!     Inner spline (s<=1):

#if defined (CRESET)
      do i=1,2
#else
      do i=1,1
#endif

      allocate( fspl( 2, NT_in ) )

      do M=1,M0

         fspl( 1, 1:NT_in ) = RFH_Re(1:NT_in,M)

         call mkspline( SH(1:NT_in), NT_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         RFH_Re_fspl_in( 1:2, 1:NT_in, M ) = fspl( 1:2, 1:NT_in )

         fspl( 1, 1:NT_in ) = RFH_Im(1:NT_in,M)

         call mkspline( SH(1:NT_in), NT_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )



         RFH_Im_fspl_in( 1:2, 1:NT_in, M ) = fspl( 1:2, 1:NT_in )

         fspl( 1, 1:NT_in ) = ZFH_Re(1:NT_in,M)
      
         call mkspline( SH(1:NT_in), NT_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         ZFH_Re_fspl_in( 1:2, 1:NT_in, M ) = fspl( 1:2, 1:NT_in )

         fspl( 1, 1:NT_in ) = ZFH_Im(1:NT_in,M)
      
         call mkspline( SH(1:NT_in), NT_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         ZFH_Im_fspl_in( 1:2, 1:NT_in, M ) = fspl( 1:2, 1:NT_in )

      enddo

      deallocate( fspl              )

!     Outer spline (s>1):

      allocate  ( fspl( 2, NT_out ) )

      do M=1,M0

         fspl( 1, 1:NT_out ) = RFH_Re(NT_in:NT*2,M)

         call mkspline( SH(NT_in:NT*2), NT_out, fspl, 0, 0.0_gpu, 0,          &
                         0.0_gpu, ilinx, iflag )

         RFH_Re_fspl_out( 1:2, 1:NT_out, M ) = fspl( 1:2, 1:NT_out )

         fspl( 1, 1:NT_out ) = RFH_Im(NT_in:NT*2,M)
      
         call mkspline( SH(NT_in:NT*2), NT_out, fspl, 0, 0.0_gpu, 0,          &
                        0.0_gpu, ilinx, iflag )

         RFH_Im_fspl_out( 1:2, 1:NT_out, M ) = fspl( 1:2, 1:NT_out )

         fspl( 1, 1:NT_out ) = ZFH_Re(NT_in:NT*2,M)
      
         call mkspline( SH(NT_in:NT*2), NT_out, fspl, 0, 0.0_gpu, 0,          &
                        0.0_gpu, ilinx, iflag )

         ZFH_Re_fspl_out( 1:2, 1:NT_out, M ) = fspl( 1:2, 1:NT_out )

         fspl( 1, 1:NT_out ) = ZFH_Im(NT_in:NT*2,M)
      
         call mkspline( SH(NT_in:NT*2), NT_out, fspl, 0, 0.0_gpu, 0,          &
                        0.0_gpu, ilinx, iflag )

         ZFH_Im_fspl_out( 1:2, 1:NT_out, M ) = fspl( 1:2, 1:NT_out )

      enddo

#if defined (CRESET)

      if( i==1 )then

      write(io(1),'(A41,2F18.15)') ':mars_read : Rmag(file), Zmag(file)   : ',&
                                     RFH_Re(1,1), ZFH_Re(1,1)

!     Reset Rmag,Zmag:

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.8_gpu,            &
                  0.5_gpu*pi, R, Z, 0, iflag )

      dat(1) = R(1)

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.8_gpu,            &
                  1.5_gpu*pi, R, Z, 0, iflag )

      dat(2) = R(1)

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.8_gpu,            &
                  0.0_gpu, R, Z, 0, iflag )

      dat(3) = Z(1)

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.8_gpu,            &
                  pi, R, Z, 0, iflag )

      dat(4) = Z(1)

      Rmag        =  (dat(1)+dat(2))*0.5_gpu
      Zmag        =  (dat(3)+dat(4))*0.5_gpu

      RFH_Re(1,1) = Rmag
      ZFH_Re(1,1) = Zmag

      S (1)       = 0.0_gpu
      SH(1)       = 0.0_gpu

      write(io(1),'(A41,2F18.15)') ':mars_read : Rmag(tuned), Zmag(tuned) : ',&
                    RFH_Re(1,1), ZFH_Re(1,1) 

      endif

#endif

      deallocate( fspl )

      enddo

!     Determine splines B1, B2 and B3:

!     Determine inner and outer zones for B1, B2 & B3:

      do i=1,NT
         if( S(i)>1.0_gpu )then
             N1_in  = i-1
             N1_out = NT-N1_in+1 ! Include separatrix
             exit
         endif
      enddo
      do i=1,NT
         if( SM(i)>1.0_gpu )then
             N2_in  = i-1
             N2_out = NT-N2_in   ! Exclude separatrix (half grid)
             exit
         endif
      enddo
      
      N3_in  = N2_in
      N3_out = N2_out
      
      write(io(1),*) ':mars_read : [B1 ]   Inner zone range ',N1_in
      write(io(1),*) ':mars_read : [B1 ]   Outer     -"-    ',N1_out
      write(io(1),*) ':mars_read : [B2,B3] Inner zone range ',N2_in
      write(io(1),*) ':mars_read : [B2,B3] Outer     -"-    ',N2_out
      
      allocate( B1_Re_fspl_in ( 2, N1_in,  M2-M1+1 ),                         &
                B1_Im_fspl_in ( 2, N1_in,  M2-M1+1 ) )
      allocate( B2_Re_fspl_in ( 2, N2_in,  M2-M1+1 ),                         &
                B2_Im_fspl_in ( 2, N2_in,  M2-M1+1 ) )
      allocate( B3_Re_fspl_in ( 2, N3_in,  M2-M1+1 ),                         &
                B3_Im_fspl_in ( 2, N3_in,  M2-M1+1 ) )
      allocate( B1_Re_fspl_out( 2, N1_out, M2-M1+1 ),                         &
                B1_Im_fspl_out( 2, N1_out, M2-M1+1 ) )
      allocate( B2_Re_fspl_out( 2, N2_out, M2-M1+1 ),                         &
                B2_Im_fspl_out( 2, N2_out, M2-M1+1 ) )
      allocate( B3_Re_fspl_out( 2, N3_out, M2-M1+1 ),                         &
                B3_Im_fspl_out( 2, N3_out, M2-M1+1 ) )

!     Inner spline (s<=1):

      do M=1,M2-M1+1

         allocate( fspl( 2, N1_in ) )
      
         fspl( 1, 1:N1_in ) = B1_Re(1:N1_in,M)
      
         call mkspline( S(1:N1_in),  N1_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B1_Re_fspl_in( 1:2, 1:N1_in, M ) = fspl( 1:2, 1:N1_in )

         fspl( 1, 1:N1_in ) = B1_Im(1:N1_in,M)
      
         call mkspline( S(1:N1_in),  N1_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B1_Im_fspl_in( 1:2, 1:N1_in, M ) = fspl( 1:2, 1:N1_in )

         deallocate( fspl              )

         allocate( fspl( 2, N2_in ) )
      
         fspl( 1, 1:N2_in ) = B2_Re(1:N2_in,M)
      
         call mkspline( SM(1:N2_in), N2_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B2_Re_fspl_in( 1:2, 1:N2_in, M ) = fspl( 1:2, 1:N2_in )

         fspl( 1, 1:N2_in ) = B2_Im(1:N2_in,M)
      
         call mkspline( SM(1:N2_in), N2_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B2_Im_fspl_in( 1:2, 1:N2_in, M ) = fspl( 1:2, 1:N2_in )

         fspl( 1, 1:N3_in ) = B3_Re(1:N3_in,M)
      
         call mkspline( SM(1:N3_in), N3_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B3_Re_fspl_in( 1:2, 1:N3_in, M ) = fspl( 1:2, 1:N3_in )

         fspl( 1, 1:N3_in ) = B3_Im(1:N3_in,M)
      
         call mkspline( SM(1:N3_in), N3_in, fspl, 0, 0.0_gpu, 0, 0.0_gpu,     &
                        ilinx, iflag )

         B3_Im_fspl_in( 1:2, 1:N3_in, M ) = fspl( 1:2, 1:N3_in )

         deallocate( fspl              )

      enddo

!     Outer spline (s>1):

      do M=1,M2-M1+1

         allocate( fspl( 2, N1_out ) )
      
         fspl( 1, 1:N1_out ) = B1_Re(N1_in:NT,M)
      
         call mkspline( S(N1_in:NT),  N1_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu,   &
                        ilinx, iflag )

         B1_Re_fspl_out( 1:2, 1:N1_out, M ) = fspl( 1:2, 1:N1_out )

         fspl( 1, 1:N1_out ) = B1_Im(N1_in:NT,M)
      
         call mkspline( S(N1_in:NT),  N1_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu,   &
                        ilinx, iflag )

         B1_Im_fspl_out( 1:2, 1:N1_out, M ) = fspl( 1:2, 1:N1_out )

         deallocate( fspl              )

         allocate( fspl( 2, N2_out ) )

         fspl( 1, 1:N2_out ) = B2_Re(N2_in+1:NT,M)
      
         call mkspline( SM(N2_in+1:NT), N2_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu, &
                        ilinx, iflag )

         B2_Re_fspl_out( 1:2, 1:N2_out, M ) = fspl( 1:2, 1:N2_out )

         fspl( 1, 1:N2_out ) = B2_Im(N2_in+1:NT,M)
      
         call mkspline( SM(N2_in+1:NT), N2_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu, &
                        ilinx, iflag )

         B2_Im_fspl_out( 1:2, 1:N2_out, M ) = fspl( 1:2, 1:N2_out )

         fspl( 1, 1:N3_out ) = B3_Re(N3_in+1:NT,M)
      
         call mkspline( SM(N3_in+1:NT), N3_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu, &
                        ilinx, iflag )

         B3_Re_fspl_out( 1:2, 1:N3_out, M ) = fspl( 1:2, 1:N3_out )

         fspl( 1, 1:N3_out ) = B3_Im(N3_in+1:NT,M)
      
         call mkspline( SM(N3_in+1:NT), N3_out, fspl, 0, 0.0_gpu, 0, 0.0_gpu, &
                        ilinx, iflag )

         B3_Im_fspl_out( 1:2, 1:N3_out, M ) = fspl( 1:2, 1:N3_out )

         deallocate( fspl              )
         
      enddo

!     Find magnetic axis:

      s_     = 0.0_gpu ! SH(1)    ! Note, not quite == 0.0
      chi_   = 0.0_gpu

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_, chi_,           &
                  R, Z, 0, iflag )

      R_mag_ = R(1)
      Z_mag_ = Z(1)

!     Find LCFS:

      write(io(1),*) ':mars_read : Locate LCFS....'

      allocate( rho_2(NL+1), tht(NL+1) )

      do i=0,NL

         s_   =  1.0_gpu
         chi_ = -pi + 2.0_gpu*pi*real(i,gpu)/real(NL,gpu)

         call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,          &
                     ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_, chi_,        &
                     R, Z, 0, iflag )

         rho_2(i+1) = ( R(1)-R_mag_ )**2 + ( Z(1)-Z_mag_ )**2

         if( Z(1)==Z_mag_ .and. R(1)==R_mag_ )then
             tht(i+1) = 0.0_gpu
         else
             tht(i+1) = atan2( Z(1)-Z_mag_, R(1)-R_mag_ )
         endif

!        Ensure tht is strictly ascending:

         if( i>=1 )then
             if( tht(i+1) < tht(i) )tht(i+1) = tht(i+1) + 2.0_gpu*pi
         endif
      enddo

      allocate  ( fspl( 2, NL+1 ) )

      fspl( 1, 1:NL+1 ) = rho_2( 1:NL+1 )
      
      call mkspline( tht, NL+1, fspl, 0, 0.0_gpu, 0, 0.0_gpu, ilinx, iflag )
      
!     Build grid:
      
      allocate( cells_h( nR_f, nZ_f, 7 ) )
      
      cells_h( 1:nR_f,1:nZ_f,7 ) = HUGE(0.0_gpu)
      
      iter = 0

!     Use M.C. scatter to get "close" to solution:

      write(io(1),*) ':mars_read : M.C. scatter stage....'

!----- Main GRID finding loop -----

      write( STR_R, '(I4)' ) nR_f
      write( STR_Z, '(I4)' ) nZ_f

!     n mode number string:

      write( STR_N, '(I2)' ) NC

#if defined (COILROW)

!     Coil row tag:

#if (TOKAMAK==1)

      if( COILROW<1 .or. COILROW>3 )then     
          write(io(1),*) ':mars_read : ERROR : ITER has 3 coil rows!'
          stop
      else

          STR_C = TAG_I(COILROW)

      endif
#else

      write( STR_C, '(I1)' ) COILROW
#endif

#endif

      open( unit=lun, file=TRIM(ADJUSTL(root))//'mars_read_'//                &
                                         TRIM(ADJUSTL(STR_R))//'_'//          &
                                         TRIM(ADJUSTL(STR_Z))//'_'//          &
                                         TRIM(ADJUSTL(STR_N))//'.CACHE',      &
                      status='old', form='unformatted', iostat=ios )

      if( ios==0 )then

          write(io(1),*) ':mars_read : CACHE file identified : '//            &
                          TRIM(ADJUSTL(root))//                               &
                          'mars_read_'//                                      &
                                        TRIM(ADJUSTL(STR_R))//'_'//           &
                                        TRIM(ADJUSTL(STR_Z))//'.CACHE'

          read( lun ) nR_f_, nZ_f_
          read( lun ) R0C_, R1C_, Z0C_, Z1C_

          if( nR_f_ /= nR_f .or. nZ_f_ /= nZ_f .or.                           &
              R0C_  /= R0F_ .or. R1C_  /= R1F_ .or.                           &
              Z0C_  /= Z0F_ .or. Z1C_  /= Z1F_ )then

              write(io(1),*) ':mars_read : CACHE does not match : ', nR_f_,   &
                                                                     nZ_f_
              ios = 1
              close( unit=lun )
          endif
      endif

      if( ios==0 )then

          write(io(1),*) ':mars_read : Load spatial pre-mapping from CACHE'

          allocate( cells( nR_f, nZ_f, 7 ) )
          read( lun ) cells_h

          close( unit=lun )

      else

!$OMP PARALLEL PRIVATE( i,j, s_, chi_, R, Z, iflag )                          &
!$OMP          SHARED ( DUM )

!$OMP DO
      do i=0,1000
         do j=0,360

         s_   = 1.5_gpu*real(i,gpu)/1000.0_gpu
         chi_ = 2.0_gpu*pi*real(j,gpu)/360.0_gpu


         if( s_ <= 1.0_gpu )then
             call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,      &
                         ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_, chi_, R, &
                         Z, 0, iflag )
         else
             call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,ZFH_Re_fspl_out,     &
                         ZFH_Im_fspl_out, SH(NT_in:NT_out), NT_out, M0, s_,   &
                         chi_, R, Z, 0, iflag )
         endif

         DUM(i+1,j+1,1) = R(1)
         DUM(i+1,j+1,2) = Z(1)

         enddo
      enddo
!$OMP END DO

!$OMP END PARALLEL

      allocate(RAN (2,NS))
      allocate(RAN_(NS))

      do

      iter = iter + 1

      write(io(1),*) ':mars_read : iteration : ', iter

      call random_number( RAN  )
      call random_number( RAN_ )

!$OMP PARALLEL PRIVATE( i,cells, s_, chi_, R, Z, iR, iZ, R_p, Z_p, d_p, j, k, &
!$OMP                   det, ds, dchi,tht_,f, iflag )                         &
!$OMP          SHARED ( RAN,RAN_,ict,tht,fspl )

      allocate( cells( nR_f, nZ_f, 7 ) )
      
      cells( 1:nR_f,1:nZ_f,7 ) = HUGE(0.0_gpu)
      
!$OMP DO

      do i=1,NS

         if( mod(i,NS/100)==0 )write(io(1),*) i

         if( RAN_(i) < 0.5_gpu )then

#if (TOKAMAK==1)
#ifndef BPLASMAOLD
             s_   = min(max((RAN(1,i)**0.5_gpu)*2.0_gpu,SH(1)),SH(2*NT))
#else
             s_   = min(max((RAN(1,i)**0.5_gpu)*2.5_gpu,SH(1)),SH(2*NT))
#endif
#else
             s_   = min(max((RAN(1,i)**0.5_gpu)*3.5_gpu,SH(1)),SH(2*NT))
#endif
             chi_ = 2.0_gpu*dpi*RAN(2,i)

         else

#if (TOKAMAK==1)
             s_   = 0.99_gpu + RAN(1,i)*0.02_gpu
             chi_ = 4.35_gpu + RAN(2,i)*0.15_gpu
#else
             s_   = 0.80_gpu + RAN(1,i)*0.05_gpu
             chi_ = 4.15_gpu + RAN(2,i)*0.55_gpu
#endif
         endif

         if( s_ <= 1.0_gpu )then
             call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,      &
                         ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_, chi_, R, &
                         Z, 0, iflag )
         else
             call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,ZFH_Re_fspl_out,     &
                         ZFH_Im_fspl_out, SH(NT_in:NT_out), NT_out, M0, s_,   &
                         chi_, R, Z, 0, iflag )
         endif

         iR  = nint( real(nR_f-1,gpu)*(R(1)-R0F_)/(R1F_-R0F_) ) + 1
         iZ  = nint( real(nZ_f-1,gpu)*(Z(1)-Z0F_)/(Z1F_-Z0F_) ) + 1

         if( iR >= 1 .and. iR <= nR_f .and.                                   &
             iZ >= 1 .and. iZ <= nZ_f )then
             
             R_p = R0F_ + (R1F_-R0F_)*real(iR-1,gpu)/real(nR_f-1,gpu)
             Z_p = Z0F_ + (Z1F_-Z0F_)*real(iZ-1,gpu)/real(nZ_f-1,gpu)
             d_p = (R(1)-R_p)**2 + (Z(1)-Z_p)**2

             if( d_p < cells(iR,iZ,7) )then
                 if( cells(iR,iZ,7)==HUGE(0.0_gpu) )then

!                    Only fill if s<=1 and PSIn<=1 or s>0 and PSIn>0....

                     if( Z_p==Z_mag_ .and. R_p==R_mag_ )then
                         tht_ = 0.0_gpu
                     else
                         tht_ = atan2( Z_p-Z_mag_, R_p-R_mag_ )
                     endif

                     if(tht_ < tht(1))tht_=tht_+2.0_gpu*pi

                     ict  = [1,0,0]

                     call evspline_( tht_, tht, NL+1, 2, fspl( 1:2, 1:NL+1 ), &
                                     ict, f, ierr )
  
                     f(1) = ( (R_p-R_mag_)**2+(Z_p-Z_mag_)**2 )/f(1)

                     det  = R(2)*Z(3) - R(3)*Z(2)

                     if( ( (s_ <= 1.0_gpu .and. f(1) <= 1.0_gpu) .or.         &
                           (s_ >  1.0_gpu .and. f(1) >  1.0_gpu) ) .and.      &
                            abs(det) /= 0.0_gpu )then

                          cells(iR,iZ,1) = s_
                          cells(iR,iZ,2) = chi_
                          cells(iR,iZ,3) = R_p-R(1)
                          cells(iR,iZ,4) = Z_p-Z(1)



                          ds   = ( Z(3)*cells(iR,iZ,3)-R(3)*cells(iR,iZ,4))/det
                          dchi = (-Z(2)*cells(iR,iZ,3)+R(2)*cells(iR,iZ,4))/det

                          cells(iR,iZ,5) = ds
                          cells(iR,iZ,6) = dchi
                          cells(iR,iZ,7) = d_p
                     endif

!                    Also fill in viscinity....

                     do j=max(iR-iRF,1),min(iR+10,nR_f)
                        do k=max(iZ-iZF,1),min(iZ+10,nZ_f)
                        
                           R_p = R0F_ + (R1F_-R0F_)*real(j-1,gpu)/            &
                                 real(nR_f-1,gpu)
                           Z_p = Z0F_ + (Z1F_-Z0F_)*real(k-1,gpu)/            &
                                 real(nZ_f-1,gpu)

!                          Only fill if s<=1 and PSIn<=1 or s>0 and PSIn>0....

                           if( Z_p==Z_mag_ .and. R_p==R_mag_ )then
                               tht_ = 0.0_gpu
                           else
                               tht_ = atan2( Z_p-Z_mag_, R_p-R_mag_ )
                           endif

                           if(tht_ < tht(1))tht_=tht_+2.0_gpu*pi

                           ict  = [1,0,0]

                           call evspline_( tht_, tht, NL+1, 2,                &
                                           fspl( 1:2, 1:NL+1 ), ict, f, ierr )
  
                           f(1) = ( (R_p-R_mag_)**2+(Z_p-Z_mag_)**2 )/f(1)

                           if( (s_ <= 1.0_gpu .and. f(1) <= 1.0_gpu) .or.     &
                               (s_ >  1.0_gpu .and. f(1) >  1.0_gpu) )then

                               det  = R(2)*Z(3) - R(3)*Z(2)

                               if( cells(j,k,7)==HUGE(0.0_gpu) .and.          &
                                   abs(det) /= 0.0_gpu )then
                               
                                   d_p = (R(1)-R_p)**2 + (Z(1)-Z_p)**2  

                                   cells(j,k,1) = s_
                                   cells(j,k,2) = chi_
                                   cells(j,k,3) = R_p-R(1)
                                   cells(j,k,4) = Z_p-Z(1)

                                   ds   = ( Z(3)*cells(j,k,3)-                &
                                            R(3)*cells(j,k,4))/det
                                   dchi = (-Z(2)*cells(j,k,3)+                &
                                            R(2)*cells(j,k,4))/det

                                   cells(j,k,5) = ds
                                   cells(j,k,6) = dchi
                                   cells(j,k,7) = d_p
                               endif

                           endif
                        enddo
                     enddo

                 else

!                    Only fill if s<=1 and PSIn<=1 or s>0 and PSIn>0....

                     if( Z_p==Z_mag_ .and. R_p==R_mag_ )then
                         tht_ = 0.0_gpu
                     else
                         tht_ = atan2( Z_p-Z_mag_, R_p-R_mag_ )
                     endif

                     if( tht_ < tht(1) )tht_=tht_+2.0_gpu*pi

                     ict  = [1,0,0]

                     call evspline_( tht_, tht, NL+1, 2, fspl( 1:2, 1:NL+1 ), &
                                     ict, f, ierr )
  
                     f(1) = ((R_p-R_mag_)**2+(Z_p-Z_mag_)**2)/f(1)

                     det  = R(2)*Z(3) - R(3)*Z(2)

                     if( ( (s_ <= 1.0_gpu .and. f(1) <= 1.0_gpu) .or.         &
                           (s_ >  1.0_gpu .and. f(1) >  1.0_gpu) ) .and.      &
                            abs(det) /= 0.0_gpu )then

                          cells(iR,iZ,1) = s_
                          cells(iR,iZ,2) = chi_
                          cells(iR,iZ,3) = R_p-R(1)
                          cells(iR,iZ,4) = Z_p-Z(1)

                          ds   = ( Z(3)*cells(iR,iZ,3)-R(3)*cells(iR,iZ,4))/det
                          dchi = (-Z(2)*cells(iR,iZ,3)+R(2)*cells(iR,iZ,4))/det

                          cells(iR,iZ,5) = ds
                          cells(iR,iZ,6) = dchi
                          cells(iR,iZ,7) = d_p

                     endif
                 endif
             endif
         endif
      
      enddo

!$OMP END DO

!     Reduction:

!$OMP CRITICAL

      do iR=1,nR_f
         do iZ=1,nZ_f
            if( cells(iR,iZ,7) < cells_h(iR,iZ,7) )then
                cells_h(iR,iZ,1:7) = cells(iR,iZ,1:7)
            endif
         enddo
      enddo

!$OMP END CRITICAL

      deallocate( cells )

!$OMP END PARALLEL

      if( maxval( cells_h(:,:,7) ) < HUGE(0.0_gpu) )exit

      write(io(1),*) ':mars_read: Data contains gaps - iterate....',          &
                       maxval( cells_h(:,:,7) )

      enddo

      endif

      write(io(1),*) ':mars_read : Refining....'

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.0_gpu,            &
                  0.0_gpu, R, Z, 0, iflag )
      Rmag = R(1)
      Zmag = Z(1)

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.7_gpu,            &
                  0.5_gpu*pi, R, Z, 0, iflag )

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.7_gpu,            &
                  1.5_gpu*pi, R, Z, 0, iflag )

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.7_gpu,            &
                  0.0_gpu*pi, R, Z, 0, iflag )

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,             &
                  ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.7_gpu,            &
                  1.0_gpu*pi, R, Z, 0, iflag )

!     Storage for Re and Im parts of B:

#if defined (FXXYY)
      allocate( B( nR_f, nZ_f, 7, 7, 8 )  )
      B       ( 1:nR_f,1:nZ_f,1:7,1:7,1:8 ) = 0.0_gpu
#else
      allocate( B( nR_f, nZ_f, 1, 1, 8 )  )
      B       ( 1:nR_f,1:nZ_f,1,  1,  1:8 ) = 0.0_gpu
#endif

!     ...and for derivatives:

      allocate( FDER( nR_f, nZ_f, 18 ) )

!$OMP PARALLEL PRIVATE( i, j, s_, chi_, R, Z, iR, iZ, R_p, Z_p, d_p, det, ds, &
!$OMP                   dchi, B1, B2, B3, iter, JAC, jR, jZ,                  &
!$OMP                   fxx, fyy, fxxyy, RAX, ZAX, B_, iflag )                &
!$OMP          SHARED ( cells_h, B, FDER )

      allocate( RAX(7), ZAX(7) )
      allocate( B_(7,7), fxx(7,7), fyy(7,7), fxxyy(7,7) )

      do jR=jR0,jR1
      do jZ=jZ0,jZ1

      write(io(1),*) ':mars_read : jR, jZ :', jR, jZ, omp_get_wtime()

      dR_ = (jR-4)*dROFF
      dZ_ = (jZ-4)*dZOFF

!$OMP DO

      do iR=1,nR_f
         do iZ=1,nZ_f
            do iter=1,NMAX

               s_   = max(cells_h(iR,iZ,1) + cells_h(iR,iZ,5), 0.0_gpu)
999            chi_ = cells_h(iR,iZ,2) + cells_h(iR,iZ,6)

               R_p  = R0F_ + dR_ + (R1F_-R0F_)*real(iR-1,gpu)/real(nR_f-1,gpu)
               Z_p  = Z0F_ + dZ_ + (Z1F_-Z0F_)*real(iZ-1,gpu)/real(nZ_f-1,gpu)

!chi_ = atan2( Z_p-Zmag, R_p-Rmag )

               sin_ = sin(chi_)
               cos_ = cos(chi_)
               chi_ = atan2(sin_,cos_)

               if( chi_ < 0.0_gpu    )chi_ = chi_ + 2.0_gpu*pi

               B1(1:2) = 0.0_gpu
               B2(1:2) = 0.0_gpu
               B3(1:2) = 0.0_gpu
               
               if( s_ <= 1.0_gpu )then
                   call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,&
                               ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_,    &
                               chi_, R, Z, 0, iflag )
               else
                   call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,               &
                               ZFH_Re_fspl_out,ZFH_Im_fspl_out,               &
                               SH(NT_in:NT_out), NT_out, M0, s_,              &
                               chi_, R, Z, 0, iflag )
               endif

               det  = R(2)*Z(3) - R(3)*Z(2)

               if( abs(det)==0.0_gpu )then
                   s_ = s_ + 1.0E-08_gpu
                   goto 999
               endif

               d_p = (R(1)-R_p)**2 + (Z(1)-Z_p)**2

               cells_h(iR,iZ,1) = s_
               cells_h(iR,iZ,2) = chi_

               cells_h(iR,iZ,3) = R_p-R(1)
               cells_h(iR,iZ,4) = Z_p-Z(1)

               ds   = ( Z(3)*cells_h(iR,iZ,3)-R(3)*cells_h(iR,iZ,4))/det
               dchi = (-Z(2)*cells_h(iR,iZ,3)+R(2)*cells_h(iR,iZ,4))/det

               cells_h(iR,iZ,5) = ds
               cells_h(iR,iZ,6) = dchi
               cells_h(iR,iZ,7) = d_p

               if( sqrt(d_p) < DPTOL .or. iter==NMAX )then

                   if( iter==NMAX )write(io(1),*)                             &
                       ':mars_read : WARNING : sqrt(d_p) :', sqrt(d_p), s_,   &
                         chi_, iR, iZ

!                  Determine BR, BT and BZ....

                   if( s_ <= 1.0_gpu )then

                       call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in,            &
                                   ZFH_Re_fspl_in, ZFH_Im_fspl_in,            &
                                   SH(1:NT_in), NT_in, M0, s_,                &
                                   chi_, R, Z, 1, iflag )

                       call B_F( B1_Re_fspl_in, B1_Im_fspl_in, S (1:N1_in),   &
                                 N1_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B1,  &
                                 1, iflag )                     
                       call B_F( B2_Re_fspl_in, B2_Im_fspl_in, SM(1:N2_in),   &
                                 N2_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B2,  &
                                 1, iflag )
                       call B_F( B3_Re_fspl_in, B3_Im_fspl_in, SM(1:N3_in),   &
                                 N3_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B3,  &
                                 1, iflag )
                   else

                       call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,           &
                                   ZFH_Re_fspl_out,ZFH_Im_fspl_out,           &
                                   SH(NT_in:NT_out), NT_out, M0, s_,          &
                                   chi_, R, Z, 1, iflag )

                       call B_F( B1_Re_fspl_out, B1_Im_fspl_out,              &
                                 S (N1_in:NT), N1_out, M1, M2, nmde, s_, chi_,&
                                 0.0_gpu, B1, 1, iflag )
                       call B_F( B2_Re_fspl_out, B2_Im_fspl_out,              &
                                 SM(N2_in+1:NT), N2_out, M1, M2, nmde, s_,    &
                                 chi_,0.0_gpu, B2, 1, iflag )
                       call B_F( B3_Re_fspl_out, B3_Im_fspl_out,              &
                                 SM(N3_in+1:NT), N3_out, M1, M2, nmde, s_,    &
                                 chi_,0.0_gpu, B3, 1, iflag )
                   endif
                   
                   JAC         = (R(2)*Z(3) - R(3)*Z(2))*R(1)

#if defined (FXXYY)
                   B(iR,iZ,jR,jZ, 1) = (B1(1)*R(2)+B2(1)*R(3))/JAC
                   B(iR,iZ,jR,jZ, 2) = (B1(2)*R(2)+B2(2)*R(3))/JAC
                   B(iR,iZ,jR,jZ, 3) = (B1(1)*Z(2)+B2(1)*Z(3))/JAC
                   B(iR,iZ,jR,jZ, 4) = (B1(2)*Z(2)+B2(2)*Z(3))/JAC

!                  See header: the -ve sign here flips BT into the ITER
!                  coordinate system.

                   B(iR,iZ,jR,jZ, 5) = -B3(1)*R(1)/JAC
                   B(iR,iZ,jR,jZ, 6) = -B3(2)*R(1)/JAC
#else
                   B(iR,iZ,1, 1,  1) = (B1(1)*R(2)+B2(1)*R(3))/JAC
                   B(iR,iZ,1, 1,  2) = (B1(2)*R(2)+B2(2)*R(3))/JAC
                   B(iR,iZ,1, 1,  3) = (B1(1)*Z(2)+B2(1)*Z(3))/JAC
                   B(iR,iZ,1, 1,  4) = (B1(2)*Z(2)+B2(2)*Z(3))/JAC

!                  See header: the -ve sign here flips BT into the ITER
!                  coordinate system.

                   B(iR,iZ,1, 1,  5) = -B3(1)*R(1)/JAC
                   B(iR,iZ,1, 1,  6) = -B3(2)*R(1)/JAC
#endif
                   B(iR,iZ,1, 1,  7) =  s_
                   B(iR,iZ,1, 1,  8) =  chi_

                   exit
               endif

            enddo
         enddo
      enddo

!$OMP END DO

      enddo
      enddo

#if defined (FXXYY)

!     Determine fxx, fyy and fxxyy:

!$OMP DO

      do iR=1,nR_f
         do iZ=1,nZ_f
            do i=1,6
               B_(1:7,1:7) = B(iR,iZ,1:7,1:7,i)

               RAX(1:7) = R0F_ + (R1F_-R0F_)*real(iR-1,gpu)/real(nR_f-1,gpu)+ &
                          (/((j-4)*dROFF,j=1,7)/)
               ZAX(1:7) = Z0F_ + (Z1F_-Z0F_)*real(iZ-1,gpu)/real(nZ_f-1,gpu)+ &
                          (/((j-4)*dZOFF,j=1,7)/)

               call coeff2( RAX, ZAX, B_, fxx, fyy, fxxyy, iflag )

               FDER(iR,iZ,1+(i-1)*3) = fxx  (4,4)
               FDER(iR,iZ,2+(i-1)*3) = fyy  (4,4)
               FDER(iR,iZ,3+(i-1)*3) = fxxyy(4,4)
            enddo
         enddo
      enddo

!$OMP END DO

#endif

!$OMP END PARALLEL

!     Write out mapping CACHE file:

      if( ios /= 0 )then

      write(io(1),*) ':mars_read : Write out mapping CACHE file....'

      open ( unit=lun,  file=TRIM(ADJUSTL(root))//'mars_read_'//              &
                                           TRIM(ADJUSTL(STR_R))//'_'//        &
                                           TRIM(ADJUSTL(STR_Z))//'_'//        &
                                           TRIM(ADJUSTL(STR_N))//'.CACHE',    &
             form='unformatted', status='replace' )

      write( lun ) nR_f, nZ_f
      write( lun ) R0F_, R1F_, Z0F_, Z1F_
      write( lun ) cells_h

      close( unit=lun )

      endif

      write(io(1),*) ':mars_read : Refinement complete.'

#if defined (MATCH)

      write(io(1),*) ':mars_read : Check phase and amplitude @ S: ', SMATCH

      allocate( ddat(360) )

      ddat( 1:9 ) = 0.0_gpu

      do k=-180,179

      chi_ = real(k,gpu)*pi/180.0_gpu

      call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in,                             &
                  ZFH_Re_fspl_in, ZFH_Im_fspl_in,                             &
                  SH(1:NT_in), NT_in, M0, SMATCH, chi_, R, Z, 1, iflag )

      i       = nint( real(nR_-1,gpu)*(R(1)-Rmin_)/(Rmax_-Rmin_) ) + 1
      j       = nint( real(nZ_-1,gpu)*(Z(1)-Zmin_)/(Zmax_-Zmin_) ) + 1

!     N.B. vectors are in the same coordinate system for +|n| but the equations
!     are different.  Remember from Eqn. [2] for the MARS-F data:
!
!     BR = Re{BR}*cos( +|n|*PHI_I ) - Im{BR}*sin( +|n|*PHI_I ) [2]
!
!     and for the ITER filaments model
!
!     BR = Re{BR}*cos( +|n|*PHI_I ) + Im{BR}*sin( +|n|*PHI_I ),
!
!     hence the - sign on the Im parts of the MARS-F data:

      PHSE(1) = acos((B(i,j,1,1,1)*Re_F(i,j,1)-B(i,j,1,1,2)*Im_F(i,j,1))/     &
               (sqrt( B(i,j,1,1,2)**2+B(i,j,1,1,1)**2)*                       &
                sqrt( Im_F(i,j,1)**2+Re_F(i,j,1)**2)))*                       &
                sign( 1.0_gpu, (-B(i,j,1,1,2)*Re_F(i,j,1) -                   &
                                 B(i,j,1,1,1)*Im_F(i,j,1)) )
      PHSE(2) = acos((B(i,j,1,1,3)*Re_F(i,j,3)-B(i,j,1,1,4)*Im_F(i,j,3))/     &
               (sqrt( B(i,j,1,1,4)**2+B(i,j,1,1,3)**2)*                       &
                sqrt( Im_F(i,j,3)**2+Re_F(i,j,3)**2)))*                       &
                sign( 1.0_gpu, (-B(i,j,1,1,4)*Re_F(i,j,3) -                   &
                                 B(i,j,1,1,3)*Im_F(i,j,3)) )
      PHSE(3) = acos((B(i,j,1,1,5)*Re_F(i,j,2)-B(i,j,1,1,6)*Im_F(i,j,2))/     &
               (sqrt( B(i,j,1,1,6)**2+B(i,j,1,1,5)**2)*                       &
                sqrt( Im_F(i,j,2)**2+Re_F(i,j,2)**2)))*                       &
                sign( 1.0_gpu, (-B(i,j,1,1,6)*Re_F(i,j,2) -                   &
                                 B(i,j,1,1,5)*Im_F(i,j,2)) )

      BF(1) =  sqrt(Re_F(i,j,1)**2+Im_F(i,j,1)**2)
      BF(2) =  sqrt(Re_F(i,j,3)**2+Im_F(i,j,3)**2)
      BF(3) =  sqrt(Re_F(i,j,2)**2+Im_F(i,j,2)**2)
      BF(4) =  sqrt(BF(1)**2+BF(2)**2+BF(3)**2)

      if(BF(4) > ddat(1))then

!        Update tuning point:

         ddat(1)   = BF(4)
         ddat(2:4) = BF(1:3)
         ddat(5:7) = PHSE(1:3)
         ddat(8)   = BF(1)/sqrt(B(i,j,1,1,1)**2+B(i,j,1,1,2)**2)
         ddat(9)   = BF(2)/sqrt(B(i,j,1,1,3)**2+B(i,j,1,1,4)**2)
         ddat(10)  = BF(3)/sqrt(B(i,j,1,1,5)**2+B(i,j,1,1,6)**2)
         ddat(11)  = R(1)
         ddat(12)  = Z(1)
         i0 = i
         j0 = j
      endif

      enddo

!--- End Main GRID finding loop ---

      write(io(1),'(A33,1E16.8)')  ':mars_read : Peak field       : ',        &
                                     ddat(1)
      write(io(1),'(A33,3E16.8)')  ':mars_read : BR, BZ, BT [T]   : ',        &
                                     ddat(2:4)

!-   This is the phase correction on BR,BZ,BT required to match the Aalto field
!-   map at SMATCH and the highest field perturbation location on that flux
!-   surface.

      write(io(1),'(A33,3E16.8)')  ':mars_read : Phase offsets    : ',        &
                                     ddat(5:7)*180.0_gpu/                     &
                                     (nmde*pi)
      write(io(1),'(A33,3E16.8)')  ':mars_read : AMPL corrections : ',        &
                                     ddat(8:10)
      write(io(1),'(A33,2E16.8)')  ':mars_read : R, Z [m]         : ',        &
                                     ddat(11:12)
      write(io(1),'(A33,2I5)')     ':mars_read : Matching nodes   : ', i0, j0

      call display( bold, lun=io(1) )

!-    The suggested adjustment will be the largest out of BR,BZ,BT. One can
!-    force one component of the total field operturbation to be a best match
!-    if one overides the conditionals below. The PH1 adjustments should be
!-    added to the relevent component of PH1 to get a best fit.

      if ( maxval(ddat(2:4)) == ddat(2) )then
      write(io(1),*)               ':mars_read : Matching to BR'
      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjustment to PH1 : ',       &
                                     ddat(5)*180.0_gpu/(nmde*pi)

      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjusted IMUL     : ',       &
                                     IMUL(MATCH)* ddat(8)

      elseif ( maxval(ddat(2:4)) == ddat(3) )then
      write(io(1),*)               ':mars_read : Matching to BZ'
      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjustment to PH1 : ',       &
                                     ddat(6)*180.0_gpu/(nmde*pi)

      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjusted IMUL     : ',       &
                                     IMUL(MATCH)* ddat(9)

      else
      write(io(1),*)               ':mars_read : Matching to BT'
      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjustment to PH1 : ',       &
                                     ddat(7)*180.0_gpu/(nmde*pi)

      write(io(1),'(A34,1E16.8)')  ':mars_read : Adjusted IMUL     : ',       &
                                     IMUL(MATCH)* ddat(10)

      endif

      call display( norm, lun=io(1) )

#endif

#if defined (UMPF)

      write(io(1),*) ':mars_read : Expand onto 3D mesh....'

      if( nT_f==0 )nT_f = 37

      allocate( PHI_(nT_f) )
      allocate( BR0(nR_f,nZ_f,nT_f), BZ0(nR_f,nZ_f,nT_f), BT0(nR_f,nZ_f,nT_f) )

      PHI_ = (/(real(i,gpu),i=0,nT_f-1)/)*2.0_gpu*pi/real(nT_f-1,gpu)

!     B is for +|n| in ITER coordinate system, but nmde -s -ve here,
!     i.e. reconstruction has form Re{B}*cos(|n|*PHI_) - Im{B}*sin(|n|*PHI_)
!     == Re{B}*cos(-|n|*PHI_) + Im{B}*sin(-|n|*PHI_) with PHI_ in ITER 
!     coordinates.

      do i=1,nR_f
         do j=1,nZ_f
            BR0(i,j,1:nT_f) = B(i,j,1, 1, 1)*cos(PHI_*nmde) +                 &
                              B(i,j,1, 1, 2)*sin(PHI_*nmde)
            BZ0(i,j,1:nT_f) = B(i,j,1, 1, 3)*cos(PHI_*nmde) +                 &
                              B(i,j,1, 1, 4)*sin(PHI_*nmde)
            BT0(i,j,1:nT_f) = B(i,j,1, 1, 5)*cos(PHI_*nmde) +                 &
                              B(i,j,1, 1, 6)*sin(PHI_*nmde)
         enddo
      enddo

! -> dB_map.dat runs 0 to 2pi in PHI:

      open( unit=lun,file=TRIM(ADJUSTL(root))//'dB_map.dat_MARSEXP', form='unformatted' )

      STR = '3D Expansion of MARS-F data....'

      write( lun ) STR

      STR = 'Field is MARS-F origin'

      write( lun ) STR

      write( STR, '(A6,F12.6,2F16.6)' ) 'PH1 : ',PH1(1:3)

      write( lun ) STR

      STR = 'UNKNOWN'

      write( lun ) STR
      write( lun ) STR
      write( lun ) STR

      write( lun ) 1
      write( lun ) nR_f

#if defined (DBSINGLE)
      write( lun ) real(R0F_,single)
      write( lun ) real((R1F_-R0F_)/real(nR_f-1,gpu),single)
      write( lun ) nZ_f
      write( lun ) real(Z0F_,single)
      write( lun ) real((Z1F_-Z0F_)/real(nZ_f-1,gpu),single)
      write( lun ) nT_f
      write( lun ) real(BR0,single)
      write( lun ) real(BZ0,single)
      write( lun ) real(BT0,single)
#else
      write( lun ) R0F_
      write( lun ) (R1F_-R0F_)/real(nR_f-1,gpu)
      write( lun ) nZ_f
      write( lun ) Z0F_
      write( lun ) (Z1F_-Z0F_)/real(nZ_f-1,gpu)
      write( lun ) nT_f
      write( lun ) BR0
      write( lun ) BZ0
      write( lun ) BT0
#endif
      close( lun )

      write(io(1),*) ':mars_read : File dump complete.'

      deallocate( PHI_, BR0, BZ0, BT0 )

#endif

!     BPLASMA# file:

#ifndef COILROW
      write(io(1),*) ':mars_read : Write out file : '//TRIM(ADJUSTL(root))//'BPLASMA_n'//          &
                                                        TRIM(ADJUSTL(STR_N))

      open( unit=lun, file=TRIM(ADJUSTL(root))//'BPLASMA_n'//TRIM(ADJUSTL(STR_N)),                 &
            form='formatted', status='replace' )
#else
      write(io(1),*) ':mars_read : Write out file : '//TRIM(ADJUSTL(root))//'BPLASMA_n'//          &
                                                        TRIM(ADJUSTL(STR_N))//&
                                                   '_'//TRIM(ADJUSTL(STR_C))

      open( unit=lun, file=TRIM(ADJUSTL(root))//'BPLASMA_n'//TRIM(ADJUSTL(STR_N))//                &
                                   '_'//TRIM(ADJUSTL(STR_C)),                 &
            form='formatted', status='replace' )
#endif

            write( lun, '(A39,A2)' )                                          &
                            '**** Source mars_read.f90 **** : n = ',          &
                             TRIM(ADJUSTL(STR_N))
#ifndef COILROW
            write( lun, * ) 'Phase U,M,L : '
            write( lun, '(3E15.8)' ) PH1(1:3)
#else
            write( lun, * ) 'Phase Row '//TRIM(ADJUSTL(STR_C))//' :'
            write( lun, '(3E15.8)' ) PH1(COILROW)
#endif

#if defined (FXXYY)
            write( lun, '(2I9)' ) nR_f, nZ_f
#endif
            do iR=1,nR_f
               R_p = R0F_ + (R1F_-R0F_)*real(iR-1,gpu)/real(nR_f-1,gpu)
               do iZ=1,nZ_f
                  Z_p = Z0F_ + (Z1F_-Z0F_)*real(iZ-1,gpu)/real(nZ_f-1,gpu)

! To be consistent with old code/runs, LOCUST is going to flip BT, so flip
! it back here.  LOCUST will also use -|n| multiplied by -PHI:

#if defined (FXXYY)
                  write( lun, 200 ) R_p, ' ', Z_p, ' ',                       &
                                    B(iR,iZ,4,4,1), ' ',                      &
                                    B(iR,iZ,4,4,2), ' ',                      &
                                    B(iR,iZ,4,4,3), ' ',                      &
                                    B(iR,iZ,4,4,4), ' ',                      &
                                   -B(iR,iZ,4,4,5), ' ',                      &
                                   -B(iR,iZ,4,4,6), ' '
#else
                  write( lun, 200 ) R_p, ' ', Z_p, ' ',                       &
                                    B(iR,iZ,1,1,1), ' ',                      &
                                    B(iR,iZ,1,1,2), ' ',                      &
                                    B(iR,iZ,1,1,3), ' ',                      &
                                    B(iR,iZ,1,1,4), ' ',                      &
                                   -B(iR,iZ,1,1,5), ' ',                      &
                                   -B(iR,iZ,1,1,6), ' '
#endif
               enddo
            enddo

#if defined (FXXYY)
!           This part is extra to the traditional BPLASMA MARS-F format:

            do i=1,18
               write( lun, '(8E16.8)' ) FDER(1:nR_f,1:nZ_f,i)
            enddo
#endif

      close( unit=lun )

#if defined (SCHIDUMP)

      open ( unit=lun, file=TRIM(ADJUSTL(root))//'BPLASMA_s_chi.dat', form='unformatted' )
      write( lun ) NT*2, 361, nR_f, nZ_f  ! Grid dimensions
      write( lun ) SH                     ! S
      write( lun ) R0F_, R1F_, Z0F_, Z1F_ ! Grid size
      write( lun ) Rmag, Zmag             ! Magnetic axis
      write( lun ) B(1:nR_f,1:nZ_f,1,1,7) ! s[R,Z]
      write( lun ) B(1:nR_f,1:nZ_f,1,1,8) ! chi[R,Z]

#endif

      deallocate( B, FDER )

#if defined (SCHIDUMP)

!     Generate plot vs. s and chi (to check continuity across s=1):
      
      allocate( B( 6, NT*2, 361, 1, 1 ) )
      
      do j=1,NT*2
         do k=1,361
            s_   = SH(j)
            chi_ = real(k-1,gpu)*2.0_gpu*pi/360.0_gpu
            
            if( s_ <= 1.0_gpu )then
                call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,   &
                            ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_,       &
                            chi_, R, Z, 0, iflag )
            else
                call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,                  &
                            ZFH_Re_fspl_out,ZFH_Im_fspl_out,                  &
                            SH(NT_in:NT_out), NT_out, M0, s_,                 &
                            chi_, R, Z, 0, iflag )
            endif

            if( s_ <= 1.0_gpu )then
                call B_F( B1_Re_fspl_in, B1_Im_fspl_in, S (1:N1_in),          &
                          N1_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B1,         &
                          0, iflag )                     
                call B_F( B2_Re_fspl_in, B2_Im_fspl_in, SM(1:N2_in),          &
                          N2_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B2,         &
                          0, iflag )
                call B_F( B3_Re_fspl_in, B3_Im_fspl_in, SM(1:N3_in),          &
                          N3_in, M1, M2, nmde, s_, chi_, 0.0_gpu, B3,         &
                          0, iflag )
           else
                call B_F( B1_Re_fspl_out, B1_Im_fspl_out,                     &
                          S (N1_in:NT), N1_out, M1, M2, nmde, s_, chi_,       &
                          0.0_gpu, B1, 0, iflag )
                call B_F( B2_Re_fspl_out, B2_Im_fspl_out,                     &
                          SM(N2_in+1:NT), N2_out, M1, M2, nmde, s_, chi_,     &
                          0.0_gpu, B2, 0, iflag )
                call B_F( B3_Re_fspl_out, B3_Im_fspl_out,                     &
                          SM(N3_in+1:NT), N3_out, M1, M2, nmde, s_, chi_,     &
                          0.0_gpu, B3, 0, iflag )
            endif
                   
            JAC       = (R(2)*Z(3) - R(3)*Z(2))*R(1)
                   
            B( 1,j,k,1,1) = R(2) ! (B1(1)*R(2)+B2(1)*R(3))/JAC
            B( 2,j,k,1,1) = Z(3) ! (B1(2)*R(2)+B2(2)*R(3))/JAC
            B( 3,j,k,1,1) = R(3) ! (B1(1)*Z(2)+B2(1)*Z(3))/JAC
            B( 4,j,k,1,1) = Z(2) ! (B1(2)*Z(2)+B2(2)*Z(3))/JAC
            B( 5,j,k,1,1) = R(1) !  B3(1)*R(1)/JAC
            B( 6,j,k,1,1) =  JAC ! B3(2)*R(1)/JAC

         enddo
      enddo

      write( lun ) B
      close( unit=lun )

      deallocate( B )

#endif

#if defined (SCHIPLOT)

!     Generate plot vs. s and chi (to check continuity across s=1):
      
      allocate( B( 10, NT*2, 361, 1, 1 ) )
      
      do j=1,NT*2
         do k=1,361
            s_   = SH(j)
            chi_ = real(k-1,gpu)*2.0_gpu*pi/360.0_gpu
            
            if( s_ <= 1.0_gpu )then
                call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,   &
                            ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, s_,       &
                            chi_, R, Z, 0, iflag )
            else
                call rspos( RFH_Re_fspl_out,RFH_Im_fspl_out,                  &
                            ZFH_Re_fspl_out,ZFH_Im_fspl_out,                  &
                            SH(NT_in:NT_out), NT_out, M0, s_,                 &
                            chi_, R, Z, 0, iflag )
            endif

            if( s_ <= 1.0_gpu )then
                call B_F( B1_Re_fspl_in, B1_Im_fspl_in, S (1:N1_in),          &
                          N1_in, M1, M2, n, s_, chi_, 0.0_gpu, B1,            &
                          0, iflag )                     
                call B_F( B2_Re_fspl_in, B2_Im_fspl_in, SM(1:N2_in),          &
                          N2_in, M1, M2, n, s_, chi_, 0.0_gpu, B2,            &
                          0, iflag )
                call B_F( B3_Re_fspl_in, B3_Im_fspl_in, SM(1:N3_in),          &
                          N3_in, M1, M2, n, s_, chi_, 0.0_gpu, B3,            &
                          0, iflag )
           else
                call B_F( B1_Re_fspl_out, B1_Im_fspl_out,                     &
                          S (N1_in:NT), N1_out, M1, M2, n, s_, chi_,          &
                          0.0_gpu, B1, 0, iflag )
                call B_F( B2_Re_fspl_out, B2_Im_fspl_out,                     &
                          SM(N2_in+1:NT), N2_out, M1, M2, n, s_, chi_,        &
                          0.0_gpu, B2, 0, iflag )
                call B_F( B3_Re_fspl_out, B3_Im_fspl_out,                     &
                          SM(N3_in+1:NT), N3_out, M1, M2, n, s_, chi_,        &
                          0.0_gpu, B3, 0, iflag )
            endif
                   
            JAC       = (R(2)*Z(3) - R(3)*Z(2))*R(1)
                   
            B( 1,j,k,1,1) = (B1(1)*R(2)+B2(1)*R(3))/JAC
            B( 2,j,k,1,1) = (B1(1)*Z(2)+B2(1)*Z(3))/JAC
            B( 3,j,k,1,1) =  B3(1)*R(1)/JAC
            B( 4,j,k,1,1) =  R(1)
            B( 5,j,k,1,1) =  Z(1)
            B( 6,j,k,1,1) =  B1(1)
            B( 7,j,k,1,1) =  Z(2)
            B( 8,j,k,1,1) =  B2(1)
            B( 9,j,k,1,1) =  Z(3)
            B(10,j,k,1,1) =  B3(1) ! JAC
         enddo
      enddo
      
      open ( unit=lun, file=TRIM(ADJUSTL(root))//'BSCHI.dat', form='unformatted' )
      write( lun ) NT*2, 361
      write( lun ) SH
      write( lun ) S
      write( lun ) q
      write( lun ) B
      close( unit=lun )

#endif

#if defined (DIAGNOSE)

    write(io(1),*) ':mars_read : DIAGNOSE'

    call B_F( B1_Re_fspl_in, B1_Im_fspl_in, S (1:N1_in),                      &
              N1_in, M1, M2, nmde, 0.497511d0, 0.0d0, 0.0_gpu, B1,            &
              -1, iflag )                     
    call B_F( B2_Re_fspl_in, B2_Im_fspl_in, SM(1:N2_in),                      &
              N2_in, M1, M2, nmde, 0.497511d0, 0.0d0, 0.0_gpu, B2,            &
              0, iflag )
    call B_F( B3_Re_fspl_in, B3_Im_fspl_in, SM(1:N3_in),                      &
              N3_in, M1, M2, nmde, 0.497511d0, 0.0d0, 0.0_gpu, B3,            &
              0, iflag )

    call rspos( RFH_Re_fspl_in, RFH_Im_fspl_in, ZFH_Re_fspl_in,               &
                ZFH_Im_fspl_in, SH(1:NT_in), NT_in, M0, 0.497511d0,           &
                0.0d0, R, Z, 0, iflag )

    JAC       = (R(2)*Z(3) - R(3)*Z(2))*R(1)

    write(io(1),*) ':mars_read : s       [-] :', 0.497511_gpu
    write(io(1),*) ':mars_read : chi     [-] :', 0.0_gpu
    write(io(1),*) ':mars_read : R       [m] :', R(1)
    write(io(1),*) ':mars_read : Z       [m] :', Z(1)
    write(io(1),*) ':mars_read : Re{B1}  [?] :', B1(1)
    write(io(1),*) ':mars_read : Im{B1}  [?] :', B1(2)
    write(io(1),*) ':mars_read : Re{B2}  [?] :', B2(1)
    write(io(1),*) ':mars_read : Im{B2}  [?] :', B2(2)
    write(io(1),*) ':mars_read : Re{B3}  [?] :', B3(1)
    write(io(1),*) ':mars_read : Im{B3}  [?] :', B3(2)
    write(io(1),*) ':mars_read : dR/ds   [m] :', R(2)
    write(io(1),*) ':mars_read : dR/dchi [m] :', R(3)
    write(io(1),*) ':mars_read : dZ/ds   [m] :', Z(2)
    write(io(1),*) ':mars_read : dZ/dchi [m] :', Z(3)
    write(io(1),*) ':mars_read : JAC         :', JAC

    write(io(1),*) ':mars_read : BR      [T] :', (B1(1)*R(2)+B2(1)*R(3))/JAC
    write(io(1),*) ':mars_read : BZ      [T] :', (B1(1)*Z(2)+B2(1)*Z(3))/JAC
    write(io(1),*) ':mars_read : BT      [T] :',  B3(1)*R(1)/JAC
    write(io(1),*) ':mars_read : |B|     [T] :', sqrt( (B1(1)*R(2)+           &
                                                        B2(1)*R(3))**2 +      &
                                                       (B1(1)*Z(2)+           &
                                                        B2(1)*Z(3))**2 +      &
                                                       (B3(1)*R(1))**2 )/JAC
#endif

      write(io(1),*) ':mars_read : File dumps complete.'

200  format(8(E15.8,A1))

end program

