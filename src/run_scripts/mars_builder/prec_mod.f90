module prec_mod
!------------------------------------------------------------------------------
!
!     MASTER parameter module for LOCUST-GPU.
!
! *** Author:
!
!     Rob Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
!
! *** Creation Date:
!
!     31/03/2010 - present
!
!------------------------------------------------------------------------------

#if defined (MOVIE)
#ifndef NTMAX
# define NTMAX 100000
#endif
#ifndef STATIC
# define STATIC
#endif
#ifndef TEST
#endif
#ifndef NTMIN
# define NTMIN 0
#endif
#endif
#ifndef TOKAMAK
# define TOKAMAK 1
#endif

implicit none

! I/O

integer, dimension(40), parameter :: io   = [ 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
                                              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
                                              6, 6, 6, 6, 6, 6, 6, 6, 6, 6, &
                                              6, 6, 6, 6, 6, 6, 6, 6, 6, 6 ]
integer, dimension(40), parameter :: incl = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
                                              1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ]

character(len=13),      parameter :: path_out = '/home/rakers/'
character(len=48)                 :: runfile
!character(len=100)               :: file_tet = '/home/rakers/&
!                                              &full_model_40mm.dat_THICKELM'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                              &test_cf.attilamesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                              &test_cf_3.mesh.inp_ANSYS_CR'
character(len=100)               :: file_tet = '/home/rakers/&
                                              &ITER_meshC.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                              &MAST-U_VII.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                    &MAST-U_VII_cleaned_v5.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                    &MAST_Final_XIII_v4_bc2.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                    &MAST-U_VII_cleaned_v6.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                    &MAST-U_VII_cleaned_v7.mesh.inp_ANSYS'
!character(len=100)               :: file_tet = '/home/rakers/&
!                                    &MAST-U_VII_cleaned_v9.mesh.inp_ANSYS'
                                                                             
! Numerical precision etc.

integer,       parameter :: single = selected_real_kind( precision( 1.0e0 ) )
integer,       parameter :: double = selected_real_kind( precision( 1.0d0 ) )

#if defined (BP)
 integer,      parameter :: gpu    = double
#else
 integer,      parameter :: gpu    = single
#endif

integer,       parameter :: long   = selected_int_kind ( 18 )

! Physical constants:

real( gpu ),   parameter :: echg   = 1.60217646e-19_gpu
real( gpu ),   parameter :: echg_  = 1.60217646e+00_gpu
real( gpu ),   parameter :: c      = 2.9979246e+08_gpu
real( gpu ),   parameter :: pi     = 3.1415926535897932_gpu
real(double),  parameter :: dpi    = 3.1415926535897932_double
real( gpu ),   parameter :: mp     = 1.6725966e-27_gpu
real( gpu ),   parameter :: mp_    = 1.6725966e+00_gpu
real( gpu ),   parameter :: me     = 9.10938215e-31_gpu
real( gpu ),   parameter :: An     = 1.0013783_gpu
real( gpu ),   parameter :: AD     = 1.9990024_gpu
real( gpu ),   parameter :: AT     = 2.9937181_gpu
real( gpu ),   parameter :: AC12   = 12.0_gpu
real( gpu ),   parameter :: AHe3   = 2.9931529_gpu
real( gpu ),   parameter :: AHe4   = 4.0015062_gpu
real( gpu ),   parameter :: Ae     = 1.0_gpu/1836.15266_gpu ! me-/mp+
real( gpu ),   parameter :: eps0   = 8.85418782e-12_gpu
real( gpu ),   parameter :: eps0_  = 8.85418782e+00_gpu
real( gpu ),   parameter :: sig_St = 5.670367e-08_gpu
real( gpu ),   parameter :: T_0_K  = 273.15_gpu

real( gpu )              :: Pdep   = 33.0e6_gpu    ! Injected power [W]

! CUDA parameters:

integer,       parameter :: threadsPerBlock = 32 ! 32
integer,       parameter :: blocksPerGrid   = 1024 ! 1024 ! 128  ! 128 ! 1024 ! 1024 ! 1024 ! threadsPerBlock *
                                                   ! blocksPerGrid must =
                                                   ! N * 4096

! Mersenne Twister needs to know # gpus:

integer                  :: ngpu   = 0             ! # available matching GPUs

!------------------------------------------------------------------------------

! Used by neutrons:

!integer,      parameter :: stype  = 1  ! Final state ptcl. type (0=n,1=p)
!integer,      parameter :: itype  = 1  ! Allow differential cross section asym?
!integer,      parameter :: nsamp  = 1  ! # neutrons per calcn.
!integer                 :: mran        ! # RNG points (neutrons)
!real( gpu ),  parameter :: fudge_fac1 = 1.0_gpu ! Must be 1.0 for neutrons

! Used by Wrapper

!-> DD fusion product simulation:

!integer,      parameter :: isym   = 1  ! Simulation type (0=fast ion, 1=p)      
!integer,      parameter :: iterm  = 0  ! Terminate if ptcl. leaves plasma etc.
!integer                 :: niter  = 20  ! # iterations for isym=1 simulation
!integer,      parameter :: nchunks= 50           ! # particle chunks
!real( gpu )             :: dt0    = 1.0e-9_gpu ! t-step   [s]
!real( gpu ),  parameter :: tstp   = 1.0e-9_gpu ! Max intgn step size        [s]
!real( gpu ),  parameter :: Ab    = 1.0_gpu ! AD
!real( gpu ),  parameter :: Ai    = AD

!-> D NBI: ------

integer,      parameter :: ihag   = 0        ! Generate HAGIS test DFn.                 [-]
integer,      parameter :: filtyp = 1        ! Filtering algorithm for -DFNLOAD         [-]
integer,      parameter :: imap   = 0        ! Generate orbit topology map              [-]
integer,      parameter :: isym   = 0        ! Simulation type   (0=fast ion, 1=p)      [-]
integer,      parameter :: iterm  = 0        ! Terminate if ptcl. leaves plasma         [-]
integer,      parameter :: ithrm  = 0        ! Fast ion (0) or thermal (1) mode         [-]
integer                 :: niter  = 64      ! 280 ! 24 ! 16       ! # iterations for isym=1 simulation       [-]
integer,      parameter :: i3dr   =-1        ! MARS-F BT flip factor                    [-]
integer                 :: nbuf   = 0        ! Include a slow ptcl. cycle               [-]
integer,      parameter :: nscan  = 10000    ! # points for topology calc               [-]
integer,      parameter :: nchunks= 2 ! 20 ! 4        ! # particle chunks                        [-]
integer,      parameter :: ncon   = 30       ! # concentric layers for tet search       [-]
integer(long),parameter :: nphc   = 11_long   ! 7_long   ! # levels in -DSPLIT split cache          [-]

#if defined (POINGYR)
integer,      parameter :: npoin  = 400000
integer,      parameter :: npois  = 10000
#else
integer,      parameter :: npoin  = 200
integer,      parameter :: npois  = 10
#endif

real( gpu ),  parameter :: PSI_rej= 0.70_gpu ! Inner PSI rejection for -DPOINCARE       [-]
real( gpu ),  parameter :: PSI_max= 1.00_gpu ! Max s for Poincare leakage               [-]
integer,      parameter :: ngyr   = 200000   ! Steps before return (Poincare)           [-]

! PFC data:

! Components to strip out:

integer,   dimension(1) :: istrip=[-1]            ! Component IDs to change to vacuum   [-]
integer,   dimension(1) :: ipfmsk=[ 2]            !       -"-     to reject from full model

! FILD instrumentation:

integer,   dimension(1) :: ifld = [-1]         ! PFCs to output ptcl. data           [-]
integer,   dimension(1) :: jfld = [-1]      ! PFC grouping for -DPFCF             [-]
real( gpu ),  parameter :: amin = 5.0e-06_gpu     ! Triangle refinement limit         [m^2]
integer                 :: ntet                   ! # tetrahedra in model               [-]
integer                 :: ntri                   ! Triangle grid "dimension"           [-]
integer                 :: ntri_c                 ! # triangles in entire model         [-]
integer                 :: ncmp                   ! # components in model               [-]
integer                 :: ndes                   ! # nodes in tetrahedral mesh         [-]
integer                 :: nhit                   ! # PFC hits                          [-]

#if defined (TEST)
real( gpu )             :: dt0    = 2.0e-9_gpu    ! t-step requested (2e-08)            [s]
real( gpu ),  parameter :: tstp   = 2.0e-9_gpu    ! Max intgn step size                 [s]
real( gpu ),  parameter :: tdsp   = 2.0e-05_gpu   ! Diagnostic display period           [s]
real( gpu ),  parameter :: V_0    = 9789.5934_gpu ! Test mode injection vel           [m/s]
#else
integer,      parameter :: nbor   = 50            ! # Boris steps per dt0               [-]
real( gpu )             :: dt0    = 1.0e-07_gpu ! 1.0e-07_gpu   ! t-step requested (2e-08)            [s]
real( gpu ),  parameter :: tstp   = 2.0e-08_gpu ! 1.0e-08_gpu ! 1.0e-08_gpu   ! 2.0e-08_gpu - GUIDE3 time sub-step  [s]
real( gpu ),  parameter :: tstp_0 = 5.0e-10_gpu   ! Min step for adaptive               [s]
real( gpu ),  parameter :: tdsp   = 1024.0e-07_gpu! Caching/compression period          [s]  
real( gpu ),  parameter :: t_dspl = 2048.0e-07_gpu! Diagnostic display period           [s]
real( gpu ),  parameter :: V_0    = 2.5e6_gpu     ! Test mode injection vel           [m/s]
#endif

real( gpu ),  parameter :: dtol   = 1.0e-04_gpu   ! TET cell tracking tolerance         [m]
real( gpu ),  parameter :: dLPP   = 0.2e-01_gpu   ! Poincare plot spatial step          [m]
real( gpu ),  parameter :: dLPP_0 = 1.0e-06_gpu   ! Min spatial step                    [m]




integer,      parameter :: iscan  = 0             ! Scan mode                 (0=off, 1=on)
integer,      parameter :: iaddth = 1             ! Add on thermal ions       (1=on, 0=off)
integer,      parameter :: icoll  = 1             ! Collisions                (1=on, 0=off)
integer,      parameter :: iscat  = 1             ! Scattering                (1=on, 0=off)
integer,      parameter :: idiff  = 1             ! Energy diffusion          (1=on, 0=off)
integer,      parameter :: iloss  = 0             ! CX losses                 (1=on, 0=off)
integer,      parameter :: ncyc   = 5             ! # cycles for integrator #1 (imccah) [-]
integer,      parameter :: norb   = 10            ! # poloidal transits to find LCFS    [-]
real( gpu ),  parameter :: dorb   = 90.0_gpu      ! Toroidal transit pi*norb/dorb

real( gpu ),  parameter :: bthr_v = 0.02_gpu      ! Bounce V|| threshold                [-]
real( gpu ),  parameter :: mcyc   = 1.0_gpu       ! * tgyro(GC) b.t.  thresh            [-]
real( gpu ),  parameter :: lcyc   = 2.0_gpu       ! * tgyro(GC) B_r=0 thresh            [-]
real( gpu ),  parameter :: tcut   = 0.1_gpu       ! Min step cut for GJRK45             [s]
integer,      parameter :: nran   = 5             ! # RNG points per step               [-]
integer,      parameter :: npnt   = 60            ! # points per gyration               [-]
real( gpu ),  parameter :: RKmult = 1.0_gpu       ! Relax R.K. error (0.1_gpu)          [-]
real( gpu ),  parameter :: pcore  = 1.0_gpu       ! Jacobian core enhancement           [-]
real( gpu ),  parameter :: R0C    = 3.5_gpu
real( gpu ),  parameter :: R1C    = 8.9_gpu
real( gpu ),  parameter :: Z0C    =-5.0_gpu
real( gpu ),  parameter :: Z1C    =+5.0_gpu

integer,      parameter :: njac   = 7             ! # RNG points for Jac calc           [-]
integer,      parameter :: mjac   = blocksPerGrid*threadsPerBlock*10
integer,      parameter :: mgen   = 10000         ! Jacobian display interval           [-]
integer,      parameter :: ngen   = mgen*1        ! J inventory multiplier              [-]
integer,      parameter :: mhag   = blocksPerGrid*threadsPerBlock*10
integer,      parameter :: mrun   = 10000   
integer,      parameter :: nrun   = mrun*1        ! 200 default

! TEST mode:

#if defined (MOVIE)
real( gpu ),  parameter :: siglg  = 0.02_gpu      ! r.m.s. width of test src            [-]
#else
real( gpu ),  parameter :: siglg  = 0.02_gpu      ! r.m.s. width of test src            [-]
#endif

! Used by neutrons:

integer,      parameter :: stype  = 1  ! Final state ptcl. type (0=n,1=p)
integer,      parameter :: atype  = 1  ! Allow differential cross section asym?
integer,      parameter :: nsamp  = 1  ! # neutrons per calcn.
integer                 :: mran        ! # RNG points (neutrons)

! Artificial scale factors (for testing):

real( gpu ),  parameter :: fudge_fac1 = 1.0_gpu ! x factor for ion velocity
real( gpu ),  parameter :: fudge_fac2 = 1.0_gpu ! x factor for BR & BZ
real( gpu ),  parameter :: fudge_fac3 = 1.0_gpu ! x factor for BT

!------------------------------------------------------------------------------

! Any number of bulk ion species including impurities can be added here:

integer,      parameter                  :: nion = 3
real( gpu ),  parameter, dimension(nion) :: fi   = [ 0.42385000_gpu,          &
                                                     0.42385000_gpu,          &
                                                     0.02538333_gpu ]  ! Ion fractions              [-]
real( gpu ),  parameter, dimension(nion) :: Ai   = [ AD,AT,AC12 ]      ! Background ion mass      [amu]
real( gpu ),  parameter, dimension(nion) :: Mi   =   Ai*mp             ! Bulk ion mass             [kg]
real( gpu ),  parameter, dimension(nion) :: Zi   = [+1.0_gpu,                 &
                                                    +1.0_gpu,                 &
                                                    +6.0_gpu]          ! Bulk ion Z                 [-]

real( gpu ),  parameter                  :: Ab   = AD                  ! Fast ion mass            [amu]
real( gpu ),  parameter                  :: Zb   = +1.0_gpu            ! Trace ptcl. Z              [-]

real( gpu ),  parameter :: Mb     = Ab*mp       ! Trace ptcl. mass          [kg]
real( gpu ),  parameter :: Mb_    = Ab*mp_      ! Trace ptcl. mass          [kg]

real( gpu ),  parameter :: gam_t  = 1.0e2_gpu*Zb**2*echg_**4/                  &
                                            (Mb_**2*4.0_gpu*pi*eps0_**2)

real( gpu ),  parameter :: gnormh = Zb*echg/                                   &
                                    (2*pi*Mb)   !  tgyro_GC normalization
real( gpu ),  parameter :: pnormh = npnt*gnormh !  t-step   normalization
real( gpu ),  parameter :: eomh   = Zb*echg/Mb  !  e/m

! DD Fusion Q values [J]:

real( gpu ),  parameter :: Q_DDn  = (AD+AD-AHe3-An     )*mp*(c**2)
real( gpu ),  parameter :: Q_DDp  = (AD+AD-AT  -1.0_gpu)*mp*(c**2)

! General numerics:

#if defined   (AR_GAUGE)
integer,      parameter :: nc_3D  = 16*6+6*4
#elif defined (AZ_GAUGE)
integer,      parameter :: nc_3D  = 16*6+6*4
#else
integer,      parameter :: nc_3D  = 16*6
#endif
integer,      parameter :: nc_2D  = 16
#ifndef NOCOMPRESS
integer,      parameter :: npt    = 37
#else
integer,      parameter :: npt    = 52
#endif
integer,      parameter :: npt_   = 11
integer,      parameter :: n1d    = 6

! Grid cell dimensions for wrapper and neutrons:

real( gpu ),  parameter :: tau_ang= 0.35_gpu*pi
real( gpu ),  parameter :: POWPSI = 2.0_gpu ! Set to 2.0_gpu to make contrours and H vs. sqrt(PSIn)
integer,      parameter :: nEQ_R  = 257 ! 150 ! 256 ! 301 ! 257 ! 150
integer,      parameter :: nEQ_Z  = 513 ! 200 ! 257 ! 401 ! 513 ! 200
integer,      parameter :: nG     = 400  ! # R,Z master grid cells
integer,      parameter :: nF     = 49   ! # R,Z cell boundary
integer,      parameter :: nLF    = 181  ! # Vphi/V cell boundaries
integer,      parameter :: nPP    = 345 
integer,      parameter :: nVF    = 45   ! # V cell boundaries
integer,      parameter :: nEF    = 45   ! # E cell boundaries
integer,      parameter :: nPF    = 2    ! # Poloidal gyro-phase
integer,      parameter :: nMU    = 360  !
integer,      parameter :: nPSIF  = 101   !
integer,      parameter :: nPOLF  = 361  ! Advisable for this to be odd
integer,      parameter :: nPP_f  = 345  !
integer,      parameter :: nMU_f  = 360  !
integer,      parameter :: nEF_f  = 45   !
integer,      parameter :: npn    = 20   ! Cell granularity for IDFTYP=4
integer,      parameter :: nMs    = 1000

#if defined (MOVIE)
integer,      parameter :: nMV    = NTMAX/nMs
#else
integer,      parameter :: nMV    = 1
#endif

integer,      parameter :: nxP    = 1080 !
integer,      parameter :: nyP    = 2000 !
integer,      parameter :: ntP    = 240  ! 
integer,      parameter :: npln   = 720  ! 
real( gpu ),  parameter :: x0P    = 3.5_gpu
real( gpu ),  parameter :: x1P    = 8.9_gpu
real( gpu ),  parameter :: y0P    =-5.0_gpu
real( gpu ),  parameter :: y1P    =+5.0_gpu

integer,      parameter :: nRn    = 30   ! # R   cell boundaries (neutron grid)
integer,      parameter :: nZn    = 30   ! # Z   cell boundaries (neutron grid)

integer,      parameter :: iRS    = 14   ! First radial   cell of DFn. grid.
integer,      parameter :: iZS    = 11   ! First vertical cell of DFn. grid.

real( gpu ),  parameter :: V1n    = 2.40e7_gpu*1.1_gpu*fudge_fac1
real( gpu ),  parameter :: V0n    = 1.95e7_gpu*1.1_gpu*fudge_fac1
real( gpu ),  parameter :: dVn    = (V1n-V0n)/(nVF -1)

integer,      parameter :: nPSI   = 1000       ! # flux surface for Te etc.
integer,      parameter :: nIMG   = 65         ! Fusion product cam resolution
real( gpu ),  parameter :: FGmax  = 150.0_gpu  ! (150) Goosing factor @ Vscl m/s
real( gpu ),  parameter :: FGstp  = 1.0_gpu    ! Max increase in FG per orbit
real( gpu ),  parameter :: goocn1 = 10.00_gpu  ! # transits per V bin traverse
real( gpu ),  parameter :: goocn2 = 50.0_gpu   ! 5.0_gpu    ! # transits per pi/2 scatter
real( gpu ),  parameter :: goocn3 = 50.0_gpu   ! # transits per CX loss
real( gpu ),  parameter :: osafe  = 0.0_gpu    ! # initial FG=1 crossings
integer,      parameter :: ncom   = 0          ! Compression threshold
real( gpu ),  parameter :: Cf     = 1.0_gpu    ! 30 Collision operator split
real( gpu ),  parameter :: ICf    = 1.0_gpu/Cf ! 1/Cf
real( gpu ),  parameter :: mu_cut = 2.0_gpu    ! Cut off for HAGIS test DFn.

! Master grid boundary (data are located at cell nodes):

integer,      parameter :: nchx   = 30
integer,      parameter :: nchy   = 30
integer,      parameter :: ieq_L  = 5
integer,      parameter :: ieq_U  = 17
integer,      parameter :: jeq_L  = 14
integer,      parameter :: jeq_U  = 15

! Dfn. storage grid boundary (data are located at cell centres):

! ASCOT 3D tricubic spline grid:
!real( gpu ),  parameter :: R0F    = 3.52_gpu ! Assumption now is that 
!                                             ! collisions.dat is matched to 
!                                             ! this grid, res nG
!real( gpu ),  parameter :: R1F    = 8.98_gpu
!real( gpu ),  parameter :: Z0F    =-5.01_gpu
!real( gpu ),  parameter :: Z1F    = 5.01_gpu

! Old style BPLASMA files:

!real( gpu ),  parameter :: R0F    = 3.4_gpu
!real( gpu ),  parameter :: R1F    = 9.1_gpu 
!real( gpu ),  parameter :: Z0F    =-5.2_gpu
!real( gpu ),  parameter :: Z1F    = 5.2_gpu

! Settings best matched to ITER void:

#if (TOKAMAK==1)
real( gpu ),  parameter :: R0F    = 4.0_gpu
real( gpu ),  parameter :: R1F    = 8.6_gpu
real( gpu ),  parameter :: Z0F    =-4.8_gpu
real( gpu ),  parameter :: Z1F    = 4.8_gpu
#else
#ifndef ASCOTAUG
real( gpu ),  parameter :: R0F    = 0.50_gpu ! Main D.Fn. grid
real( gpu ),  parameter :: R1F    = 2.50_gpu ! 
real( gpu ),  parameter :: Z0F    =-1.50_gpu ! 
real( gpu ),  parameter :: Z1F    = 1.50_gpu !
#else
real( gpu ),  parameter :: R0F    = 0.75_gpu ! Main D.Fn. grid
real( gpu ),  parameter :: R1F    = 2.67_gpu ! 
real( gpu ),  parameter :: Z0F    =-1.504_gpu ! 
real( gpu ),  parameter :: Z1F    = 1.504_gpu !
#endif
#endif
real( gpu ),  parameter :: dRG    = (R1F-R0F)/real(nG-1,gpu)
real( gpu ),  parameter :: dZG    = (Z1F-Z0F)/real(nG-1,gpu)

real( gpu ),  parameter :: P0F    =-pi*(1.0_gpu+(1.0_gpu/real(nPF-1, gpu)))
real( gpu ),  parameter :: P1F    = P0F + 2*pi
real( gpu ),  parameter :: dPF    = (P1F-P0F) /real( nPF  -1, gpu)
real( gpu ),  parameter :: dRF    = (R1F-R0F) /real( nF   -1, gpu)
real( gpu ),  parameter :: dZF    = (Z1F-Z0F) /real( nF   -1, gpu)
real( gpu ),  parameter :: dLF    = 2.0_gpu   /real( nLF  -1, gpu)
real( gpu ),  parameter :: dPSIF  = 1.0_gpu   /real( nPSIF-1, gpu)
real( gpu ),  parameter :: dPOLF  = 2.0_gpu*pi/real( nPOLF-1, gpu)
real( gpu ),  parameter :: dOM    = 4.0_gpu*pi/real((nPF-1)*(nLF-1), gpu)

real( gpu ),  parameter :: Rcam   = 1.40_gpu
real( gpu ),  parameter :: Zcam   = 0.00_gpu
real( gpu ),  parameter :: dZcam  = 0.01_gpu
real( gpu )             :: sfac

! Accuracy graph:

real( double ), parameter :: E_0_t   = 70.0e3_double*real(echg,double)
real( double ), parameter :: E_c_t   = 10.0e3_double*real(echg,double)
real( double ), parameter :: dE_t    = 10.0e3_double*real(echg,double)

end module prec_mod
