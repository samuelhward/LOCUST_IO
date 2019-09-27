program dB_map_collate
!
!     Compbine two 3D field maps (i.e. to assimilate multiple n's).
!
!     Compile with:
!
!     pgfortran -o dB_map_collate prec_mod.f90 dB_map_collate.f90 -DBP 
!                  -Mpreprocess
!
      use prec_mod

      implicit none

      character( len=100 )             :: STR1, STR2, STR3, STR4, STR5, STR6
      integer                          :: lun   = 5
      integer                          :: S, nR_f, nZ_f, nT_f, i, n
      real( gpu )                      :: R0F_, dRF_, Z0F_, dZF_
      real( gpu ),        allocatable  :: BR0(:,:,:), BZ0(:,:,:), BT0(:,:,:)
      real( gpu ),        allocatable  :: DAT(:,:,:)
      character( len=2 ), dimension(4) :: TAIL = ['_1','_2','_3','_4']

      write(io(1),*) ':dB_map_collate : Combine 3D field maps....'

      read(*,*) n

      write(io(1),*) ':dB_map_collate : # harmonics : ', n

      do i=1,n

      write(io(1),*) ':dB_map_collate : Load file : '//'dB_map.dat'//TAIL(i)

      open( unit=lun,file='dB_map.dat'//TAIL(i), form='unformatted' )

      read( lun ) STR1
      read( lun ) STR2
      read( lun ) STR3
      read( lun ) STR4
      read( lun ) STR5
      read( lun ) STR6

      read( lun ) S
      read( lun ) nR_f

      read( lun ) R0F_
      read( lun ) dRF_
      read( lun ) nZ_f
      read( lun ) Z0F_
      read( lun ) dZF_
      read( lun ) nT_f

      if( i==1 )then

          allocate( BR0(nR_f,nZ_f,nT_f) )
          allocate( BZ0(nR_f,nZ_f,nT_f) )
          allocate( BT0(nR_f,nZ_f,nT_f) )
          allocate( DAT(nR_f,nZ_f,nT_f) )

          BR0(:,:,:) = 0.0_gpu
          BZ0(:,:,:) = 0.0_gpu
          BT0(:,:,:) = 0.0_gpu

      endif

      read( lun ) DAT

      BR0(1:nR_f,1:nZ_f,1:nT_f) = BR0(1:nR_f,1:nZ_f,1:nT_f) +                 &
                                  DAT(1:nR_f,1:nZ_f,1:nT_f)

      read( lun ) DAT

      BZ0(1:nR_f,1:nZ_f,1:nT_f) = BZ0(1:nR_f,1:nZ_f,1:nT_f) +                 &
                                  DAT(1:nR_f,1:nZ_f,1:nT_f)

      read( lun ) DAT

      BT0(1:nR_f,1:nZ_f,1:nT_f) = BT0(1:nR_f,1:nZ_f,1:nT_f) +                 &
                                  DAT(1:nR_f,1:nZ_f,1:nT_f)

      close( lun )

      enddo

      write(io(1),*) ':dB_map_collate : Write out combined map....'

      open( unit=lun,file='dB_map.dat_assimilated', form='unformatted' )

      STR1 = 'COMBINED HARMONICS from dB_map_collate.f90'

      write( lun ) STR1
      write( lun ) STR2
      write( lun ) STR3
      write( lun ) STR4
      write( lun ) STR5
      write( lun ) STR6

      write( lun ) S
      write( lun ) nR_f

      write( lun ) R0F_
      write( lun ) dRF_
      write( lun ) nZ_f
      write( lun ) Z0F_
      write( lun ) dZF_
      write( lun ) nT_f
      write( lun ) BR0
      write( lun ) BZ0
      write( lun ) BT0

      close( lun )

      deallocate( BR0, BZ0, BT0, DAT )

end program
