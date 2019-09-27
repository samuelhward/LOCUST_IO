module disp_mod

character (len=4) :: norm = 'norm'
character (len=4) :: bold = 'bold'

contains

subroutine display(font,lun)
!------------------------------------------------------------------------------
!
! *** Description:
!
!     Change screen display font from normal to bold or back.
!
! *** Calling variables :
!     
!     font          C      host   i/p   : Font type              [bold or norm]
!     lun           I      host   i/p   : Output unit number                [-]
!
! *** Author:
!
!     R.Akers, D3/1.36, Culham Centre for Fusion Energy, x6323
!
! *** Creation date
!
!     LOCUST-CPU        : 01/02/00
!     Converted to PGF  : 16/09/10
!
!     Modification log  :
!
!------------------------------------------------------------------------------

      implicit none

      character (len=5)             :: a_format = "(aNN)"
      character (len=4)             :: bold     = achar(27)//'[7m'
      character (len=4)             :: normal   = achar(27)//'[m '
      character (len=4), intent(in) :: font
      integer, optional, intent(in) :: lun
      integer                       :: UNIT

      if(present(lun))then
         UNIT = lun
      else
         UNIT = 6
      endif

      select case (font)
             case ( 'norm' ) 

                write(a_format(3:4), fmt="(i2.2)") len_trim(normal)
                write(UNIT,advance='NO',fmt=a_format) normal

             case ( 'bold' )
 
                write(a_format(3:4), fmt="(i2.2)") len_trim(bold  )
                write(UNIT,advance='NO',fmt=a_format) bold

             case default

                write(a_format(3:4), fmt="(i2.2)") len_trim(normal)
                write(UNIT,advance='NO',fmt=a_format) normal

      end select

end subroutine display

end module disp_mod

