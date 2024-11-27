
      module law129_upd_mod
      contains
        subroutine law129_upd(mat_id   ,titr     ,matparam ,table    ,ntable  )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
          use constant_mod
          use matparam_def_mod
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
          implicit none
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
          integer, intent(in) :: mat_id
          character(len=nchartitle), intent(in) :: titr
          type(matparam_struct_), intent(inout), target :: matparam
          integer, intent(in) :: ntable
          type(ttable), dimension(ntable) ,intent(in) ::  table
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
          integer :: numtabl,i
          integer, dimension(:), pointer :: itable
! 
          ! Get the number of material tables
          numtabl = matparam%ntable
          itable => matparam%table(1:numtabl)%notable
          write(*,*) 'numtabl = ',numtabl
          write(*,*) 'itable = ',itable(1:numtabl)
!
          ! Conversion of user material id to system material id
          call mattab_usr2sys(titr    ,mat_id  ,ntable  ,table   ,numtabl ,itable  )
          do i = 1, numtabl
            matparam%table(i)%ndim = 0
            allocate(matparam%table(i)%x(0))
            allocate(matparam%table(i)%y1d(0))
            allocate(matparam%table(i)%y2d(0,0))
            allocate(matparam%table(i)%y3d(0,0,0))
            allocate(matparam%table(i)%y4d(0,0,0,0))
          enddo
          write(*,*) 'itable = ',itable(1:numtabl)
!
        end subroutine law129_upd
      end module law129_upd_mod

