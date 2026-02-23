!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2026 Altair Engineering Inc.
!Copyright>
!Copyright>        This program is free software: you can redistribute it and/or modify
!Copyright>        it under the terms of the GNU Affero General Public License as published by
!Copyright>        the Free Software Foundation, either version 3 of the License, or
!Copyright>        (at your option) any later version.
!Copyright>
!Copyright>        This program is distributed in the hope that it will be useful,
!Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
!Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!Copyright>        GNU Affero General Public License for more details.
!Copyright>
!Copyright>        You should have received a copy of the GNU Affero General Public License
!Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
!Copyright>
!Copyright>
!Copyright>        Commercial Alternative: Altair Radioss Software
!Copyright>
!Copyright>        As an alternative to this open-source version, Altair also offers Altair Radioss
!Copyright>        software under a commercial license.  Contact Altair to discuss further if the
!Copyright>        commercial version may interest you: https://www.altair.com/radioss/.
      module work_hardening_tabulated_mod
      contains
      subroutine work_hardening_tabulated(                                     &
        matparam ,nel      ,sigy     ,pla      ,epsd     ,dsigy_dpla,nvartmp  ,&
        vartmp   ,vpflag   ,timestep )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use constant_mod
        use table_mat_vinterp_mod
        use precision_mod, only : WP
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        type(matparam_struct_),          intent(in)    :: matparam   !< Material parameters data
        integer,                         intent(in)    :: nel        !< Number of elements in the group
        real(kind=WP),   dimension(nel), intent(inout) :: sigy       !< Equivalent stress
        real(kind=WP),   dimension(nel), intent(inout) :: pla        !< Cumulated plastic strain
        real(kind=WP),   dimension(nel), intent(inout) :: epsd       !< Strain rate
        real(kind=WP),   dimension(nel), intent(inout) :: dsigy_dpla !< Derivative of eq. stress w.r.t. cumulated plastic strain
        integer,                         intent(in)    :: nvartmp    !< Number of variables used in tabulated hardening
        integer, dimension(nel,nvartmp), intent(inout) :: vartmp     !< Temporary variables for tabulated hardening
        integer,                         intent(in)    :: vpflag     !< Strain rate dependency flag
        real(kind=WP),                   intent(in)    :: timestep   !< Time step
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,offset,ipos(nel,2)
        real(kind=WP) :: xvec(nel,2),dsigy_depsd(nel)
        logical :: flag_extrap
!===============================================================================
!
        !=======================================================================
        !< - Tabulated work hardening model
        !=======================================================================
        offset = matparam%iparam(4)
        !< Recover flat extrapolation flag from work hardening parameters
        flag_extrap = (matparam%uparam(offset + 1) == 0)
        !< Prepare input vectors for interpolation
        do i = 1,nel
          xvec(i,1)   = pla(i)
          xvec(i,2)   = epsd(i)
          ipos(i,1:2) = vartmp(i,1:2)
        enddo
        !< Interpolate to get sigy and dsigy_dpla
        call table_mat_vinterp(matparam%table(1),nel,nel,ipos,xvec,sigy,       &
          dsigy_dpla,flag_extrap,.false.)
        !< Update temporary variables
        do i = 1,nel
          vartmp(i,1:2) = ipos(i,1:2)
        enddo
        !< For full viscoplasticity, add strain rate contribution to dsigy_dpla
        if (vpflag == 4) then
          do i = 1,nel
            dsigy_dpla(i) = dsigy_dpla(i) + dsigy_depsd(i)*(one/timestep)
          enddo
        endif
!
      end subroutine work_hardening_tabulated
      end module work_hardening_tabulated_mod
