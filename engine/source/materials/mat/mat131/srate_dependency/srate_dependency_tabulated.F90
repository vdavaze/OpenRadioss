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
      module srate_dependency_tabulated_mod
      contains
      subroutine srate_dependency_tabulated(                                   &
        matparam ,nel      ,sigy     ,epsd     ,dsigy_dpla,vpflag   ,timestep, &
        nvartmp  ,vartmp   )
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
        real(kind=WP),   dimension(nel), intent(inout) :: epsd       !< Strain rate
        real(kind=WP),   dimension(nel), intent(inout) :: dsigy_dpla !< Derivative of eq. stress w.r.t. cumulated plastic strain
        integer,                         intent(in)    :: vpflag     !< Strain rate dependency flag
        real(kind=WP),                   intent(in)    :: timestep   !< Time step
        integer,                         intent(in)    :: nvartmp    !< Number of variables used in tabulated strain rate dependency
        integer, dimension(nel,nvartmp), intent(inout) :: vartmp     !< Temporary variables for tabulated strain rate dependency
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,offset,offset_var,ipos(nel,1)
        real(kind=WP) :: xvec(nel,1),dfact_depsd(nel),srate_fac(nel)
!===============================================================================
!
        !=======================================================================
        !< - Tabulated strain rate dependency model
        !=======================================================================
        !< Table offset
        offset = matparam%iparam(6)
        offset_var = matparam%iparam(8)
        !< Prepare input vectors for interpolation
        do i = 1,nel
          xvec(i,1) = epsd(i)
          ipos(i,1) = vartmp(i,offset_var + 1)
        enddo
        !< Interpolate to get srate_fac and dfact_depsd
        call table_mat_vinterp(matparam%table(offset+1),nel,nel,ipos,xvec,     &
          srate_fac,dfact_depsd)
        !< Update temporary variables
        do i = 1,nel
          sigy(i) = sigy(i)*srate_fac(i)
          dsigy_dpla(i) = dsigy_dpla(i)*srate_fac(i)
          vartmp(i,offset_var + 1) = ipos(i,1)
        enddo
        !< For full viscoplasticity, add strain rate contribution to dsigy_dpla
        if (vpflag == 4) then
          do i = 1,nel
            dsigy_dpla(i) = dsigy_dpla(i) + sigy(i)*dfact_depsd(i)*(one/timestep)
          enddo
        endif
!
      end subroutine srate_dependency_tabulated
      end module srate_dependency_tabulated_mod
