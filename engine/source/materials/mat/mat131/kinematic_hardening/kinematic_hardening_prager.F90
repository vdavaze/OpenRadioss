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
      module kinematic_hardening_prager_mod
      contains
      subroutine kinematic_hardening_prager(                                   &
        nel      ,nindx    ,indx     ,l_sigb    ,dsigb_dlam,dsigy_dpla,        &
        normxx   ,normyy   ,normzz   ,normxy    ,normyz    ,normzx    ,        &
        chard    )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use constant_mod
        use precision_mod, only : WP
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        integer,                       intent(in)    :: nel          !< Number of elements in the group
        integer,                       intent(in)    :: nindx        !< Number of elements to consider in the computation (for partial updates)
        integer,       dimension(nel), intent(in)    :: indx         !< Indices of the elements to consider in the computation (for partial updates)
        integer,                       intent(in)    :: l_sigb       !< Number of backstress components
        real(kind=WP), dimension(nel,l_sigb),intent(inout) :: dsigb_dlam !< Backstress components derivative w.r.t plastic multiplier
        real(kind=WP), dimension(nel), intent(in)    :: dsigy_dpla   !< Derivative of yield stress wrt equivalent plastic strain
        real(kind=WP), dimension(nel), intent(in)    :: normxx       !< 1st derivative of equivalent stress wrt stress xx
        real(kind=WP), dimension(nel), intent(in)    :: normyy       !< 1st derivative of equivalent stress wrt stress yy
        real(kind=WP), dimension(nel), intent(in)    :: normzz       !< 1st derivative of equivalent stress wrt stress zz
        real(kind=WP), dimension(nel), intent(in)    :: normxy       !< 1st derivative of equivalent stress wrt stress xy
        real(kind=WP), dimension(nel), intent(in)    :: normyz       !< 1st derivative of equivalent stress wrt stress yz
        real(kind=WP), dimension(nel), intent(in)    :: normzx       !< 1st derivative of equivalent stress wrt stress zx
        real(kind=WP),                 intent(in)    :: chard        !< Mixed hardening parameter
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,ii
!===============================================================================
!
        !=======================================================================
        !< - Prager kinematic hardening model
        !=======================================================================
        !< Compute the backstress components derivative w.r.t plastic multiplier
#include "vectorize.inc"
        do ii = 1,nindx
          i = indx(ii)
          dsigb_dlam(i,1) = two_third*chard*dsigy_dpla(i)*normxx(i)
          dsigb_dlam(i,2) = two_third*chard*dsigy_dpla(i)*normyy(i)
          dsigb_dlam(i,3) = two_third*chard*dsigy_dpla(i)*normzz(i)
          dsigb_dlam(i,4) = two_third*chard*dsigy_dpla(i)*normxy(i)
          dsigb_dlam(i,5) = two_third*chard*dsigy_dpla(i)*normyz(i)
          dsigb_dlam(i,6) = two_third*chard*dsigy_dpla(i)*normzx(i)
        enddo
!
      end subroutine kinematic_hardening_prager
      end module kinematic_hardening_prager_mod
