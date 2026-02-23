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
      module kinematic_hardening_chaboche_mod
      contains
      subroutine kinematic_hardening_chaboche(                                 &
        matparam ,nel      ,l_sigb   ,dsigb_dlam,sigb      ,chard    ,         &
        normxx   ,normyy   ,normzz   ,normxy    ,normyz    ,normzx   ,         &
        dpla_dlam)
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
        type(matparam_struct_),        intent(in)    :: matparam     !< Material parameters data
        integer,                       intent(in)    :: nel          !< Number of elements in the group
        integer,                       intent(in)    :: l_sigb       !< Number of backstress components
        real(kind=WP), dimension(nel,l_sigb),intent(inout) :: dsigb_dlam !< Backstress components derivative w.r.t plastic multiplier
        real(kind=WP),                 intent(in)    :: chard        !< Mixed hardening parameter
        real(kind=WP), dimension(nel,l_sigb),intent(in) :: sigb      !< Backstress components for kinematic hardening
        real(kind=WP), dimension(nel), intent(in)    :: normxx       !< 1st derivative of equivalent stress wrt stress xx
        real(kind=WP), dimension(nel), intent(in)    :: normyy       !< 1st derivative of equivalent stress wrt stress yy
        real(kind=WP), dimension(nel), intent(in)    :: normzz       !< 1st derivative of equivalent stress wrt stress zz
        real(kind=WP), dimension(nel), intent(in)    :: normxy       !< 1st derivative of equivalent stress wrt stress xy
        real(kind=WP), dimension(nel), intent(in)    :: normyz       !< 1st derivative of equivalent stress wrt stress yz
        real(kind=WP), dimension(nel), intent(in)    :: normzx       !< 1st derivative of equivalent stress wrt stress zx
        real(kind=WP),                 intent(in)    :: dpla_dlam    !< Derivative of equivalent plastic strain w.r.t plastic multiplier
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,offset
        real(kind=WP) :: akh(4),ckh(4)
!===============================================================================
!
        !=======================================================================
        !< - Chaboche-Rousselier kinematic hardening model
        !=======================================================================
        offset = matparam%iparam(20) + 1
        !< Recover kinematic hardening parameters
        ckh(1) = matparam%uparam(offset + 1) !< 
        akh(1) = matparam%uparam(offset + 2) !<
        ckh(2) = matparam%uparam(offset + 3) !< 
        akh(2) = matparam%uparam(offset + 4) !<
        ckh(3) = matparam%uparam(offset + 5) !< 
        akh(3) = matparam%uparam(offset + 6) !<
        ckh(4) = matparam%uparam(offset + 7) !< 
        akh(4) = matparam%uparam(offset + 8) !<
        !< Compute the backstress components derivative w.r.t plastic multiplier
        do i = 1,nel
          !< 1st Chaboche-Rousselier back stress
          dsigb_dlam(i, 1) = chard*(akh(1)*ckh(1)*normxx(i) - ckh(1)*sigb(i,1)*dpla_dlam)
          dsigb_dlam(i, 2) = chard*(akh(1)*ckh(1)*normyy(i) - ckh(1)*sigb(i,2)*dpla_dlam)
          dsigb_dlam(i, 3) = chard*(akh(1)*ckh(1)*normzz(i) - ckh(1)*sigb(i,3)*dpla_dlam)
          dsigb_dlam(i, 4) = chard*(akh(1)*ckh(1)*normxy(i) - ckh(1)*sigb(i,4)*dpla_dlam)
          dsigb_dlam(i, 5) = chard*(akh(1)*ckh(1)*normyz(i) - ckh(1)*sigb(i,5)*dpla_dlam)
          dsigb_dlam(i, 6) = chard*(akh(1)*ckh(1)*normzx(i) - ckh(1)*sigb(i,6)*dpla_dlam)
          !< 2nd Chaboche-Rousselier back stress
          dsigb_dlam(i, 7) = chard*(akh(2)*ckh(2)*normxx(i) - ckh(2)*sigb(i, 7)*dpla_dlam)
          dsigb_dlam(i, 8) = chard*(akh(2)*ckh(2)*normyy(i) - ckh(2)*sigb(i, 8)*dpla_dlam)
          dsigb_dlam(i, 9) = chard*(akh(2)*ckh(2)*normzz(i) - ckh(2)*sigb(i, 9)*dpla_dlam)
          dsigb_dlam(i,10) = chard*(akh(2)*ckh(2)*normxy(i) - ckh(2)*sigb(i,10)*dpla_dlam)
          dsigb_dlam(i,11) = chard*(akh(2)*ckh(2)*normyz(i) - ckh(2)*sigb(i,11)*dpla_dlam)
          dsigb_dlam(i,12) = chard*(akh(2)*ckh(2)*normzx(i) - ckh(2)*sigb(i,12)*dpla_dlam)
          !< 3rd Chaboche-Rousselier back stress
          dsigb_dlam(i,13) = chard*(akh(3)*ckh(3)*normxx(i) - ckh(3)*sigb(i,13)*dpla_dlam)
          dsigb_dlam(i,14) = chard*(akh(3)*ckh(3)*normyy(i) - ckh(3)*sigb(i,14)*dpla_dlam)
          dsigb_dlam(i,15) = chard*(akh(3)*ckh(3)*normzz(i) - ckh(3)*sigb(i,15)*dpla_dlam)
          dsigb_dlam(i,16) = chard*(akh(3)*ckh(3)*normxy(i) - ckh(3)*sigb(i,16)*dpla_dlam)
          dsigb_dlam(i,17) = chard*(akh(3)*ckh(3)*normyz(i) - ckh(3)*sigb(i,17)*dpla_dlam)
          dsigb_dlam(i,18) = chard*(akh(3)*ckh(3)*normzx(i) - ckh(3)*sigb(i,18)*dpla_dlam)  
          !< 4th Chaboche-Rousselier back stress
          dsigb_dlam(i,19) = chard*(akh(4)*ckh(4)*normxx(i) - ckh(4)*sigb(i,19)*dpla_dlam)
          dsigb_dlam(i,20) = chard*(akh(4)*ckh(4)*normyy(i) - ckh(4)*sigb(i,20)*dpla_dlam)
          dsigb_dlam(i,21) = chard*(akh(4)*ckh(4)*normzz(i) - ckh(4)*sigb(i,21)*dpla_dlam)
          dsigb_dlam(i,22) = chard*(akh(4)*ckh(4)*normxy(i) - ckh(4)*sigb(i,22)*dpla_dlam)
          dsigb_dlam(i,23) = chard*(akh(4)*ckh(4)*normyz(i) - ckh(4)*sigb(i,23)*dpla_dlam)
          dsigb_dlam(i,24) = chard*(akh(4)*ckh(4)*normzx(i) - ckh(4)*sigb(i,24)*dpla_dlam)
        enddo
!
      end subroutine kinematic_hardening_chaboche
      end module kinematic_hardening_chaboche_mod
