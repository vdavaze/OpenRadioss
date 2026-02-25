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
      module yield_criterion_vonmises_mod
      contains
      subroutine yield_criterion_vonmises(                                     &
              nel      ,nindx    ,indx     ,seq      ,eltype   ,               &
              signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,     &
              normxx   ,normyy   ,normzz   ,normxy   ,normyz   ,normzx   )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use constant_mod
        use mvsiz_mod
        use precision_mod, only : WP
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        integer,                       intent(in)    :: nel      !< Number of elements in the group
        integer,                       intent(in)    :: nindx    !< Number of elements to consider in the computation (for partial updates)
        integer,       dimension(nel), intent(in)    :: indx     !< Indices of the elements to consider in the computation (for partial updates)
        real(kind=WP), dimension(nel), intent(inout) :: seq      !< Equivalent stress
        integer,                       intent(in)    :: eltype   !< Element type
        real(kind=WP), dimension(nel), intent(in)    :: signxx   !< Current stress xx
        real(kind=WP), dimension(nel), intent(in)    :: signyy   !< Current stress yy
        real(kind=WP), dimension(nel), intent(in)    :: signzz   !< Current stress zz
        real(kind=WP), dimension(nel), intent(in)    :: signxy   !< Current stress xy
        real(kind=WP), dimension(nel), intent(in)    :: signyz   !< Current stress yz
        real(kind=WP), dimension(nel), intent(in)    :: signzx   !< Current stress zx
        real(kind=WP), dimension(nel), intent(inout) :: normxx   !< 1st derivative of equivalent stress wrt stress xx
        real(kind=WP), dimension(nel), intent(inout) :: normyy   !< 1st derivative of equivalent stress wrt stress yy
        real(kind=WP), dimension(nel), intent(inout) :: normzz   !< 1st derivative of equivalent stress wrt stress zz
        real(kind=WP), dimension(nel), intent(inout) :: normxy   !< 1st derivative of equivalent stress wrt stress xy
        real(kind=WP), dimension(nel), intent(inout) :: normyz   !< 1st derivative of equivalent stress wrt stress yz
        real(kind=WP), dimension(nel), intent(inout) :: normzx   !< 1st derivative of equivalent stress wrt stress zx
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,ii
!===============================================================================
!
        !=======================================================================
        !< - Von Mises yield criterion and its derivatives
        !=======================================================================
        !< Solid element
        if (eltype == 1) then 
#include "vectorize.inc"
          do ii = 1,nindx
            i = indx(ii)
            !< Equivalent stress
            seq(i) =  half*(signyy(i)-signzz(i))**2 +                          &                 
                      half*(signzz(i)-signxx(i))**2 +                          & 
                      half*(signxx(i)-signyy(i))**2 +                          &
                     three*(signxy(i))**2           +                          &
                     three*(signyz(i))**2           +                          &
                     three*(signzx(i))**2 
            seq(i) = sqrt(seq(i))
            !< First order derivative of eq. stress
            normxx(i) = half*(one/max(seq(i),em20))*                           &
                               (two*signxx(i) - signyy(i) - signzz(i))
            normyy(i) = half*(one/max(seq(i),em20))*                           &
                               (two*signyy(i) - signxx(i) - signzz(i))
            normzz(i) = half*(one/max(seq(i),em20))*                           &
                               (two*signzz(i) - signyy(i) - signxx(i))
            normxy(i) = (one/max(seq(i),em20))*three*signxy(i)
            normyz(i) = (one/max(seq(i),em20))*three*signyz(i)
            normzx(i) = (one/max(seq(i),em20))*three*signzx(i)              
          enddo 
        !< Shell element
        elseif (eltype == 2) then 
#include "vectorize.inc"
          do ii = 1,nindx
            i = indx(ii)
            !< Equivalent stress
            seq(i) = signxx(i)**2 + signyy(i)**2 -  signxx(i)*signyy(i)        &
                                                 + three*(signxy(i)**2)
            seq(i) = sqrt(seq(i))
            !< First order derivative of eq. stress
            normxx(i) =  half*(one/max(seq(i),em20))*(two*signxx(i) - signyy(i))
            normyy(i) =  half*(one/max(seq(i),em20))*(two*signyy(i) - signxx(i))
            normzz(i) = - normxx(i) - normyy(i)
            normxy(i) = three*(one/max(seq(i),em20))*signxy(i)
            normyz(i) = zero
            normzx(i) = zero
          enddo
        endif
!
      end subroutine yield_criterion_vonmises
      end module yield_criterion_vonmises_mod
