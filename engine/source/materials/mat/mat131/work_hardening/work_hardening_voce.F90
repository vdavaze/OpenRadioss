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
      module work_hardening_voce_mod
      contains
      subroutine work_hardening_voce(                                          &
        matparam ,nel      ,sigy     ,pla      ,dsigy_dpla)
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
        type(matparam_struct_),        intent(in)    :: matparam   !< Material parameters data
        integer,                       intent(in)    :: nel        !< Number of elements in the group
        real(kind=WP), dimension(nel), intent(inout) :: sigy       !< Equivalent stress
        real(kind=WP), dimension(nel), intent(inout) :: pla        !< Cumulated plastic strain
        real(kind=WP), dimension(nel), intent(inout) :: dsigy_dpla !< Derivative of eq. stress w.r.t. cumulated plastic strain
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: offset,i
        real(kind=WP) :: r0,q1,b1,q2,b2,q3,b3
!===============================================================================
!
        !=======================================================================
        !< - Voce work hardening model
        !=======================================================================
        offset = matparam%iparam(4)
        !< Recover work hardening parameters
        r0 = matparam%uparam(offset + 1) !< Initial yield stress
        q1 = matparam%uparam(offset + 2) !< Voce 1 saturation stress
        b1 = matparam%uparam(offset + 3) !< Voce 1 saturation rate
        q2 = matparam%uparam(offset + 4) !< Voce 2 saturation stress
        b2 = matparam%uparam(offset + 5) !< Voce 2 saturation rate
        q3 = matparam%uparam(offset + 6) !< Voce 3 saturation stress
        b3 = matparam%uparam(offset + 7) !< Voce 3 saturation rate
        do i = 1,nel
          sigy(i) = r0 + q1*(one - exp(-b1*pla(i)))                            &
                       + q2*(one - exp(-b2*pla(i)))                            &
                       + q3*(one - exp(-b3*pla(i)))
          dsigy_dpla(i) = q1*b1*exp(-b1*pla(i)) +                              &
                          q2*b2*exp(-b2*pla(i)) +                              &
                          q3*b3*exp(-b3*pla(i))
        enddo
!
      end subroutine work_hardening_voce
      end module work_hardening_voce_mod
