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
      module self_heating_taylor_mod
      contains
      subroutine self_heating_taylor(                                          &
        matparam ,nel      ,sigy     ,dtemp_dpla,dpla     ,epsd     ,vpflag   ,&
        timestep )
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
        real(kind=WP), dimension(nel), intent(inout) :: dtemp_dpla !< Derivative of temperature w.r.t. cumulated plastic strain
        real(kind=WP), dimension(nel), intent(in)    :: dpla       !< Increment of cumulated plastic strain
        real(kind=WP), dimension(nel), intent(in)    :: epsd       !< Equivalent strain rate
        integer,                       intent(in)    :: vpflag     !< Viscoplastic formulation flag
        real(kind=WP),                 intent(in)    :: timestep   !< Current time step
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: offset,i
        real(kind=WP) :: eta,cp,deis,dead,rho
        real(kind=WP), dimension(nel) :: weight,dweight
!===============================================================================
!
        !=======================================================================
        !< - Taylor-Quinney (extended) self-heating model
        !=======================================================================
        offset = matparam%iparam(16)
        !< Recover self heating parameters
        eta  = matparam%uparam(offset + 1) !< Taylor-Quinney coefficient
        cp   = matparam%uparam(offset + 2) !< Thermal massic capacity
        deis = matparam%uparam(offset + 3) !< Strain rates for the beginning of adiabatic transition
        dead = matparam%uparam(offset + 4) !< Strain rates for the end of adiabatic transition
        rho  = matparam%rho0               !< Material initial density
        !< Strain rate weight factor computation
        if (vpflag /= 4) then
          do i = 1,nel
            if (epsd(i) < deis) then
              weight(i) = zero
            elseif (epsd(i) > dead) then
              weight(i) = one
            else
              weight(i) = ((epsd(i)-deis)**2)*(three*dead - two*epsd(i) - deis)&
                                                              /((dead-deis)**3)
            endif
            dweight(i) = zero
          enddo
        else
          do i = 1,nel
            if (epsd(i) < deis) then
              weight(i)  = zero
              dweight(i) = zero
            elseif (epsd(i) > dead) then
              weight(i)  = one
              dweight(i) = zero
            else
              weight(i) = ((epsd(i)-deis)**2)*(three*dead - two*epsd(i) - deis)&
                                                              /((dead-deis)**3)
              dweight(i) = (six*(epsd(i)-deis)*(dead - epsd(i)))/              &
                                     (timestep*((dead-deis)**3))
            endif
          enddo
        endif
        !< Update derivative of temperature w.r.t. cumulated plastic strain
        do i = 1,nel
          dtemp_dpla(i) = (eta/(rho*cp))*sigy(i)*(dweight(i)*dpla(i) + weight(i))  
        enddo  
!
      end subroutine self_heating_taylor
      end module self_heating_taylor_mod
