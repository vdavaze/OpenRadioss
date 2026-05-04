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
      module sigeps136g_mod
      contains
        subroutine sigeps136g(                                                 &
          nel    ,matparam,rho    ,dmg     ,                                   &
          depsxx ,depsyy ,depsxy  ,depsyz  ,depszx ,                           &                             
          depbxx ,depbyy ,depbxy  ,                                            &
          sigoxx ,sigoyy ,sigoxy  ,sigoyz  ,sigozx ,                           &
          momoxx ,momoyy ,momoxy  ,                                            &
          signxx ,signyy ,signxy  ,signyz  ,signzx ,                           &
          momnxx ,momnyy ,momnxy  ,                                            &
          ssp    ,et     ,gs      )

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
          integer,                       intent(in)    :: nel      !< Number of elements in the group
          type(matparam_struct_),        intent(in)    :: matparam !< Material parameters data
          real(kind=WP), dimension(nel), intent(in)    :: rho      !< Density at current time
          real(kind=WP), dimension(nel,2),intent(inout):: dmg      !< Bending damage at current time
          real(kind=WP), dimension(nel), intent(in)    :: depsxx   !< Membrane strain increment xx
          real(kind=WP), dimension(nel), intent(in)    :: depsyy   !< Membrane strain increment yy
          real(kind=WP), dimension(nel), intent(in)    :: depsxy   !< Membrane strain increment xy
          real(kind=WP), dimension(nel), intent(in)    :: depsyz   !< Membrane strain increment yz
          real(kind=WP), dimension(nel), intent(in)    :: depszx   !< Membrane strain increment zx
          real(kind=WP), dimension(nel), intent(in)    :: depbxx   !< Bending curvature increment xx
          real(kind=WP), dimension(nel), intent(in)    :: depbyy   !< Bending curvature increment yy
          real(kind=WP), dimension(nel), intent(in)    :: depbxy   !< Bending curvature increment xy
          real(kind=WP), dimension(nel), intent(in)    :: sigoxx   !< Membrane force xx at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: sigoyy   !< Membrane force yy at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: sigoxy   !< Membrane force xy at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: sigoyz   !< Membrane force yz at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: sigozx   !< Membrane force zx at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: momoxx   !< Bending moment xx at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: momoyy   !< Bending moment yy at previous time step
          real(kind=WP), dimension(nel), intent(in)    :: momoxy   !< Bending moment xy at previous time step
          real(kind=WP), dimension(nel), intent(inout) :: signxx   !< Membrane force xx at current time step
          real(kind=WP), dimension(nel), intent(inout) :: signyy   !< Membrane force yy at current time step
          real(kind=WP), dimension(nel), intent(inout) :: signxy   !< Membrane force xy at current time step
          real(kind=WP), dimension(nel), intent(inout) :: signyz   !< Membrane force yz at current time step
          real(kind=WP), dimension(nel), intent(inout) :: signzx   !< Membrane force zx at current time step
          real(kind=WP), dimension(nel), intent(inout) :: momnxx   !< Bending moment xx at current time step
          real(kind=WP), dimension(nel), intent(inout) :: momnyy   !< Bending moment yy at current time step
          real(kind=WP), dimension(nel), intent(inout) :: momnxy   !< Bending moment xy at current time
          real(kind=WP), dimension(nel), intent(inout) :: ssp      !< Current sound speed
          real(kind=WP), dimension(nel), intent(inout) :: et       !< Hourglass stabilization variable
          real(kind=WP), dimension(nel), intent(in)    :: gs       !< Correction factor for transverse shear
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
          integer :: i
          real(kind=WP) :: young,nu,lambda_m, mu_m, lambda_b, mu_b, mu_tsh
          real(kind=WP) :: k0_1, k0_2, dmax1, dmax2
          real(kind=WP) :: tr_kb, phi_b, Y1, Y2, d1_trial, d2_trial
          real(kind=WP), dimension(nel) :: xi_tr, xi_xy
!
          !=====================================================================
          !< - Initialisation of computation on time step
          !=====================================================================
!
          !< Recovering real model parameters
          young    = matparam%young
          nu       = matparam%nu
          lambda_m = matparam%uparam(2)
          mu_m     = matparam%uparam(3)
          lambda_b = matparam%uparam(4)
          mu_b     = matparam%uparam(5)
          mu_tsh   = matparam%uparam(6)
          k0_1     = matparam%uparam(7)
          k0_2     = matparam%uparam(8)
          dmax1    = matparam%uparam(9)
          dmax2    = matparam%uparam(10)
!
          !=====================================================================
          !< - Computation of the trial stress tensor (elastic predictor)
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            ! -> Membrane forces
            signxx(i) = sigoxx(i) + lambda_m*(depsxx(i) + depsyy(i)) +         &
                                                  two*mu_m*depsxx(i)
            signyy(i) = sigoyy(i) + lambda_m*(depsxx(i) + depsyy(i)) +         &
                                                  two*mu_m*depsyy(i)
            signxy(i) = sigoxy(i) +   mu_m*depsxy(i)
            signyz(i) = sigoyz(i) + mu_tsh*depsyz(i)*gs(i)
            signzx(i) = sigozx(i) + mu_tsh*depszx(i)*gs(i)
            ! -> Bending moments
            momnxx(i) = momoxx(i) + lambda_b*(depbxx(i) + depbyy(i)) +         &
                                                  two*mu_b*depbxx(i)
            momnyy(i) = momoyy(i) + lambda_b*(depbxx(i) + depbyy(i)) +         &
                                                  two*mu_b*depbyy(i)
            momnxy(i) = momoxy(i) + mu_b*depbxy(i)
            ! -> Hourglass stabilization variable
            et(i) = one
            ! -> Sound speed
            ssp(i) = sqrt((young/(one - nu*nu))/rho(i))
          enddo
!
          !=====================================================================
          !< - Computation of the bending damage and softening factor
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            !< Trace of the bending curvature increment tensor
            tr_kb = depbxx(i) + depbyy(i)
            !< Bending damage energy release rate
            phi_b = half*lambda_b * tr_kb*tr_kb                                &
              + mu_b * (depbxx(i)*depbxx(i) + depbyy(i)*depbyy(i)              &
              + half*half*depbxy(i)*depbxy(i))
            !< Thermodynamic force conjugate to the bending damages
            Y1 = phi_b / ((one + dmg(i,1)) * (one + dmg(i,1)))
            Y2 = phi_b / ((one + dmg(i,2)) * (one + dmg(i,2)))
            !< Evolution law of the bending damage
            ! --> Positive bending, inner face
            if (tr_kb > zero) then 
              if (Y1 > k0_1) then
                d1_trial = sqrt(phi_b / k0_1) - one
                d1_trial = max(zero, min(d1_trial, dmax1))
                dmg(i,1) = max(dmg(i,1), d1_trial)
              end if
            ! --> Negative bending, outer face
            else 
              if (Y2 > k0_2) then
                d2_trial = sqrt(phi_b / k0_2) - one
                d2_trial = max(zero, min(d2_trial, dmax2))
                dmg(i,2) = max(dmg(i,2), d2_trial)
              end if
            endif
            !< Softening factor for the moment Mxx and Myy
            if (tr_kb >= zero) then
              xi_tr(i) = one / (one + dmg(i,1))
            else
              xi_tr(i) = one / (one + dmg(i,2))
            end if
            !< Softening factor for the moment Mxy
            if (depbxy(i) >= zero) then
              xi_xy(i) = one / (one + dmg(i,1))
            else
              xi_xy(i) = one / (one + dmg(i,2))
            endif

          enddo
!
          !=========================================================================
        end subroutine sigeps136g
      end module sigeps136g_mod
