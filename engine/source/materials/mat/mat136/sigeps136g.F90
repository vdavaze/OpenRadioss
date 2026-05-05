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
          nel    ,matparam,rho    ,dmg     ,thk0   ,                           &
          depsxx ,depsyy ,depsxy  ,depsyz  ,depszx ,                           &                             
          depbxx ,depbyy ,depbxy  ,                                            &
          sigoxx ,sigoyy ,sigoxy  ,sigoyz  ,sigozx ,                           &
          signxx ,signyy ,signxy  ,signyz  ,signzx ,                           &
          momnxx ,momnyy ,momnxy  ,                                            &
          ssp    ,et     ,gs      ,nuvar   ,uvar   ,                           &
          shf    )
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
          real(kind=WP), dimension(nel,3),intent(inout):: dmg      !< Bending damage at current time
          real(kind=WP), dimension(nel), intent(in)    :: thk0     !< Initial thickness
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
          real(kind=WP), dimension(nel), intent(inout) :: gs       !< Correction factor for transverse shear
          integer,                       intent(in)    :: nuvar    !< Number of user variables
          real(kind=WP), dimension(nel,nuvar), intent(inout) :: uvar !< User variables at current time step
          real(kind=WP), dimension(nel), intent(in)    :: shf      !< Shear correction factor force coefficient
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
          integer :: i
          real(kind=WP) :: young,nu,shear,lambda_m, mu_m, gamma
          real(kind=wp), dimension(nel) :: lambda_b, mu_b, k0_1, k0_2
          real(kind=WP) :: dmax, f_t, f_c, center, radius, kappa1, kappa2
          real(kind=WP) :: tr_kb, phi_b, Y1, Y2, d1_trial, d2_trial
          real(kind=WP), dimension(nel) :: xi_tr, epbxx, epbyy, epbxy, xi_kap1,&
            xi_kap2
          real(kind=WP) :: dkap1dkxx, dkap1dkyy, dkap1dkxy, dkap2dkxx,         &
            dkap2dkyy, dkap2dkxy, phi_b_pos, phi_b_neg
!
          !=====================================================================
          !< - Initialisation of computation on time step
          !=====================================================================
!
          !< Recovering real model parameters
          young    = matparam%young
          nu       = matparam%nu
          shear    = matparam%shear
          lambda_m = matparam%uparam(1)
          mu_m     = matparam%uparam(2)
          lambda_b(1:nel) = thk0(1:nel)*matparam%uparam(3)
          mu_b(1:nel) = thk0(1:nel)*matparam%uparam(4)
          gs(1:nel) = shear*shf(1:nel)
          f_t      = matparam%uparam(5)
          f_c      = matparam%uparam(6)
          k0_1(1:nel) = f_t*f_t*(one - nu*nu) / (six*young*thk0(1:nel)*thk0(1:nel))
          k0_2(1:nel) = f_c*f_c*(one - nu*nu) / (six*young*thk0(1:nel)*thk0(1:nel))
          gamma    = matparam%uparam(7)
          dmax     = matparam%uparam(8)
!
          !< Recover the damage variable at the beginning of the time step
          do i = 1, nel
            uvar(i,1) = uvar(i,1) + depbxx(i)
            uvar(i,2) = uvar(i,2) + depbyy(i)
            uvar(i,3) = uvar(i,3) + depbxy(i)
            epbxx(i)  = uvar(i,1)
            epbyy(i)  = uvar(i,2)
            epbxy(i)  = uvar(i,3)
          enddo
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
            signyz(i) = sigoyz(i) +  gs(i)*depsyz(i)
            signzx(i) = sigozx(i) +  gs(i)*depszx(i)
            ! -> Hourglass stabilization variable
            et(i) = one
            ! -> Sound speed
            ssp(i) = sqrt((young/(one - nu*nu))/rho(i))
          enddo
!
          ! !=====================================================================
          ! !< - Computation of the plasticity
          ! !=====================================================================
          ! ! !< Loop over elements
          ! ! do i = 1, nel
          !  Mf_pos = (thk0(i)/four)*f_c*(one - (signxx(i)/(f_c*thk0(i)))**2) - &
          !           Om_sup*sigy*(thk0(i)/two - dsup) -                        &
          !           Om_inf*sigy*(thk0(i)/two - dinf)
          ! ! enddo
!
          !=====================================================================
          !< - Computation of the bending damage and softening factor
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            !< Trace of the bending curvature increment tensor
            tr_kb = epbxx(i) + epbyy(i)
            !< Eigenvalues of the bending curvature tensor
            center = (epbxx(i) + epbyy(i))/two
            radius = sqrt(((epbxx(i) - epbyy(i))*half)**2 + (epbxy(i)*half)**2)
            kappa1 = center + radius
            kappa2 = center - radius
            !< Free elastic bending energy (without damage)
            phi_b = half*lambda_b(i) * tr_kb*tr_kb                             &
                           + mu_b(i) * ( kappa1*kappa1  + kappa2*kappa2 )
            !< Thermodynamic force Y1 (positive bending - tension inner face)
            Y1 = zero
            if (tr_kb  >= zero) Y1 = Y1 + half*lambda_b(i) * tr_kb*tr_kb
            if (kappa1 >= zero) Y1 = Y1 + mu_b(i) * kappa1*kappa1
            if (kappa2 >= zero) Y1 = Y1 + mu_b(i) * kappa2*kappa2
            phi_b_pos = Y1
            Y1 = ((one - gamma)/(one + dmg(i,2))**2) * Y1
            !< Thermodynamic force Y2 (negative bending - tension outer face)
            Y2 = zero
            if (tr_kb  < zero) Y2 = Y2 + half*lambda_b(i) * tr_kb*tr_kb
            if (kappa1 < zero) Y2 = Y2 + mu_b(i) * kappa1*kappa1
            if (kappa2 < zero) Y2 = Y2 + mu_b(i) * kappa2*kappa2
            phi_b_neg = Y2
            Y2 = ((one - gamma)/(one + dmg(i,3))**2) * Y2
            !< Evolution law of the bending damage
            ! --> Positive bending, inner face
            if (Y1 >= k0_1(i)) then
              d1_trial = sqrt(phi_b_pos / k0_1(i)) - one
              d1_trial = max(zero, min(d1_trial, dmax))
              dmg(i,2) = max(dmg(i,2), d1_trial)
            end if
            ! --> Negative bending, outer face
            if (Y2 >= k0_2(i)) then
              d2_trial = sqrt(phi_b_neg / k0_2(i)) - one
              d2_trial = max(zero, min(d2_trial, dmax))
              dmg(i,3) = max(dmg(i,3), d2_trial)
            end if
            !< Softening factor for the moment Mxx and Myy
            if (tr_kb >= zero) then
              xi_tr(i) = (one + gamma*dmg(i,2)) / (one + dmg(i,2))
            else
              xi_tr(i) = (one + gamma*dmg(i,3)) / (one + dmg(i,3))
            end if
            !< Softening factor for kappa1 and kappa2
            if (kappa1 >= zero) then
              xi_kap1(i) = (one + gamma*dmg(i,2)) / (one + dmg(i,2))
            else
              xi_kap1(i) = (one + gamma*dmg(i,3)) / (one + dmg(i,3))
            endif
            if (kappa2 >= zero) then
              xi_kap2(i) = (one + gamma*dmg(i,2)) / (one + dmg(i,2))
            else
              xi_kap2(i) = (one + gamma*dmg(i,3)) / (one + dmg(i,3))
            endif
            dmg(i,1) = max(dmg(i,2), dmg(i,3)) !< Maximum of the positive and negative bending damage
            !< Update the bending moments with the softening factors
            dkap1dkxx = half + (epbxx(i) - epbyy(i))/(four*max(radius,em20))
            dkap1dkyy = half - (epbxx(i) - epbyy(i))/(four*max(radius,em20))
            dkap1dkxy = (epbxy(i)*half)/(two*max(radius,em20))
            dkap2dkxx = half - (epbxx(i) - epbyy(i))/(four*max(radius,em20))
            dkap2dkyy = half + (epbxx(i) - epbyy(i))/(four*max(radius,em20))
            dkap2dkxy = -(epbxy(i)*half)/(two*max(radius,em20))
            momnxx(i) = lambda_b(i)*tr_kb*xi_tr(i) + two*mu_b(i)*(             &
                            xi_kap1(i)*kappa1*dkap1dkxx +                      &
                            xi_kap2(i)*kappa2*dkap2dkxx ) 
            momnyy(i) = lambda_b(i)*tr_kb*xi_tr(i) + two*mu_b(i)*(             & 
                            xi_kap1(i)*kappa1*dkap1dkyy +                      &
                            xi_kap2(i)*kappa2*dkap2dkyy )
            momnxy(i) = two*mu_b(i)*(                                          &
                            xi_kap1(i)*kappa1*dkap1dkxy +                      &
                            xi_kap2(i)*kappa2*dkap2dkxy )
          enddo       
!
          !=====================================================================
        end subroutine sigeps136g
      end module sigeps136g_mod
