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
          real(kind=WP) :: young,nu,shear,lambda_m, mu_m, gamma,Dm11,Dm12,     &
          real(kind=wp), dimension(nel) :: lambda_b, mu_b, k0_1, k0_2
          real(kind=WP) :: dmax, f_t, f_c, center, radius, kappa1, kappa2
          real(kind=WP) :: tr_kb, phi_b, Y1, Y2, d1_trial, d2_trial
          real(kind=WP), dimension(nel) :: xi_tr, epbxx, epbyy, epbxy, xi_kap1,&
            xi_kap2, Db11, Db12, Db21, Db22, Db33
          real(kind=WP) :: dkap1dkxx, dkap1dkyy, dkap1dkxy, dkap2dkxx,         &
            dkap2dkyy, dkap2dkxy, phi_b_pos, phi_b_neg
!
          !=====================================================================
          !< - Initialisation of computation on time step
          !=====================================================================
!
          !< Recover integer model parameters
          nlayer   = matparam%iparam(1)                      !< Number of reinforcement layers
          !< Recovering real model parameters
          young    = matparam%young                          !< Concrete Young's modulus
          nu       = matparam%nu                             !< Concrete Poisson's ratio
          shear    = matparam%shear                          !< Concrete shear modulus
          lambda_m = matparam%uparam(1)                      !< Membrane Lamé parameter (precomputed in hm_read_mat136)
          mu_m     = matparam%uparam(2)                      !< Membrane shear modulus (precomputed in hm_read_mat136)
          f_t      = matparam%uparam(5)                      !< Concrete Tensile strength
          f_c      = matparam%uparam(6)                      !< Concrete Compressive strength
          gamma    = matparam%uparam(7)                      !< Concrete 
          dmax     = matparam%uparam(8)                      !< Maximum damage
          sigy     = matparam%uparam(9)                      !< Yield stress
          omega_x  = matparam%uparam(10)                     !< Area of the steel reinforcement in x direction per unit length of concrete 
          omega_y  = matparam%uparam(11)                     !< Area of the steel reinforcement in y direction per unit length of concrete 
          rho_x    = matparam%uparam(12)                     !< Through thickness position of the steel reinforcement in x direction (normalized by the initial thickness)
          rho_y    = matparam%uparam(13)                     !< Through thickness position of the steel reinforcement in y direction (normalized by the initial thickness)
          Dm11     = lambda_m + two*mu_m                     !< Membrane stiffness D11 = lambda + 2*mu
          Dm12     = lambda_m                                !< Membrane stiffness D12 = lambda
!
          !< Computation of some real parameters
          gs(1:nel) = shear*shf(1:nel)                       !< Correction factor for transverse shear 
          lambda_b(1:nel) = thk0(1:nel)*matparam%uparam(3)   !< Bending Lamé parameter (precomputed in hm_read_mat136)
          mu_b(1:nel) = thk0(1:nel)*matparam%uparam(4)       !< Bending shear modulus (precomputed in hm_read_mat136)
          k0_1(1:nel) = f_t*f_t*(one - nu*nu) / (six*young*thk0(1:nel)*thk0(1:nel)) !< Bending damage threshold in positive bending (tension inner face)
          k0_2(1:nel) = f_c*f_c*(one - nu*nu) / (six*young*thk0(1:nel)*thk0(1:nel)) !< Bending damage threshold in negative bending (tension outer face)
!
          !< Recompute the bending curvature tensor at current time step
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
          !< - COMPUTATION OF THE MEMBRANE TRIAL STRESS TENSOR
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            ! -> Membrane forces
            signxx(i) = sigoxx(i) + lambda_m*(depsxx(i) + depsyy(i)) +         &
                                                  two*mu_m*depsxx(i)
            signyy(i) = sigoyy(i) + lambda_m*(depsxx(i) + depsyy(i)) +         &
                                                  two*mu_m*depsyy(i)
            signxy(i) = sigoxy(i) +  mu_m*depsxy(i)
            signyz(i) = sigoyz(i) + gs(i)*depsyz(i)
            signzx(i) = sigozx(i) + gs(i)*depszx(i)
            ! -> Hourglass stabilization variable
            et(i) = one
            ! -> Sound speed
            ssp(i) = sqrt((young/(one - nu*nu))/rho(i))
          enddo
!
          !=====================================================================
          !< - COMPUTATION OF THE BENDING DAMAGE AND TRIAL BENDING MOMENTS
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
!
            !< Computation of thermodynamic forces driving the bending damage
            !<------------------------------------------------------------------
            !< Trace of the bending curvature increment tensor
            tr_kb = epbxx(i) + epbyy(i)
            !< Eigenvalues of the bending curvature tensor
            center = (epbxx(i) + epbyy(i))/two
            radius = sqrt(((epbxx(i) - epbyy(i))*half)**2 + (epbxy(i)*half)**2)
            kappa1 = center + radius
            kappa2 = center - radius
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
!
            !< Damage evolution and computation of the bending moments
            !< -----------------------------------------------------------------
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
!
            !< Computation of the damaged bending stiffness matrix
            !< -----------------------------------------------------------------
            !< D_b^(1) : constant damage contribution to the tangent stiffness 
            Db11(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                     &
                                             xi_kap1(i)*dkap1dkxx*dkap1dkxx +  &
                                             xi_kap2(i)*dkap2dkxx*dkap2dkxx)
            Db12(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                     &
                                             xi_kap1(i)*dkap1dkxx*dkap1dkyy +  &
                                             xi_kap2(i)*dkap2dkxx*dkap2dkyy)
            Db21(i) = Db12(i)
            Db22(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                     &
                                             xi_kap1(i)*dkap1dkyy*dkap1dkyy +  &
                                             xi_kap2(i)*dkap2dkyy*dkap2dkyy)
            Db33(i) = two*mu_b(i)*(xi_kap1(i)*dkap1dkxy*dkap1dkxy +            &
                                   xi_kap2(i)*dkap2dkxy*dkap2dkxy)
            !< D_b^(2) : variable damage contribution to the tangent stiffness
            !< Gradient of W+ w.r.t. principal curvatures
            dW1_dk1 = zero
            dW1_dk2 = zero
            if (tr_kb  >= zero) then
              dW1_dk1 = dW1_dk1 + lambda_b(i) * tr_kb
              dW1_dk2 = dW1_dk2 + lambda_b(i) * tr_kb
            end if
            if (kappa1 >= zero) dW1_dk1 = dW1_dk1 + two*mu_b(i)*kappa1
            if (kappa2 >= zero) dW1_dk2 = dW1_dk2 + two*mu_b(i)*kappa2
            !< Gradient of W- w.r.t. principal curvatures
            dW2_dk1 = zero
            dW2_dk2 = zero
            if (tr_kb  < zero) then
              dW2_dk1 = dW2_dk1 + lambda_b(i) * tr_kb
              dW2_dk2 = dW2_dk2 + lambda_b(i) * tr_kb
            end if
            if (kappa1 < zero) dW2_dk1 = dW2_dk1 + two*mu_b(i) * kappa1
            if (kappa2 < zero) dW2_dk2 = dW2_dk2 + two*mu_b(i) * kappa2
            !< Gradient of W+/W- back in the reinforcement frame via chain rule
            !< dW/dkab = dW/dk1 * dk1/dkab + dW/dk2 * dk2/dkab
            dW1_dkxx = dW1_dk1*dkap1dkxx + dW1_dk2*dkap2dkxx
            dW1_dkyy = dW1_dk1*dkap1dkyy + dW1_dk2*dkap2dkyy
            dW1_dkxy = dW1_dk1*dkap1dkxy + dW1_dk2*dkap2dkxy
            dW2_dkxx = dW2_dk1*dkap1dkxx + dW2_dk2*dkap2dkxx
            dW2_dkyy = dW2_dk1*dkap1dkyy + dW2_dk2*dkap2dkyy
            dW2_dkxy = dW2_dk1*dkap1dkxy + dW2_dk2*dkap2dkxy
            !< Scalar prefactors: (gamma-1) / (2*(1+D)*W)
            !< Guard against W=0 (no damage active -> term vanishes anyway)
            fac1 = zero
            if (phi_b_pos > em20) fac1 = (gamma - one) / (two*(one + dmg(i,2))*phi_b_pos)
            fac2 = zero
            if (phi_b_neg > em20) fac2 = (gamma - one) / (two*(one + dmg(i,3))*phi_b_neg)
            !< Assemble D_b^(2) as rank-1 updates (outer products)
            !< Only active if the damage criterion is active (loading)
            Db11(i) = Db11(i) + fac1*dW1_dkxx*dW1_dkxx + fac2*dW2_dkxx*dW2_dkxx
            Db12(i) = Db12(i) + fac1*dW1_dkxx*dW1_dkyy + fac2*dW2_dkxx*dW2_dkyy
            Db21(i) = Db12(i)
            Db22(i) = Db22(i) + fac1*dW1_dkyy*dW1_dkyy + fac2*dW2_dkyy*dW2_dkyy
            Db33(i) = Db33(i) + fac1*dW1_dkxy*dW1_dkxy + fac2*dW2_dkxy*dW2_dkxy
          enddo     
!
          !=====================================================================
          !< - COMPUTATION OF THE LIMIT MOMENT FOR THE PLASTICITY CRITERIA
          !=====================================================================
          !< Reset tables
          mfx_pos(1:nel) = zero
          mfx_neg(1:nel) = zero
          mfy_pos(1:nel) = zero
          mfy_neg(1:nel) = zero
          allocate(rho_ext(0:nlayer+1))
!
          !< Critical moment for the plasticity criterion in bending for the 
          !  reinforcement in x direction
          !---------------------------------------------------------------------
          rho_ext(0) = -one
          do k = 1, nlayer
            rho_ext(k) = rho_x(k)
          enddo
          rho_ext(nlayer + 1) = one
          !< Loop over the elements
          do i = 1, nel
            na = two/(thk0(i)*f_c)
            !-------------------------------------------------------------------
            !< Positive bending moment 
            !-------------------------------------------------------------------
            n_f = two
            do j = 0, nlayer
              S_j = zero
              do k = 1, nlayer
                if (rho_x(k) < rho_ext(j)) then
                  xi_k =  one
                else
                  xi_k = -one
                endif
                S_j = S_j + sig_y(k) * omega_x(k) * xi_k
              enddo
              n_star = one + na * (signxx(i) - S_j)
              if ((n_star >= rho_ext(j)).and.(n_star < rho_ext(j+1))) then
                n_f = n_star
                exit
              endif
            enddo
            mfx_pos(i) = (thk0(i)**2 *f_c/four)*(one - n_f**2)
            dmfx_pos(i) = -thk0(i) * n_f
            do k = 1, nlayer
              mfx_pos(i) = mfx_pos(i)                                          &
                         + sig_y(k)*omega_x(k)*sign(one,n_f-rho_x(k))*         &
                           rho_x(k)*thk0(i)*half     
            enddo  
            !-------------------------------------------------------------------
            !< Negative bending moment
            !-------------------------------------------------------------------
            n_f = two
            do j = 0, nlayer
              S_j = zero
              do k = 1, nlayer
                if (rho_x(k) < rho_ext(j)) then
                  xi_k =  one
                else
                  xi_k = -one
                end if
                S_j = S_j + sig_y(k) * omega_x(k) * xi_k
              end do
              n_star = one + na * (-signxx(i) - S_j)
              if ((n_star >= rho_ext(j)) .and. &
                  (n_star <  rho_ext(j + 1))) then
                n_f = n_star
                exit
              end if
            end do
            mfx_neg(i) = -(thk0(i)**2 * f_c / four) * (one - n_f**2)
            dmfx_neg(i) = thk0(i) * n_f
            do k = 1, nlayer
              mfx_neg(i) = mfx_neg(i)                                          &
                         - sig_y(k)*omega_x(k)*sign(one,n_f-rho_x(k))*         &
                           rho_x(k)*thk0(i)*half
            end do                        
          enddo
!
          rho_ext(0) = -one
          do k = 1, nlayer
            rho_ext(k) = rho_y(k)
          end do
          rho_ext(nlayer + 1) = one
          do i = 1,nel
            na = two / (thk0(i)*f_c)
            !-------------------------------------------------------------------
            !< Positive bending moment 
            !-------------------------------------------------------------------
            n_f = two
            do j = 0, nlayer
              S_j = zero
              do k = 1, nlayer
                if (rho_y(k) < rho_ext(j)) then
                  xi_k =  one
                else
                  xi_k = -one
                end if
                S_j = S_j + sig_y(k) * omega_y(k) * xi_k
              end do
              n_star = one + na * (signyy(i) - S_j)
              if ((n_star >= rho_ext(j)) .and. &
                  (n_star <  rho_ext(j + 1))) then
                n_f = n_star
                exit
              end if
            end do
            mfy_pos(i) = (thk0(i)**2 * f_c / four) * (one - n_f**2)
            dmfy_pos(i) = -thk0(i) * n_f
            do k = 1, nlayer
              mfy_pos(i) = mfy_pos(i)                                          &
                         + sig_y(k)*omega_y(k)*sign(one,n_f-rho_y(k))          &
                         * rho_y(k) * thk0(i) * half
            end do
            !------------------------------------------------------------
            !< Negative bending moment
            !------------------------------------------------------------
            n_f = two
            do j = 0, nlayer
              S_j = zero
              do k = 1, nlayer
                if (rho_y(k) < rho_ext(j)) then
                  xi_k =  one
                else
                  xi_k = -one
                end if
                S_j = S_j + sig_y(k) * omega_y(k) * xi_k
              end do
              n_star = one + na * (-signyy(i) - S_j)
              if ((n_star >= rho_ext(j)) .and. &
                  (n_star <  rho_ext(j + 1))) then
                n_f = n_star
                exit
              end if
            end do
            mfy_neg(i) = -(thk0(i)**2 * f_c / four) * (one - n_f**2)
            dmfy_neg(i) = thk0(i) * n_f
            do k = 1, nlayer
              mfy_neg(i) = mfy_neg(i)                                          &
                         - sig_y(k) * omega_y(k)                               &
                         * sign(one, n_f - rho_y(k))                           &
                         * rho_y(k) * thk0(i) * half
            end do
          enddo
          deallocate(rho_ext)
!
          !=====================================================================
          !< - Computation of the membrane-bending plasticity
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            !< Computation of the plasticity criteria
            f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +      &
                                                   momnxy(i)*momnxy(i)
            f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +      &
                                                   momnxy(i)*momnxy(i)
            active_I  = (f1p(i)  > zero)
            active_II = (f2p(i)  > zero)
          enddo
!
          !< Index of the active criteria
          nindx  = 0
          nindx2 = 0
          nindx3 = 0
          do i = 1, nel
            if (f1p(i) > zero .and. f2p(i) > zero) then
              indx(nindx) = i
              nindx = nindx + 1
            else if (f1p(i) > zero) then
              indx2(nindx2) = i
              nindx2 = nindx2 + 1
            else if (f2p(i) > zero) then
              indx3(nindx3) = i
              nindx3 = nindx3 + 1
            endif
          enddo
!
          !< Both criteria active : loop over the indices of the active criteria to compute the plastic correction
          if (nindx > 0) then
            do iter = 1,3
              do ii = 1, nindx
                i = indx(ii)
                !< Gradients f_I
                dfI_dMx  = -(momnyy(i) - mfy_pos(i))
                dfI_dMy  = -(momnxx(i) - mfx_pos(i))
                dfI_dMxy =  two * momnxy(i)
                dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dmfx_pos(i)
                dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dmfy_pos(i)
                !< Gradients f_II
                dfII_dMx  = -(momnyy(i) - mfy_neg(i))
                dfII_dMy  = -(momnxx(i) - mfx_neg(i))
                dfII_dMxy =  two * momnxy(i)
                dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dmfx_neg(i)
                dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dmfy_neg(i)
                !< Système 2x2 : [A]{lam} = {b}
                !< A_11 = (df_I/dN)^T(Dm+Cm)(df_I/dN) + (df_I/dM)^T(Db+Cb)(df_I/dM)
                A11 = dfI_dNx**2*Dm11 + two*dfI_dNx*dfI_dNy*Dm12               &
                    + dfI_dNy**2*Dm22                                          &
                    + dfI_dMx**2*Db11(i) + two*dfI_dMx*dfI_dMy*Db12(i)         &
                    + dfI_dMy**2*Db22(i) + dfI_dMxy**2*Db33(i)
                A12 = dfI_dNx*dfII_dNx*Dm11                                    &
                    + (dfI_dNx*dfII_dNy+dfI_dNy*dfII_dNx)*Dm12                 &
                    + dfI_dNy*dfII_dNy*Dm22                                    &
                    + dfI_dMx*dfII_dMx*Db11(i)                                 &
                    + (dfI_dMx*dfII_dMy+dfI_dMy*dfII_dMx)*Db12(i)              &
                    + dfI_dMy*dfII_dMy*Db22(i)                                 &
                    + dfI_dMxy*dfII_dMxy*Db33(i)
                A21 = A12
                A22 = dfII_dNx**2*Dm11 + two*dfII_dNx*dfII_dNy*Dm12            &
                    + dfII_dNy**2*Dm22                                         &
                    + dfII_dMx**2*Db11(i) + two*dfII_dMx*dfII_dMy*Db12(i)      &
                    + dfII_dMy**2*Db22(i) + dfII_dMxy**2*Db33(i)
                !< Seconds membres
                b1 = fI  + dfI_dNx*(Dm11*depsxx(i)+Dm12*depsyy(i))             &
                         + dfI_dNy*(Dm21*depsxx(i)+Dm22*depsyy(i))             &
                         + dfI_dMx*(Db11(i)*depbxx(i)+Db12(i)*depbyy(i)+       &
                                                      Db13(i)*depbxy(i))       &
                         + dfI_dMy*(Db21(i)*depbxx(i)+Db22(i)*depbyy(i)+       &
                                                      Db23(i)*depbxy(i))       &
                         + dfI_dMxy*(Db31(i)*depbxx(i)+Db32(i)*depbyy(i)+      &
                                                       Db33(i)*depbxy(i))
              
                b2 = fII + dfII_dNx*(Dm11*depsxx(i)+Dm12*depsyy(i))            &
                         + dfII_dNy*(Dm21*depsxx(i)+Dm22*depsyy(i))            &
                         + dfII_dMx*(Db11(i)*depbxx(i)+Db12(i)*depbyy(i)+      &
                                                       Db13(i)*depbxy(i))      &
                         + dfII_dMy*(Db21(i)*depbxx(i)+Db22(i)*depbyy(i)+      &
                                                       Db23(i)*depbxy(i))      &
                         + dfII_dMxy*(Db31(i)*depbxx(i)+Db32(i)*depbyy(i)+     &
                                                        Db33(i)*depbxy(i))
                !< Résolution 2x2 par Cramer
                det2 = A11*A22 - A12*A21
                lam_I  = (b1*A22 - b2*A12) / det2
                lam_II = (A11*b2 - A21*b1) / det2
                !< Incréments plastiques combinés
                deps0p_x = lam_I*dfI_dNx  + lam_II*dfII_dNx
                deps0p_y = lam_I*dfI_dNy  + lam_II*dfII_dNy
                dkp_x    = lam_I*dfI_dMx  + lam_II*dfII_dMx
                dkp_y    = lam_I*dfI_dMy  + lam_II*dfII_dMy
                dkp_xy   = lam_I*dfI_dMxy + lam_II*dfII_dMxy
                !< Update of the membrane stresses and bending moments
                signxx(i) = signxx(i) - Dm11*deps0p_x - Dm12*deps0p_y
                signyy(i) = signyy(i) - Dm21*deps0p_x - Dm22*deps0p_y
                momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y - Db13(i)*dkp_xy
                momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y - Db23(i)*dkp_xy
                momnxy(i) = momnxy(i) - Db31(i)*dkp_x - Db32(i)*dkp_y - Db33(i)*dkp_xy
                !< Update the limit moments for the plasticity criteria
                
                !< Update the plasticity criteria
                f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) + 
                                                        momnxy(i)*momnxy(i)
                f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) + 
                                                        momnxy(i)*momnxy(i)
              enddo
            enddo
          endif

          !< Only criterion I active
          if (nindx2 > 0) then
            do iter = 1,3
              do ii = 1, nindx2
                i = indx2(ii)
  
                !< Gradients df_I/dM et df_I/dN (cutting plane : évalués en début)
                dfI_dMx  = -(momnyy(i) - mfy_pos(i))
                dfI_dMy  = -(momnxx(i) - mfx_pos(i))
                dfI_dMxy =  two * momnxy(i)
                dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dmfx_pos(i)
                dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dmfy_pos(i)
    
                num = f1p(i)                                                     
                    ! + dfI_dNx  * (Dm11*depsxx(i) + Dm12*depsyy(i))               &
                    ! + dfI_dNy  * (Dm21*depsxx(i) + Dm22*depsyy(i))               &
                    ! + dfI_dMx  * (Db11(i)*depbxx(i) + Db12(i)*depbyy(i) +        &
                    !                                   Db13(i)*depbxy(i))         &
                    ! + dfI_dMy  * (Db21(i)*depbxx(i) + Db22(i)*depbyy(i) +        &
                    !                                   Db23(i)*depbxy(i))         &
                    ! + dfI_dMxy * (Db31(i)*depbxx(i) + Db32(i)*depbyy(i) +        &
                    !                                   Db33(i)*depbxy(i)) 
                    
                !< Dénominateur de lambda_p
                !< terme membrane : (df/dN)^T (Dm + Cm) (df/dN)
                den = dfI_dNx**2 * Dm11                                          &
                    + two*dfI_dNx*dfI_dNy*Dm12                                   &
                    + dfI_dNy**2 * Dm22                                          &
                !< terme flexion : (df/dM)^T (Db + Cb) (df/dM)
                    + dfI_dMx**2  * Db11(i)                                      &
                    + two*dfI_dMx*dfI_dMy*Db12(i)                                &
                    + dfI_dMy**2  * Db22(i)                                      &
                    + dfI_dMxy**2 * Db33(i)
              
                lam_I = num / den            
    
                !< Correction des déformations plastiques (loi de normalité, éq. 6.13)
                depsp_xx  = lam_I * dfI_dNx
                depsp_yy  = lam_I * dfI_dNy
                dkp_x     = lam_I * dfI_dMx
                dkp_y     = lam_I * dfI_dMy
                dkp_xy    = lam_I * dfI_dMxy
  
                !< Update of the membrane stresses and bending moments
                signxx(i) = signxx(i) - Dm11*deps0p_x - Dm12*deps0p_y
                signyy(i) = signyy(i) - Dm21*deps0p_x - Dm22*deps0p_y
                momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y - Db13(i)*dkp_xy
                momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y - Db23(i)*dkp_xy
                momnxy(i) = momnxy(i) - Db31(i)*dkp_x - Db32(i)*dkp_y - Db33(i)*dkp_xy
                !< Update the limit moments for the plasticity criteria
                
                !< Update the plasticity criteria
                f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) + 
                                                        momnxy(i)*momnxy(i)
  
  
              enddo
            enddo
          endif

          !< Only criterion II active
          if (nindx3 > 0) then
            do iter = 1,3
              do ii = 1, nindx3
                i = indx3(ii)
                !< Gradients df_II/dM et df_II/dN (cutting plane : évalués en début)
                dfII_dMx  = -(momnyy(i) - mfy_neg(i))
                dfII_dMy  = -(momnxx(i) - mfx_neg(i))
                dfII_dMxy =  two * momnxy(i)
                dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dmfx_neg(i)
                dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dmfy_neg(i)
    
                num = fII + dfII_dNx*(Dm11*depsxx(i)+Dm12*depsyy(i))              &
                          + dfII_dNy*(Dm21*depsxx(i)+Dm22*depsyy(i))              &
                          + dfII_dMx*(Db11(i)*depbxx(i)+Db12(i)*depbyy(i)+        &
                                                        Db13(i)*depbxy(i))        &
                          + dfII_dMy*(Db21(i)*depbxx(i)+Db22(i)*depbyy(i)+        &
                                                        Db23(i)*depbxy(i))        &
                          + dfII_dMxy*(Db31(i)*depbxx(i)+Db32(i)*depbyy(i)+       &
                                                          Db33(i)*depbxy(i))
                    
                !< Dénominateur de lambda_p
                !< terme membrane : (df/dN)^T (Dm + Cm) (df/dN)
                den = dfII_dNx**2 * Dm11                                          &
                    + two*dfII_dNx*dfII_dNy*Dm12                                   &
                    + dfII_dNy**2 * Dm22                                          &
                !< terme flexion : (df/dM)^T (Db + Cb) (df/dM)
                    + dfII_dMx**2  * Db11(i)                                      &
                    + two*dfII_dMx*dfII_dMy*Db12(i)                                &
                    + dfII_dMy**2  * Db22(i)                                      &
                    + dfII_dMxy**2 * Db33(i)
              
                lam_II = num / den            
    
                !< Correction des déformations plastiques (loi de normalité, éq. 6.13)
                depsp_xx  = lam_II * dfII_dNx
                depsp_yy  = lam_II * dfII_dNy
                dkp_x     = lam_II * dfII_dMx
                dkp_y     = lam_II * dfII_dMy
                dkp_xy    = lam_II * dfII_dMxy
  
                !< Update of the membrane stresses and bending moments
                signxx(i) = signxx(i) - Dm11*deps0p_x - Dm12*deps0p_y
                signyy(i) = signyy(i) - Dm21*deps0p_x - Dm22*deps0p_y
                momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y - Db13(i)*dkp_xy
                momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y - Db23(i)*dkp_xy
                momnxy(i) = momnxy(i) - Db31(i)*dkp_x - Db32(i)*dkp_y - Db33(i)*dkp_xy
                !< Update the limit moments for the plasticity criteria
                
                !< Update the plasticity criteria
                f2p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) + 
                                                        momnxy(i)*momnxy(i)
  
  
              enddo
            enddo
          endif           
!
          !=====================================================================
        end subroutine sigeps136g
      end module sigeps136g_mod
