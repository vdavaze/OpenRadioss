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
      subroutine sigeps136g(                                                   &
        nel    ,matparam,rho    ,dmg     ,thk0   ,                             &
        depsxx ,depsyy ,depsxy  ,depsyz  ,depszx ,                             &                             
        depbxx ,depbyy ,depbxy  ,                                              &
        sigoxx ,sigoyy ,sigoxy  ,sigoyz  ,sigozx ,                             &
        signxx ,signyy ,signxy  ,signyz  ,signzx ,                             &
        momnxx ,momnyy ,momnxy  ,                                              &
        ssp    ,et     ,gs      ,nuvar   ,uvar   ,                             &
        shf    ,pla    ,sigb    )
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
        real(kind=WP), dimension(nel), intent(inout) :: pla      !< Array of plastic strains for post-processing and output
        real(kind=WP), dimension(nel,5),intent(inout) :: sigb  !< Array of back stresses for post-processing and output
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,k,ii,j,nindx,nindx2,nindx3,indx(nel),indx2(nel),          &
          indx3(nel), iter
        real(kind=WP) :: young,nu,shear,lambda_m,mu_m,f_t,f_c,gamma,           &
          dmax1,dmax2,Dm11,Dm12,cm,cb
        real(kind=WP), dimension(nel) :: lambda_b, mu_b, k0_1, k0_2
        real(kind=WP) :: center, radius, kappa1, kappa2, tr_kb, Y1,            &
          Y2,d1_trial,d2_trial,A11,A12,A22,dW1_dk1,dW1_dk2,dW2_dk1,            &
          dW2_dk2,dW1_dkxx,dW1_dkyy,dW1_dkxy,dW2_dkxx,dW2_dkyy,dW2_dkxy,       &
          fac1,fac2,A21,dkap1dkxx,dkap1dkyy,dkap1dkxy,dkap2dkxx,dkap2dkyy,     &
          dkap2dkxy,phi_b_pos,phi_b_neg,b1,b2,den,depsp_x,depsp_y,dkp_x,       &
          dkp_y,dkp_xy,det2,dfI_dNx,dfI_dNy,dfI_dMx,dfI_dMy,dfI_dMxy,          &
          dfII_dNx,dfII_dNy,dfII_dMx,dfII_dMy,dfII_dMxy,lam_I,lam_II
        real(kind=WP), dimension(nel) :: xi_tr, epbxx, epbyy, epbxy, xi_kap1,  &
          xi_kap2, Db11, Db12, Db21, Db22, Db33, mfx_pos, mfx_neg, mfy_pos,    &
          mfy_neg, dmfx_pos, dmfx_neg, dmfy_pos, dmfy_neg, f1p, f2p, g1p, g2p
        real(kind=WP), dimension(2) :: rho_x, rho_y, sig_y, omega_x, omega_y
        real(kind=WP), parameter :: tol = 1e-3
        real(kind=WP), parameter :: xi_pos = 1.0d0
        real(kind=WP), parameter :: xi_neg = -1.0d0
        integer, parameter :: nmax = 100
!
        !=======================================================================
        !< - Initialisation of computation on time step
        !=======================================================================
!
        !< Recovering real model parameters
        young      = matparam%young        !< Concrete Young's modulus
        nu         = matparam%nu           !< Concrete Poisson's ratio
        shear      = matparam%shear        !< Concrete shear modulus
        lambda_m   = matparam%uparam(1)    !< Membrane Lamé parameter
        mu_m       = matparam%uparam(2)    !< Membrane shear modulus
        f_t        = matparam%uparam(3)    !< Concrete Tensile strength
        f_c        = matparam%uparam(4)    !< Concrete Compressive strength
        gamma      = matparam%uparam(5)    !< Concrete 
        dmax1      = matparam%uparam(6)    !< Maximum damage for positive bending
        dmax2      = matparam%uparam(7)    !< Maximum damage for negative bending
        sig_y(1)   = matparam%uparam(8)    !< Yield stress of the reinforcement in lower layer
        omega_x(1) = matparam%uparam(9)    !< Reinforcement ratio in x direction for lower layer
        omega_y(1) = matparam%uparam(10)   !< Reinforcement ratio in y direction for lower layer
        rho_x(1)   = matparam%uparam(11)   !< Position in thickness of reinforcement in x direction for lower layer
        rho_y(1)   = matparam%uparam(12)   !< Position in thickness of reinforcement in y direction for lower layer
        sig_y(2)   = matparam%uparam(13)   !< Yield stress of the reinforcement in upper layer
        omega_x(2) = matparam%uparam(14)   !< Reinforcement ratio in x direction for upper layer
        omega_y(2) = matparam%uparam(15)   !< Reinforcement ratio in y direction for upper layer
        rho_x(2)   = matparam%uparam(16)   !< Position in thickness of reinforcement in x direction for upper layer
        rho_y(2)   = matparam%uparam(17)   !< Position in thickness of reinforcement in y direction for upper layer
        cm         = matparam%uparam(18)   !< Prager hardening parameter for membrane forces
        cb         = matparam%uparam(19)   !< Prager hardening parameter for bending moments
!
        !< Computation of some real parameters
        Dm11 = lambda_m + two*mu_m                         !< Membrane stiffness D11 = lambda + 2*mu
        Dm12 = lambda_m                                    !< Membrane stiffness D12 = lambda
        gs(1:nel) = shear*shf(1:nel)                       !< Correction factor for transverse shear 
        lambda_b(1:nel) = lambda_m*thk0(1:nel)*one_over_12 !< Bending Lamé parameter
        mu_b(1:nel) = mu_m*thk0(1:nel)*one_over_12         !< Bending shear modulus
        !< Bending damage threshold in positive bending (tension inner face)   
        k0_1(1:nel) = -(gamma-one)*half*(lambda_b(1:nel) + two*mu_b(1:nel))*   &
                       ((two*f_t*(one - nu*nu))/(young*thk0(1:nel)))**2
        !< Bending damage threshold in negative bending (tension outer face)
        k0_2(1:nel) = -(gamma-one)*half*(lambda_b(1:nel) + two*mu_b(1:nel))*   &
                       ((two*f_t*(one - nu*nu))/(young*thk0(1:nel)))**2
!
        !< Recompute the bending curvature tensor at current time step
        do i = 1, nel
          uvar(i,1) = uvar(i,1) + depbxx(i)
          uvar(i,2) = uvar(i,2) + depbyy(i)
          uvar(i,3) = uvar(i,3) + depbxy(i)
          epbxx(i)  = uvar(i,1) - uvar(i,6)
          epbyy(i)  = uvar(i,2) - uvar(i,7)
          epbxy(i)  = uvar(i,3) - uvar(i,8)
        enddo
!
        !=======================================================================
        !< - COMPUTATION OF THE MEMBRANE TRIAL STRESS TENSOR
        !=======================================================================
        !< Loop over elements
        do i = 1, nel
          ! -> Membrane forces
          signxx(i) = sigoxx(i) +  Dm11*depsxx(i) + Dm12*depsyy(i)                                      
          signyy(i) = sigoyy(i) +  Dm12*depsxx(i) + Dm11*depsyy(i)
          signxy(i) = sigoxy(i) +  mu_m*depsxy(i)
          signyz(i) = sigoyz(i) + gs(i)*depsyz(i)
          signzx(i) = sigozx(i) + gs(i)*depszx(i)
          ! -> Hourglass stabilization variable
          et(i) = one
          ! -> Sound speed
          ssp(i) = sqrt((young/(one - nu*nu))/rho(i))
        enddo
!
        !=======================================================================
        !< - COMPUTATION OF THE BENDING DAMAGE AND TRIAL BENDING MOMENTS
        !=======================================================================
        !< Loop over elements
        do i = 1, nel
!
          !< Computation of thermodynamic forces driving the bending damage
          !<--------------------------------------------------------------------
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
          !< -------------------------------------------------------------------
          !< Evolution law of the bending damage
          ! --> Positive bending, inner face
          if (Y1 >= k0_1(i)) then
            d1_trial = sqrt(phi_b_pos / k0_1(i)) - one
            d1_trial = max(zero, min(d1_trial, dmax1))
            dmg(i,2) = max(dmg(i,2), d1_trial)
          end if
          ! --> Negative bending, outer face
          if (Y2 >= k0_2(i)) then
            d2_trial = sqrt(phi_b_neg / k0_2(i)) - one
            d2_trial = max(zero, min(d2_trial, dmax2))
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
          !< Save maximum damage for output and post-processing
          dmg(i,1) = max(dmg(i,2), dmg(i,3)) 
          !< Update the bending moments with the softening factors
          dkap1dkxx = half + (epbxx(i) - epbyy(i))/(four*max(radius,em20))
          dkap1dkyy = half - (epbxx(i) - epbyy(i))/(four*max(radius,em20))
          dkap1dkxy = (epbxy(i)*half)/(two*max(radius,em20))
          dkap2dkxx = half - (epbxx(i) - epbyy(i))/(four*max(radius,em20))
          dkap2dkyy = half + (epbxx(i) - epbyy(i))/(four*max(radius,em20))
          dkap2dkxy = -(epbxy(i)*half)/(two*max(radius,em20))
          momnxx(i) = lambda_b(i)*tr_kb*xi_tr(i) + two*mu_b(i)*(               &
                          xi_kap1(i)*kappa1*dkap1dkxx +                        &
                          xi_kap2(i)*kappa2*dkap2dkxx ) 
          momnyy(i) = lambda_b(i)*tr_kb*xi_tr(i) + two*mu_b(i)*(               & 
                          xi_kap1(i)*kappa1*dkap1dkyy +                        &
                          xi_kap2(i)*kappa2*dkap2dkyy )
          momnxy(i) = two*mu_b(i)*(                                            &
                          xi_kap1(i)*kappa1*dkap1dkxy +                        &
                          xi_kap2(i)*kappa2*dkap2dkxy )
!
          !< Computation of the damaged bending stiffness matrix
          !< -------------------------------------------------------------------
          !< D_b^(1) : constant damage contribution to the tangent stiffness 
          Db11(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                       &
                                           xi_kap1(i)*dkap1dkxx*dkap1dkxx +    &
                                           xi_kap2(i)*dkap2dkxx*dkap2dkxx)
          Db12(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                       &
                                           xi_kap1(i)*dkap1dkxx*dkap1dkyy +    &
                                           xi_kap2(i)*dkap2dkxx*dkap2dkyy)
          Db21(i) = Db12(i)
          Db22(i) = lambda_b(i)*xi_tr(i) + two*mu_b(i)*(                       &
                                           xi_kap1(i)*dkap1dkyy*dkap1dkyy +    &
                                           xi_kap2(i)*dkap2dkyy*dkap2dkyy)
          Db33(i) = two*mu_b(i)*(xi_kap1(i)*dkap1dkxy*dkap1dkxy +              &
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
          if (phi_b_pos > em20) fac1 = (gamma - one) /                         &
                                       (two*(one + dmg(i,2))*phi_b_pos)
          fac2 = zero
          if (phi_b_neg > em20) fac2 = (gamma - one) /                         &
                                       (two*(one + dmg(i,3))*phi_b_neg)
          !< Assemble D_b^(2) as rank-1 updates (outer products)
          !< Only active if the damage criterion is active (loading)
          Db11(i) = Db11(i) + fac1*dW1_dkxx*dW1_dkxx + fac2*dW2_dkxx*dW2_dkxx
          Db12(i) = Db12(i) + fac1*dW1_dkxx*dW1_dkyy + fac2*dW2_dkxx*dW2_dkyy
          Db21(i) = Db12(i)
          Db22(i) = Db22(i) + fac1*dW1_dkyy*dW1_dkyy + fac2*dW2_dkyy*dW2_dkyy
          Db33(i) = Db33(i) + fac1*dW1_dkxy*dW1_dkxy + fac2*dW2_dkxy*dW2_dkxy
        enddo     
!
        !=======================================================================
        !< - COMPUTATION OF THE LIMIT MOMENT FOR THE PLASTICITY CRITERIA
        !=======================================================================
        !< Add kinematic hardening contribution 
        do i = 1, nel
          signxx(i) = signxx(i) - sigb(i,1)
          signyy(i) = signyy(i) - sigb(i,2)
          momnxx(i) = momnxx(i) - sigb(i,3)
          momnyy(i) = momnyy(i) - sigb(i,4)
          momnxy(i) = momnxy(i) - sigb(i,5)
        enddo 
!
        !< Reset tables
        mfx_pos(1:nel) = zero
        mfx_neg(1:nel) = zero
        mfy_pos(1:nel) = zero
        mfy_neg(1:nel) = zero
!
        !< Critical moment for the plasticity criterion in bending for the 
        !  reinforcement in x direction
        !-----------------------------------------------------------------------
        !< Loop over the elements
        do i = 1, nel
          !< Positive bending moment
          call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,       &
            rho_y,zero,xi_pos,mfx_pos(i),dmfx_pos(i))
          !< Negative bending moment
          call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,       &
            rho_y,zero,xi_neg,mfx_neg(i),dmfx_neg(i))
        enddo
!
        !< Critical moment for the plasticity criterion in bending for the 
        !  reinforcement in y direction
        !-----------------------------------------------------------------------        
        do i = 1,nel
          !< Positive bending moment
          call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,       &
            rho_y,pi*half,xi_pos,mfy_pos(i),dmfy_pos(i))
          !< Negative bending moment
          call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,       &
            rho_y,pi*half,xi_neg,mfy_neg(i),dmfy_neg(i))
        enddo
!
        !=======================================================================
        !< - COMPUTATION OF THE MEMBRANE-BENDING PLASTICITY
        !=======================================================================
        !< Loop over elements
        do i = 1, nel
          !< Computation of the plasticity criteria
          f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +        &
                                                 momnxy(i)*momnxy(i)
          f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +        &
                                                 momnxy(i)*momnxy(i)
          g1p(i) =   momnxx(i) - mfx_pos(i) + momnyy(i) - mfy_pos(i)
          g2p(i) = -(momnxx(i) - mfx_neg(i) + momnyy(i) - mfy_neg(i))
        enddo
!
        !< Index of the active criteria
        nindx  = 0
        nindx2 = 0
        nindx3 = 0
        do i = 1, nel
          if ((f1p(i) > zero) .and. (f2p(i) > zero) .and.                      &
              (g1p(i) > zero) .and. (g2p(i) > zero)) then
            nindx = nindx + 1
            indx(nindx) = i
          else if ((f1p(i) > zero) .and. (g1p(i) > zero)) then
            nindx2 = nindx2 + 1
            indx2(nindx2) = i
          else if ((f2p(i) > zero) .and. (g2p(i) > zero)) then
            nindx3 = nindx3 + 1
            indx3(nindx3) = i
          endif
        enddo
!
        !-----------------------------------------------------------------------
        !< Both criteria active : loop over the indices of the active criteria 
        !  to compute the plastic correction
        !-----------------------------------------------------------------------
        if (nindx > 0) then
          do ii = 1, nindx
            iter = 0
            i = indx(ii)
            do while ((abs(f1p(i)) > tol .or. abs(f2p(i)) > tol) .and.         &
                      (iter < nmax))
              iter = iter + 1
              !< Normal derivatives of the criteria f_I 
              dfI_dMx  = -(momnyy(i) - mfy_pos(i))
              dfI_dMy  = -(momnxx(i) - mfx_pos(i))
              dfI_dMxy =  two * momnxy(i)
              dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dmfx_pos(i)
              dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dmfy_pos(i)
              !< Normal derivatives of the criteria f_II
              dfII_dMx  = -(momnyy(i) - mfy_neg(i))
              dfII_dMy  = -(momnxx(i) - mfx_neg(i))
              dfII_dMxy =  two * momnxy(i)
              dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dmfx_neg(i)
              dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dmfy_neg(i)
              !< Système 2x2 : [A]{lam} = {b}
              !< A_11 = (df_I/dN)^T(Dm+Cm)(df_I/dN) + (df_I/dM)^T(Db+Cb)(df_I/dM)
              A11 = dfI_dNx*(dfI_dNx*(Dm11 + cm) + dfI_dNy*Dm12) +           &
                    dfI_dNy*(dfI_dNx*Dm12 + dfI_dNy*(Dm11 + cm)) +           &
                    dfI_dMx*(dfI_dMx*(Db11(i) + cb) + dfI_dMy*Db12(i)) +     &
                    dfI_dMy*(dfI_dMx*Db21(i) + dfI_dMy*(Db22(i) + cb)) +     &
                    dfI_dMxy*(dfI_dMxy*(Db33(i) + cb))
              A12 = dfI_dNx*(dfII_dNx*(Dm11 + cm) + dfII_dNy*Dm12) +         &
                    dfI_dNy*(dfII_dNx*Dm12 + dfII_dNy*(Dm11 + cm)) +         &
                    dfI_dMx*(dfII_dMx*(Db11(i) + cb) + dfII_dMy*Db12(i)) +   &
                    dfI_dMy*(dfII_dMx*Db21(i) + dfII_dMy*(Db22(i) + cb)) +   &
                    dfI_dMxy*(dfII_dMxy*(Db33(i) + cb))
              A21 = A12
              A22 = dfII_dNx*(dfII_dNx*(Dm11 + cm) + dfII_dNy*Dm12) +        &
                    dfII_dNy*(dfII_dNx*Dm12 + dfII_dNy*(Dm11 + cm)) +        &
                    dfII_dMx*(dfII_dMx*(Db11(i) + cb) + dfII_dMy*Db12(i)) +  &
                    dfII_dMy*(dfII_dMx*Db21(i) + dfII_dMy*(Db22(i) + cb)) +  &
                    dfII_dMxy*(dfII_dMxy*(Db33(i) + cb))
              !< Right-hand side of the system
              b1 = f1p(i) 
              b2 = f2p(i)
              !< Compute inverse of A, and lambda_I, lambda_II
              det2  = A11*A22 - A12*A21
              if (abs(det2) < em20) then
                lam_I  = b1 / max(A11, em20)
                lam_II = zero
              else
                lam_I  = (b1*A22 - b2*A12) / det2
                lam_II = (A11*b2 - A21*b1) / det2
              endif
              !< Compute the iteration plastic strain increment 
              depsp_x = lam_I*dfI_dNx  + lam_II*dfII_dNx
              depsp_y = lam_I*dfI_dNy  + lam_II*dfII_dNy
              dkp_x   = lam_I*dfI_dMx  + lam_II*dfII_dMx
              dkp_y   = lam_I*dfI_dMy  + lam_II*dfII_dMy
              dkp_xy  = lam_I*dfI_dMxy + lam_II*dfII_dMxy
              !< Update integrated plastic strains for post-processing and output
              uvar(i,4) = uvar(i,4) + depsp_x
              uvar(i,5) = uvar(i,5) + depsp_y
              uvar(i,6) = uvar(i,6) + dkp_x
              uvar(i,7) = uvar(i,7) + dkp_y
              uvar(i,8) = uvar(i,8) + dkp_xy
              !< Update equivalent plastic strain for output and post-processing
              pla(i)    = pla(i)    + lam_I + lam_II
              pla(i)    = max(pla(i), zero)
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm12*depsp_x - Dm11*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y
              momnxy(i) = momnxy(i) - Db33(i)*dkp_xy
              !< Update of the kinematic hardening variables
              ! -> Remove the previous back-stress contribution
              signxx(i) = signxx(i) + sigb(i,1)
              signyy(i) = signyy(i) + sigb(i,2)
              momnxx(i) = momnxx(i) + sigb(i,3)
              momnyy(i) = momnyy(i) + sigb(i,4)
              momnxy(i) = momnxy(i) + sigb(i,5)
              ! -> Update the back-stress variables with the new plastic strain increment
              sigb(i,1) = sigb(i,1) + cm*depsp_x
              sigb(i,2) = sigb(i,2) + cm*depsp_y
              sigb(i,3) = sigb(i,3) + cb*dkp_x
              sigb(i,4) = sigb(i,4) + cb*dkp_y
              sigb(i,5) = sigb(i,5) + cb*dkp_xy
              ! -> Add the new back-stress contribution
              signxx(i) = signxx(i) - sigb(i,1)
              signyy(i) = signyy(i) - sigb(i,2)
              momnxx(i) = momnxx(i) - sigb(i,3)
              momnyy(i) = momnyy(i) - sigb(i,4)
              momnxy(i) = momnxy(i) - sigb(i,5)
              !< Update the limit moments for the plasticity criteria
              !< - Positive bending moment for the reinforcement in x direction
              call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,zero,xi_pos,mfx_pos(i),dmfx_pos(i))
              !< Negative bending moment
              call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,zero,xi_neg,mfx_neg(i),dmfx_neg(i))
              !< Positive bending moment
              call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,pi*half,xi_pos,mfy_pos(i),dmfy_pos(i))
              !< Negative bending moment
              call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,pi*half,xi_neg,mfy_neg(i),dmfy_neg(i))      
              !< Update the plasticity criteria
              f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +    &
                                                      momnxy(i)*momnxy(i)
              f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +    &
                                                      momnxy(i)*momnxy(i)
              g1p(i) =   momnxx(i) - mfx_pos(i) + momnyy(i) - mfy_pos(i)
              g2p(i) = -(momnxx(i) - mfx_neg(i) + momnyy(i) - mfy_neg(i))
              if (g1p(i) <= zero) exit
              if (g2p(i) <= zero) exit
            enddo
          enddo
        endif
!
        !-----------------------------------------------------------------------
        !< Only criterion I active
        !-----------------------------------------------------------------------
        if (nindx2 > 0) then
          do ii = 1, nindx2
            i = indx2(ii)
            iter = 0
            do while ((abs(f1p(i)) > tol) .and. (iter < nmax))
              iter = iter + 1
              !< Normal derivatives of the criteria f_I 
              dfI_dMx  = -(momnyy(i) - mfy_pos(i))
              dfI_dMy  = -(momnxx(i) - mfx_pos(i))
              dfI_dMxy =  two * momnxy(i)
              dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dmfx_pos(i)
              dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dmfy_pos(i)
              !< Dénominateur de lambda_p
              !< terme membrane : (df/dN)^T (Dm + Cm) (df/dN)
              den = dfI_dNx*((Dm11 + cm)*dfI_dNx + Dm12*dfI_dNy) +           & 
                    dfI_dNy*(Dm12*dfI_dNx + (Dm11 + cm)*dfI_dNy) +           &
              !< terme flexion : (df/dM)^T (Db + Cb) (df/dM)
                    dfI_dMx*((Db11(i) + cb)*dfI_dMx + Db12(i)*dfI_dMy) +     &
                    dfI_dMy*(Db21(i)*dfI_dMx + (Db22(i) + cb)*dfI_dMy) +     &
                    dfI_dMxy*((Db33(i) + cb)*dfI_dMxy)
              !< Plastic multiplier
              den = sign(max(abs(den),em20),den)      
              lam_I = f1p(i) / den
              !< Correction des déformations plastiques
              depsp_x   = lam_I * dfI_dNx
              depsp_y   = lam_I * dfI_dNy
              dkp_x     = lam_I * dfI_dMx
              dkp_y     = lam_I * dfI_dMy
              dkp_xy    = lam_I * dfI_dMxy
              !< Update integrated plastic strains for post-processing and output
              uvar(i,4) = uvar(i,4) + depsp_x
              uvar(i,5) = uvar(i,5) + depsp_y
              uvar(i,6) = uvar(i,6) + dkp_x
              uvar(i,7) = uvar(i,7) + dkp_y
              uvar(i,8) = uvar(i,8) + dkp_xy
              !< Update equivalent plastic strain for output and post-processing
              pla(i)    = pla(i)    + lam_I
              pla(i)    = max(pla(i), zero)
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm12*depsp_x - Dm11*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y
              momnxy(i) = momnxy(i) - Db33(i)*dkp_xy
              !< Update kinematic hardening variables
              ! -> Remove previous contribution of the back-stress to the trial stress
              signxx(i) = signxx(i) + sigb(i,1)
              signyy(i) = signyy(i) + sigb(i,2)
              momnxx(i) = momnxx(i) + sigb(i,3)
              momnyy(i) = momnyy(i) + sigb(i,4)
              momnxy(i) = momnxy(i) + sigb(i,5)
              ! -> Update of the back-stress with the new plastic strain increment
              sigb(i,1) = sigb(i,1) + cm*depsp_x
              sigb(i,2) = sigb(i,2) + cm*depsp_y
              sigb(i,3) = sigb(i,3) + cb*dkp_x
              sigb(i,4) = sigb(i,4) + cb*dkp_y
              sigb(i,5) = sigb(i,5) + cb*dkp_xy
              ! -> Re-apply the back-stress contribution to the updated trial stress
              signxx(i) = signxx(i) - sigb(i,1)
              signyy(i) = signyy(i) - sigb(i,2)
              momnxx(i) = momnxx(i) - sigb(i,3)
              momnyy(i) = momnyy(i) - sigb(i,4)
              momnxy(i) = momnxy(i) - sigb(i,5)
              !< Update the limit moments for the plasticity criteria
              !< - Positive bending moment for the reinforcement in x direction
              call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,zero,xi_pos,mfx_pos(i),dmfx_pos(i))
              !< - Positive bending moment for the reinforcement in y direction
              call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,pi*half,xi_pos,mfy_pos(i),dmfy_pos(i))
              !< Update the plasticity criteria
              f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +    &
                                                      momnxy(i)*momnxy(i)
              g1p(i) =   momnxx(i) - mfx_pos(i) + momnyy(i) - mfy_pos(i)
              if (g1p(i) <= zero) exit
            enddo
          enddo
        endif
!
        !-----------------------------------------------------------------------
        !< Only criterion II active
        !-----------------------------------------------------------------------
        if (nindx3 > 0) then
          do ii = 1, nindx3
            i = indx3(ii)
            iter = 0
            do while ((abs(f2p(i)) > tol) .and. (iter < nmax))
              iter = iter + 1
              !< Normal derivatives of the criteria f_II
              dfII_dMx  = -(momnyy(i) - mfy_neg(i))
              dfII_dMy  = -(momnxx(i) - mfx_neg(i))
              dfII_dMxy =  two * momnxy(i)
              dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dmfx_neg(i)
              dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dmfy_neg(i)
              !< Dénominateur de lambda_p
              !< terme membrane : (df/dN)^T (Dm + Cm) (df/dN)
              den = dfII_dNx*((Dm11 + cm)*dfII_dNx + Dm12*dfII_dNy) +        &
                    dfII_dNy*(Dm12*dfII_dNx + (Dm11 + cm)*dfII_dNy) +        &
              !< terme flexion : (df/dM)^T (Db + Cb) (df/dM)
                    dfII_dMx*((Db11(i) + cb)*dfII_dMx + Db12(i)*dfII_dMy) +  &
                    dfII_dMy*(Db21(i)*dfII_dMx + (Db22(i) + cb)*dfII_dMy) +  &
                    dfII_dMxy*((Db33(i) + cb)*dfII_dMxy)
              !< Plastic multiplier
              den = sign(max(abs(den),em20),den)      
              lam_II = f2p(i) / den      
              !< Correction des déformations plastiques 
              depsp_x   = lam_II * dfII_dNx
              depsp_y   = lam_II * dfII_dNy
              dkp_x     = lam_II * dfII_dMx
              dkp_y     = lam_II * dfII_dMy
              dkp_xy    = lam_II * dfII_dMxy
              !< Update integrated plastic strains for post-processing and output
              uvar(i,4) = uvar(i,4) + depsp_x
              uvar(i,5) = uvar(i,5) + depsp_y
              uvar(i,6) = uvar(i,6) + dkp_x
              uvar(i,7) = uvar(i,7) + dkp_y
              uvar(i,8) = uvar(i,8) + dkp_xy
              !< Update equivalent plastic strain for output and post-processing
              pla(i)    = pla(i)    + lam_II
              pla(i)    = max(pla(i), zero)
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm12*depsp_x - Dm11*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x - Db12(i)*dkp_y 
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x - Db22(i)*dkp_y
              momnxy(i) = momnxy(i) - Db33(i)*dkp_xy
              !< Update kinematic hardening variables
              ! -> Remove previous contribution of the back-stress to the trial stress
              signxx(i) = signxx(i) + sigb(i,1)
              signyy(i) = signyy(i) + sigb(i,2)
              momnxx(i) = momnxx(i) + sigb(i,3)
              momnyy(i) = momnyy(i) + sigb(i,4)
              momnxy(i) = momnxy(i) + sigb(i,5)
              ! -> Update of the back-stress with the new plastic strain increment
              sigb(i,1) = sigb(i,1) + cm*depsp_x
              sigb(i,2) = sigb(i,2) + cm*depsp_y
              sigb(i,3) = sigb(i,3) + cb*dkp_x
              sigb(i,4) = sigb(i,4) + cb*dkp_y
              sigb(i,5) = sigb(i,5) + cb*dkp_xy
              ! -> Re-apply the back-stress contribution to the updated trial stress
              signxx(i) = signxx(i) - sigb(i,1)
              signyy(i) = signyy(i) - sigb(i,2)
              momnxx(i) = momnxx(i) - sigb(i,3)
              momnyy(i) = momnyy(i) - sigb(i,4)
              momnxy(i) = momnxy(i) - sigb(i,5)
              !< Update the limit moments for the plasticity criteria
              !< - Negative bending moment for the reinforcement in x direction
              call calc_M(signxx(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,zero,xi_neg,mfx_neg(i),dmfx_neg(i))
              !< - Negative bending moment for the reinforcement in y direction
              call calc_M(signyy(i),thk0(i),f_c,sig_y,omega_x,omega_y,rho_x,   &
                rho_y,pi*half,xi_neg,mfy_neg(i),dmfy_neg(i))                        
              !< Update the plasticity criteria
              f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +    &
                                                      momnxy(i)*momnxy(i)
              g2p(i) = -(momnxx(i) - mfx_neg(i) + momnyy(i) - mfy_neg(i))
              if (g2p(i) <= zero) exit
            enddo
          enddo
        endif
! 
        !< Remove the contribution of the back-stress to the stresses 
        do i = 1, nel
          signxx(i) = signxx(i) + sigb(i,1)
          signyy(i) = signyy(i) + sigb(i,2)
          momnxx(i) = momnxx(i) + sigb(i,3)
          momnyy(i) = momnyy(i) + sigb(i,4)
          momnxy(i) = momnxy(i) + sigb(i,5)
        enddo
!
        !=======================================================================
        end subroutine sigeps136g
!
        !=======================================================================
        !< Function to compute the value of N_pos for a given eta
        !=======================================================================
        function calc_N(                                                       &
            eta      ,thk0     ,f_c      ,sig_y    ,omega_x  ,omega_y  ,       &
            rho_x    ,rho_y    ,alpha    ,xi       )                           & 
            result(N_val)
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
          use precision_mod, only: WP
          use constant_mod
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------          
          implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
          real(kind=WP), intent(in) :: eta        !< Position of the neutral axis (normalized by the thickness)
          real(kind=WP), intent(in) :: thk0       !< Thickness of the section
          real(kind=WP), intent(in) :: f_c        !< Concrete compressive strength
          real(kind=WP), intent(in) :: sig_y(2)   !< Yield stress of the reinforcement 
          real(kind=WP), intent(in) :: omega_x(2) !< Area of the reinforcement in x direction
          real(kind=WP), intent(in) :: omega_y(2) !< Area of the reinforcement in y direction
          real(kind=WP), intent(in) :: rho_x(2)   !< Position of the reinforcement in x direction (normalized by the thickness)
          real(kind=WP), intent(in) :: rho_y(2)   !< Position of the reinforcement in y direction (normalized by the thickness)
          real(kind=WP), intent(in) :: alpha      !< Angle of the reinforcement
          real(kind=WP), intent(in) :: xi         !< +1 for positive bending, -1 for negative bending
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
          real(kind=WP)             :: N_val      !< Value of N for the given eta
          real(kind=WP)             :: xi_inf_x   !< Sign for the reinforcement below the neutral axis in x direction
          real(kind=WP)             :: xi_sup_x   !< Sign for the reinforcement above the neutral axis in x direction
          real(kind=WP)             :: xi_inf_y   !< Sign for the reinforcement below the neutral axis in y direction
          real(kind=WP)             :: xi_sup_y   !< Sign for the reinforcement above the neutral axis in y direction
          real(kind=WP)             :: gamma_x    !< Contribution of the reinforcement in x direction to the bending resistance
          real(kind=WP)             :: gamma_y    !< Contribution of the reinforcement in y direction to the bending resistance
!
          !< Computation of the sign of the reinforcement contribution depending 
          !  on the position of the neutral axis
          if (abs(eta - rho_x(1)) < em04) then
            xi_inf_x = zero
          else
            xi_inf_x = sign(one, eta - rho_x(1))
          endif
          if (abs(eta - rho_x(2)) < em04) then
            xi_sup_x = zero
          else
            xi_sup_x = sign(one, eta - rho_x(2))
          endif
          if (abs(eta - rho_y(1)) < em04) then
            xi_inf_y = zero
          else
            xi_inf_y = sign(one, eta - rho_y(1))
          endif
          if (abs(eta - rho_y(2)) < em04) then
            xi_sup_y = zero
          else
            xi_sup_y = sign(one, eta - rho_y(2))
          endif     
!
          !< Computation of the contribution of the reinforcement to the bending resistance
          gamma_x = (xi_sup_x*omega_x(2)*sig_y(2) +                            &
                     xi_inf_x*omega_x(1)*sig_y(1)) / (f_c*thk0)
          gamma_y = (xi_sup_y*omega_y(2)*sig_y(2) +                            &
                     xi_inf_y*omega_y(1)*sig_y(1)) / (f_c*thk0)                                                               
!
          !< Computation of N for the given eta
          N_val = xi*gamma_x*(cos(alpha))**2 + xi*gamma_y*(sin(alpha))**2 -    &
                                                   (one - xi*eta)*half
          N_val = N_val * f_c * thk0
        end function calc_N
!
        !=======================================================================
        !< Computation of pos. bending moment M_pos and its derivative w.r.t N
        !=======================================================================
        subroutine calc_M(                                                     &
          sigma    ,thk0     ,f_c      ,sig_y    ,omega_x  ,omega_y  ,rho_x   ,&
          rho_y    ,alpha    ,xi       ,M_val    ,dM_dN    )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
          use precision_mod, only: WP
          use constant_mod
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
          implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
          real(kind=WP), intent(in) :: sigma       !< Membrane stress in the direction of the bending moment
          real(kind=WP), intent(in) :: thk0        !< Thickness of the section
          real(kind=WP), intent(in) :: f_c         !< Concrete compressive strength
          real(kind=WP), intent(in) :: sig_y(2)    !< Yield stress of the reinforcement
          real(kind=WP), intent(in) :: omega_x(2)  !< Area of the reinforcement in x direction
          real(kind=WP), intent(in) :: omega_y(2)  !< Area of the reinforcement in y direction
          real(kind=WP), intent(in) :: rho_x(2)    !< Position of the reinforcement in x direction (normalized by the thickness)
          real(kind=WP), intent(in) :: rho_y(2)    !< Position of the reinforcement in y direction (normalized by the thickness)
          real(kind=WP), intent(in) :: alpha       !< Angle of the reinforcement
          real(kind=WP), intent(in) :: xi          !< +1 for positive bending, -1 for negative bending
          real(kind=WP), intent(out) :: M_val      !< Value of the bending moment for the given sigma
          real(kind=WP), intent(out) :: dM_dN      !< Derivative of the bending moment w.r.t. N for the given sigma
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
          real(kind=WP) :: eta_lo, eta_hi, eta_mid
          real(kind=WP) :: N_lo, N_hi, N_mid, N
          real(kind=WP) :: n_f, xi_inf_x, xi_sup_x, xi_inf_y, xi_sup_y
          real(kind=WP) :: gamma_x, gamma_y, rhox, rhoy, denom_x, denom_y
          integer :: iter
          logical :: found
          real(kind=WP), parameter :: tol = 1.0d-8   ! tolerance sur eta
          real(kind=WP), parameter :: eta_min = -1.0d0
          real(kind=WP), parameter :: eta_max = 1.0d0
          integer,       parameter :: max_iter = 100
!
          n_f = zero
          N_lo = calc_N(eta_min,thk0,f_c,sig_y,omega_x,omega_y,rho_x,rho_y,alpha,xi)
          N_hi = calc_N(eta_max,thk0,f_c,sig_y,omega_x,omega_y,rho_x,rho_y,alpha,xi)
          N = sigma*thk0
          found = .false.
          eta_mid = half * (eta_min + eta_max)
!
          !< Check if N is within the range of N(eta) for eta in [eta_min, eta_max]
          if (N*xi < N_lo*xi .or. N*xi > N_hi*xi) then
            if (N*xi < N_lo*xi) then
              n_f = eta_min
            else
              n_f = eta_max
            end if
            found = .true.
          endif
!
          !< Bisection to find the value of eta such that N(eta) = sigma*thk0
          eta_lo = eta_min
          eta_hi = eta_max
          if (.not.found) then
            do iter = 1, max_iter
              eta_mid = half * (eta_lo + eta_hi)
              N_mid   = calc_N(eta_mid,thk0,f_c,sig_y,omega_x,omega_y,rho_x,rho_y,alpha,xi)
              !< Check convergence
              if (abs(N_mid - N) < tol * abs(N) + tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
              !< Check if the interval is too small
              if ((eta_hi - eta_lo) < tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
              !< Narrow the interval
              if (N_mid*xi < N*xi) then
                eta_lo = eta_mid
              else
                eta_hi = eta_mid
              endif
            enddo
          endif
!
          ! Convergence atteinte apres max_iter (ne devrait pas arriver)
          if (.not.found) then 
            n_f = eta_mid
            found = .true.
          endif     
!
          !< Computation of the sign of the reinforcement contribution depending 
          !  on the position of the neutral axis
          if (abs(n_f - rho_x(1)) < em04) then
            xi_inf_x = zero
          else
            xi_inf_x = sign(one, n_f - rho_x(1))
          endif
          if (abs(n_f - rho_x(2)) < em04) then
            xi_sup_x = zero
          else
            xi_sup_x = sign(one, n_f - rho_x(2))
          endif
          if (abs(n_f - rho_y(1)) < em04) then
            xi_inf_y = zero
          else
            xi_inf_y = sign(one, n_f - rho_y(1))
          endif
          if (abs(n_f - rho_y(2)) < em04) then
            xi_sup_y = zero
          else            
            xi_sup_y = sign(one, n_f - rho_y(2))
          endif
!
          !< Computation of the contribution of the reinforcement to the bending resistance
          gamma_x = (xi_sup_x*omega_x(2)*sig_y(2) +                            &
                     xi_inf_x*omega_x(1)*sig_y(1)) / (f_c*thk0)
          gamma_y = (xi_sup_y*omega_y(2)*sig_y(2) +                            &
                     xi_inf_y*omega_y(1)*sig_y(1)) / (f_c*thk0)    
!
          !< Computation of the position of the equivalent reinforcement for the concrete contribution
          denom_x = xi_sup_x*omega_x(2) + xi_inf_x*omega_x(1)
          if (abs(denom_x) < em20) then
            rhox = zero
          else
            rhox = (xi_sup_x*omega_x(2)*rho_x(2) + xi_inf_x*omega_x(1)*rho_x(1)) / &
                   (xi_sup_x*omega_x(2) + xi_inf_x*omega_x(1))
          endif
          denom_y = xi_sup_y*omega_y(2) + xi_inf_y*omega_y(1)
          if (abs(denom_y) < em20) then
            rhoy = zero
          else
            rhoy = (xi_sup_y*omega_y(2)*rho_y(2) + xi_inf_y*omega_y(1)*rho_y(1)) / &
                 (xi_sup_y*omega_y(2) + xi_inf_y*omega_y(1))
          endif
!
          !< Computation of the bending moment for the given sigma
          M_val = -half*xi*rhox*gamma_x*(cos(alpha))**2 -                      &
                   half*xi*rhoy*gamma_y*(sin(alpha))**2 +                      &
                   xi*fourth*(one - n_f**2)
          M_val = M_val * f_c * thk0**2
!
          !< Computation of the derivative of the bending moment w.r.t. N for the given sigma
          dM_dN = - thk0 * n_f
!
        end subroutine calc_M
!
      end module sigeps136g_mod
