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
        integer :: i, iter, k, nlayer, ii, j, nindx, nindx2, nindx3,           &
          indx(nel), indx2(nel), indx3(nel), ierr
        real(kind=WP) :: young,nu,shear,lambda_m,mu_m,f_t,f_c,gamma,           &
          dmax,Dm11,Dm12,Dm21,Dm22,cn1x,cn1y,cn1xy,cn2x,cn2y,cn2xy,            &
          cm1x,cm1y,cm1xy,cm2x,cm2y,cm2xy
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
          mfy_neg, dmfx_pos, dmfx_neg, dmfy_pos, dmfy_neg, f1p, f2p, Db32,     &
          Db31, Db13, Db23
        real(kind=WP), dimension(:), allocatable ::rho_x, rho_y, sig_y,        &
          omega_x, omega_y
        real(kind=WP), parameter :: tol = 1e-8
!
        !=======================================================================
        !< - Initialisation of computation on time step
        !=======================================================================
!
        !< Recover integer model parameters
        nlayer   = matparam%iparam(1)                      !< Number of reinforcement layers
        allocate(sig_y(nlayer),omega_x(nlayer),                                &
          omega_y(nlayer),rho_x(nlayer),rho_y(nlayer))
        !< Recovering real model parameters
        young    = matparam%young                          !< Concrete Young's modulus
        nu       = matparam%nu                             !< Concrete Poisson's ratio
        shear    = matparam%shear                          !< Concrete shear modulus
        lambda_m = matparam%uparam(1)                      !< Membrane Lamé parameter
        mu_m     = matparam%uparam(2)                      !< Membrane shear modulus
        f_t      = matparam%uparam(5)                      !< Concrete Tensile strength
        f_c      = matparam%uparam(6)                      !< Concrete Compressive strength
        gamma    = matparam%uparam(7)                      !< Concrete 
        dmax     = matparam%uparam(8)                      !< Maximum damage
        cn1x     = matparam%uparam(9) 
        cn1y     = matparam%uparam(10)                     
        cn1xy    = matparam%uparam(11)                     
        cn2x     = matparam%uparam(12)
        cn2y     = matparam%uparam(13)                     
        cn2xy    = matparam%uparam(14)                     
        cm1x     = matparam%uparam(15)
        cm1y     = matparam%uparam(16)                     
        cm1xy    = matparam%uparam(17)                     
        cm2x     = matparam%uparam(18)
        cm2y     = matparam%uparam(19)                    
        cm2xy    = matparam%uparam(20)                     
        do k = 1,nlayer
          sig_y(k)   = matparam%uparam(20 + 5*(k-1) + 1)   !< Yield stress of the steel reinforcement in layer k
          omega_x(k) = matparam%uparam(20 + 5*(k-1) + 2)   !< Area of the steel reinforcement in x direction in layer k per unit length of concrete
          omega_y(k) = matparam%uparam(20 + 5*(k-1) + 3)   !< Area of the steel reinforcement in y direction in layer k per unit length of concrete
          rho_x(k)   = matparam%uparam(20 + 5*(k-1) + 4)   !< Through thickness position of the steel reinforcement in x direction in layer k (normalized by the initial thickness)
          rho_y(k)   = matparam%uparam(20 + 5*(k-1) + 5)   !< Through thickness position of the steel reinforcement in y direction in layer k (normalized by the initial thickness)
        enddo
!
        !< Computation of some real parameters
        Dm11 = lambda_m + two*mu_m                         !< Membrane stiffness D11 = lambda + 2*mu
        Dm12 = lambda_m                                    !< Membrane stiffness D12 = lambda
        Dm21 = Dm12                                        !< Membrane stiffness D21 = lambda
        Dm22 = Dm11                                        !< Membrane stiffness D22 = lambda + 2*mu 
        gs(1:nel) = shear*shf(1:nel)                       !< Correction factor for transverse shear 
        lambda_b(1:nel) = thk0(1:nel)*matparam%uparam(3)   !< Bending Lamé parameter
        mu_b(1:nel) = thk0(1:nel)*matparam%uparam(4)       !< Bending shear modulus
        k0_1(1:nel) = f_t*f_t*(one - nu*nu) / &            !< Bending damage threshold in positive bending (tension inner face)                  
                      (six*young*thk0(1:nel)*thk0(1:nel)) 
        k0_2(1:nel) = f_c*f_c*(one - nu*nu) / &            !< Bending damage threshold in negative bending (tension outer face)
                      (six*young*thk0(1:nel)*thk0(1:nel)) 
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
        !=======================================================================
        !< - COMPUTATION OF THE MEMBRANE TRIAL STRESS TENSOR
        !=======================================================================
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
          Db13(i) = zero
          Db21(i) = Db12(i)
          Db22(i) = Db22(i) + fac1*dW1_dkyy*dW1_dkyy + fac2*dW2_dkyy*dW2_dkyy
          Db23(i) = zero
          Db31(i) = zero
          Db32(i) = zero
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
          call calc_M_pos(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,      &
            mfx_pos(i), dmfx_pos(i))
          !< Negative bending moment
          call calc_M_neg(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,      &
            mfx_neg(i), dmfx_neg(i)) 
        enddo
!
        !< Critical moment for the plasticity criterion in bending for the 
        !  reinforcement in y direction
        !-----------------------------------------------------------------------        
        do i = 1,nel
          !< Positive bending moment
          call calc_M_pos(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,      &
            mfy_pos(i), dmfy_pos(i))
          !< Negative bending moment
          call calc_M_neg(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,      &
            mfy_neg(i), dmfy_neg(i))
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
        enddo
!
        !< Index of the active criteria
        nindx  = 0
        nindx2 = 0
        nindx3 = 0
        do i = 1, nel
          if (f1p(i) > zero .and. f2p(i) > zero) then
            nindx = nindx + 1
            indx(nindx) = i
          else if (f1p(i) > zero) then
            nindx2 = nindx2 + 1
            indx2(nindx2) = i
          else if (f2p(i) > zero) then
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
            i = indx(ii)
            do while (abs(f1p(i)) > tol .or. abs(f2p(i)) > tol)
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
              A11 = dfI_dNx**2*Dm11 + two*dfI_dNx*dfI_dNy*Dm12                 &
                  + dfI_dNy**2*Dm22                                            &
                  + dfI_dMx**2*Db11(i) + two*dfI_dMx*dfI_dMy*Db12(i)           &
                  + dfI_dMy**2*Db22(i) + dfI_dMxy**2*Db33(i)
              A12 = dfI_dNx*dfII_dNx*Dm11                                      &
                  + (dfI_dNx*dfII_dNy+dfI_dNy*dfII_dNx)*Dm12                   &
                  + dfI_dNy*dfII_dNy*Dm22                                      &
                  + dfI_dMx*dfII_dMx*Db11(i)                                   &
                  + (dfI_dMx*dfII_dMy+dfI_dMy*dfII_dMx)*Db12(i)                &
                  + dfI_dMy*dfII_dMy*Db22(i)                                   &
                  + dfI_dMxy*dfII_dMxy*Db33(i)
              A21 = A12
              A22 = dfII_dNx**2*Dm11 + two*dfII_dNx*dfII_dNy*Dm12              &
                  + dfII_dNy**2*Dm22                                           &
                  + dfII_dMx**2*Db11(i) + two*dfII_dMx*dfII_dMy*Db12(i)        &
                  + dfII_dMy**2*Db22(i) + dfII_dMxy**2*Db33(i)
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
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm21*depsp_x - Dm22*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x -                          &
                                      Db12(i)*dkp_y -                          &
                                      Db13(i)*dkp_xy
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x -                          &
                                      Db22(i)*dkp_y -                          &
                                      Db23(i)*dkp_xy
              momnxy(i) = momnxy(i) - Db31(i)*dkp_x -                          &
                                      Db32(i)*dkp_y -                          &
                                      Db33(i)*dkp_xy
              !< Update the limit moments for the plasticity criteria
              !< - Positive bending moment for the reinforcement in x direction
              call calc_M_pos(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,  &
                mfx_pos(i), dmfx_pos(i))
              !< - Negative bending moment for the reinforcement in x direction
              call calc_M_neg(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,  &
                mfx_neg(i), dmfx_neg(i))
              !< - Positive bending moment for the reinforcement in y direction
              call calc_M_pos(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,  &
                mfy_pos(i), dmfy_pos(i))
              !< - Negative bending moment for the reinforcement in y direction
              call calc_M_neg(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,  &
                mfy_neg(i), dmfy_neg(i))         
              !< Update the plasticity criteria
              f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +    &
                                                      momnxy(i)*momnxy(i)
              f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +    &
                                                      momnxy(i)*momnxy(i)
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
            do while (abs(f1p(i)) > tol)
              !< Normal derivatives of the criteria f_I 
              dfI_dMx  = -(momnyy(i) - mfy_pos(i))
              dfI_dMy  = -(momnxx(i) - mfx_pos(i))
              dfI_dMxy =  two * momnxy(i)
              dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dmfx_pos(i)
              dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dmfy_pos(i)
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
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm21*depsp_x - Dm22*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x -                          &
                                      Db12(i)*dkp_y -                          &
                                      Db13(i)*dkp_xy
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x -                          &
                                      Db22(i)*dkp_y -                          &
                                      Db23(i)*dkp_xy
              momnxy(i) = momnxy(i) - Db31(i)*dkp_x -                          &
                                      Db32(i)*dkp_y -                          &
                                      Db33(i)*dkp_xy
              !< Update the limit moments for the plasticity criteria
              !< - Positive bending moment for the reinforcement in x direction
              call calc_M_pos(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,  &
                mfx_pos(i), dmfx_pos(i))
              !< - Positive bending moment for the reinforcement in y direction
              call calc_M_pos(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,  &
                mfy_pos(i), dmfy_pos(i))
              !< Update the plasticity criteria
              f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +    &
                                                      momnxy(i)*momnxy(i)
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
            do while (abs(f2p(i)) > tol)
              !< Normal derivatives of the criteria f_II
              dfII_dMx  = -(momnyy(i) - mfy_neg(i))
              dfII_dMy  = -(momnxx(i) - mfx_neg(i))
              dfII_dMxy =  two * momnxy(i)
              dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dmfx_neg(i)
              dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dmfy_neg(i)
              !< Dénominateur de lambda_p
              !< terme membrane : (df/dN)^T (Dm + Cm) (df/dN)
              den = dfII_dNx**2 * Dm11                                         &
                  + two*dfII_dNx*dfII_dNy*Dm12                                 &
                  + dfII_dNy**2 * Dm22                                         &
              !< terme flexion : (df/dM)^T (Db + Cb) (df/dM)
                  + dfII_dMx**2  * Db11(i)                                     &
                  + two*dfII_dMx*dfII_dMy*Db12(i)                              &
                  + dfII_dMy**2  * Db22(i)                                     &
                  + dfII_dMxy**2 * Db33(i)
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
              !< Update of the membrane stresses and bending moments
              signxx(i) = signxx(i) - Dm11*depsp_x - Dm12*depsp_y
              signyy(i) = signyy(i) - Dm21*depsp_x - Dm22*depsp_y
              momnxx(i) = momnxx(i) - Db11(i)*dkp_x -                          &
                                      Db12(i)*dkp_y -                          &
                                      Db13(i)*dkp_xy
              momnyy(i) = momnyy(i) - Db21(i)*dkp_x -                          &
                                      Db22(i)*dkp_y -                          &
                                      Db23(i)*dkp_xy
              momnxy(i) = momnxy(i) - Db31(i)*dkp_x -                          &
                                      Db32(i)*dkp_y -                          &
                                      Db33(i)*dkp_xy
              !< Update the limit moments for the plasticity criteria
              !< - Negative bending moment for the reinforcement in x direction
              call calc_M_neg(signxx(i), thk0(i), f_c, sig_y, omega_x, rho_x,  &
                mfx_neg(i), dmfx_neg(i))
              !< - Negative bending moment for the reinforcement in y direction
              call calc_M_neg(signyy(i), thk0(i), f_c, sig_y, omega_y, rho_y,  &
                mfy_neg(i), dmfy_neg(i))                        
              !< Update the plasticity criteria
              f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +    &
                                                      momnxy(i)*momnxy(i)
            enddo
          enddo
        endif
!
        !< Deallocation of the extended reinforcement position tables
        deallocate(sig_y, omega_x, omega_y, rho_x, rho_y)
!
        !=======================================================================
        end subroutine sigeps136g
!
        !=======================================================================
        !< Function to compute the value of N_pos for a given eta
        !=======================================================================
        function N_pos(eta, coef_c, sig_y, omega_x, rho_x) result(N_val)
          implicit none
          real(kind=8), intent(in) :: eta, coef_c
          real(kind=8), intent(in) :: sig_y(2), omega_x(2), rho_x(2)
          real(kind=8)             :: N_val
          real(kind=8)             :: xi_inf, xi_sup
        
          ! sign continu : sign(1.d0, 0.d0) = +1 en Fortran => coherent avec thèse
          if (abs(eta - rho_x(1)) < 1.0d-4) then
            xi_inf = 0.0d0
          else
            xi_inf = sign(1.0d0, eta - rho_x(1))
          end if
          if (abs(eta - rho_x(2)) < 1.0d-4) then
            xi_sup = 0.0d0
          else
            xi_sup = sign(1.0d0, eta - rho_x(2))
          end if
        
          N_val = coef_c * (eta - 1.0d0)          &
                + sig_y(1) * omega_x(1) * xi_inf &
                + sig_y(2) * omega_x(2) * xi_sup
        end function N_pos
!
        !=======================================================================
        !< Computation of pos. bending moment M_pos and its derivative w.r.t N
        !=======================================================================
        subroutine calc_M_pos(sigma, thk0, f_c, sig_y, omega_x, rho_x, M_pos, dM_pos_dN)
          use precision_mod, only: WP
          use constant_mod
          implicit none
          real(kind=WP), intent(in) :: sigma, thk0, f_c
          real(kind=WP), intent(in) :: sig_y(2), omega_x(2), rho_x(2)
          real(kind=WP), intent(out) :: M_pos, dM_pos_dN
!
          real(kind=WP) :: eta_lo, eta_hi, eta_mid
          real(kind=WP) :: N_lo, N_hi, N_mid, N
          real(kind=WP) :: n_f, coef_c, xi_inf, xi_sup
          integer :: iter
          logical :: found
          real(kind=WP), parameter :: tol = 1.0d-8   ! tolerance sur eta
          real(kind=WP), parameter :: eta_min = -1.0d0
          real(kind=WP), parameter :: eta_max = 1.0d0
          integer,       parameter :: max_iter = 100
!
          n_f = zero
          coef_c = thk0*f_c*half
          N_lo = N_pos(eta_min, coef_c, sig_y, omega_x, rho_x)
          N_hi = N_pos(eta_max, coef_c, sig_y, omega_x, rho_x)
          N = sigma*thk0
          found = .false.
          eta_mid = half * (eta_min + eta_max)
          if (N < N_lo .or. N > N_hi) then
            ! N hors domaine : on sature a la borne la plus proche
            if (N < N_lo) then
              n_f = eta_min
            else
              n_f = eta_max
            end if
            found = .true.
          end if
!
          ! --- Bisection ---
          eta_lo = eta_min
          eta_hi = eta_max
          if (.not.found) then
            do iter = 1, max_iter
              eta_mid = half * (eta_lo + eta_hi)
              N_mid   = N_pos(eta_mid, coef_c, sig_y, omega_x, rho_x)
          
              if (abs(N_mid - N) < tol * abs(N) + tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
          
              if ((eta_hi - eta_lo) < tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
          
              ! Resserrer l'intervalle
              if (N_mid < N) then
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
            write(*,*) 'AVERTISSEMENT : bisection non convergee apres ', iter, ' iterations' 
          endif     
          M_pos = thk0*half*(-N*n_f + f_c*thk0*fourth*(n_f - one)**2           & 
                             + sig_y(1)*omega_x(1)*abs(n_f - rho_x(1))         &
                             + sig_y(2)*omega_x(2)*abs(n_f - rho_x(2)))
          dM_pos_dN = - half* thk0 * n_f
!
        end subroutine calc_M_pos  
!
        !=======================================================================
        !< Function to compute the value of N_neg for a given eta
        !=======================================================================
        function N_neg(eta, coef_c, sig_y, omega_x, rho_x) result(N_val)
          implicit none
          real(kind=8), intent(in) :: eta, coef_c
          real(kind=8), intent(in) :: sig_y(2), omega_x(2), rho_x(2)
          real(kind=8)             :: N_val
          real(kind=8)             :: xi_inf, xi_sup
        
          ! sign continu : sign(1.d0, 0.d0) = +1 en Fortran => coherent avec thèse
          if (abs(eta - rho_x(1)) < 1.0d-4) then
            xi_inf = 0.0d0
          else
            xi_inf = sign(1.0d0, eta - rho_x(1))
          end if
          if (abs(eta - rho_x(2)) < 1.0d-4) then
            xi_sup = 0.0d0
          else
            xi_sup = sign(1.0d0, eta - rho_x(2))
          end if
        
          N_val = - coef_c * (eta + 1.0d0)                                     &
                - sig_y(1) * omega_x(1) * xi_inf                               &
                - sig_y(2) * omega_x(2) * xi_sup
        end function N_neg
!        
        !=======================================================================
        !< Computation of neg. bending moment M_neg and its derivative w.r.t N
        !=======================================================================
        subroutine calc_M_neg(sigma, thk0, f_c, sig_y, omega_x, rho_x, M_neg, dM_neg_dN)
          use precision_mod, only: WP
          use constant_mod
          implicit none
          real(kind=WP), intent(in) :: sigma, thk0, f_c
          real(kind=WP), intent(in) :: sig_y(2), omega_x(2), rho_x(2)
          real(kind=WP), intent(out) :: M_neg, dM_neg_dN
!
          real(kind=WP) :: eta_lo, eta_hi, eta_mid
          real(kind=WP) :: N_lo, N_hi, N_mid, N
          real(kind=WP) :: n_f, coef_c, xi_inf, xi_sup
          integer :: iter
          logical :: found
          real(kind=WP), parameter :: tol = 1.0d-8   ! tolerance sur eta
          real(kind=WP), parameter :: eta_min = -1.0d0
          real(kind=WP), parameter :: eta_max = 1.0d0
          integer,       parameter :: max_iter = 100
!
          n_f = zero
          coef_c = thk0*f_c*half
          N_lo = N_neg(eta_min, coef_c, sig_y, omega_x, rho_x)
          N_hi = N_neg(eta_max, coef_c, sig_y, omega_x, rho_x)
          N = sigma*thk0
          found  = .false.
          eta_mid = half * (eta_min + eta_max)
          if (N > N_lo .or. N < N_hi) then
            ! N hors domaine : on sature a la borne la plus proche
            if (N > N_lo) then
              n_f = eta_min
            else
              n_f = eta_max
            endif
            found = .true.
          endif
!
          ! --- Bisection ---
          eta_lo = eta_min
          eta_hi = eta_max
          if (.not.found) then
            do iter = 1, max_iter
              eta_mid = half * (eta_lo + eta_hi)
              N_mid   = N_neg(eta_mid, coef_c, sig_y, omega_x, rho_x)
          
              if (abs(N_mid - N) < tol * abs(N) + tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
          
              if ((eta_hi - eta_lo) < tol) then
                n_f   = eta_mid
                found = .true.
                exit
              endif
          
              ! Resserrer l'intervalle
              if (N_mid > N) then
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
            write(*,*) 'AVERTISSEMENT : bisection non convergee apres ', iter, ' iterations' 
          endif     
!
          M_neg = -thk0*half*(N*n_f + f_c*thk0*fourth*(n_f + one)**2           & 
                             + sig_y(1)*omega_x(1)*abs(n_f - rho_x(1))         &
                             + sig_y(2)*omega_x(2)*abs(n_f - rho_x(2)))
          dM_neg_dN = - half* thk0 * n_f
!
        end subroutine calc_M_neg 
!
      end module sigeps136g_mod
