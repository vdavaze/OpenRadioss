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
          Dm11     = lambda_m + two*mu_m
          Dm12     = lambda_m
          Dm21     = lambda_m
          Dm22     = lambda_m + two*mu_m 
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
          !< - Computation of the membrane trial stress tensor
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
            !< Computation of the damaged bending stiffness matrix
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
            !< -----------------------------------------------------------------------
            !< D_b^(2) : damage-variable contribution to the tangent stiffness (eq. 6.22)
            !< -----------------------------------------------------------------------
            !< Gradient of W+ w.r.t. principal curvatures (repère principal)
            dW1_dk1 = zero
            dW1_dk2 = zero
            if (tr_kb  >= zero) then
              dW1_dk1 = dW1_dk1 + lambda_b(i) * tr_kb   ! dW+/dkappa1 : trace term
              dW1_dk2 = dW1_dk2 + lambda_b(i) * tr_kb   ! dW+/dkappa2 : trace term
            end if
            if (kappa1 >= zero) dW1_dk1 = dW1_dk1 + two*mu_b(i)*kappa1
            if (kappa2 >= zero) dW1_dk2 = dW1_dk2 + two*mu_b(i)*kappa2
            !< Gradient of W- w.r.t. principal curvatures (repère principal)
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
            !< Scalar prefactors (eq. 6.22) : (gamma-1) / (2*(1+D)*W)
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
          !< - Computation of the plasticity
          !=====================================================================
          !< Loop over elements
          do i = 1, nel
            f1p(i) = -(momnxx(i) - mfx_pos(i))*(momnyy(i) - mfy_pos(i)) +      &
                                                   momnxy(i)*momnxy(i)
            f2p(i) = -(momnxx(i) - mfx_neg(i))*(momnyy(i) - mfy_neg(i)) +      &
                                                   momnxy(i)*momnxy(i)
            active_I  = (f1p(i)  > zero)
            active_II = (f2p(i)  > zero)

            !< =================================================================
            !< CAS 1 : Un seul critère actif (éq. 6.26 scalaire)
            !< =================================================================
            if (active_I .and. .not. active_II) then

              !< Gradients df_I/dM et df_I/dN (cutting plane : évalués en début)
              dfI_dMx  = -(momnyy(i) - mfy_pos(i))
              dfI_dMy  = -(momnxx(i) - mfx_pos(i))
              dfI_dMxy =  two * momnxy(i)
              dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dMfx_pos(i)
              dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dMfy_pos(i)
  
              num = fI                                                         &
                  + dfI_dNx  * (Dm11*depsxx(i) + Dm12*depsyy(i))               &
                  + dfI_dNy  * (Dm21*depsxx(i) + Dm22*depsyy(i))               &
                  + dfI_dMx  * (Db11(i)*depbxx(i) + Db12(i)*depbyy(i) +        &
                                                    Db13(i)*depbxy(i))         &
                  + dfI_dMy  * (Db21(i)*depbxx(i) + Db22(i)*depbyy(i) +        &
                                                    Db23(i)*depbxy(i))         &
                  + dfI_dMxy * (Db31(i)*depbxx(i) + Db32(i)*depbyy(i) +        &
                                                    Db33(i)*depbxy(i)) 
                  
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
            endif

            if (active_I .and. active_II) then
            
              !< Gradients f_I
              dfI_dMx  = -(momnyy(i) - mfy_pos(i))
              dfI_dMy  = -(momnxx(i) - mfx_pos(i))
              dfI_dMxy =  two * momnxy(i)
              dfI_dNx  =  (momnyy(i) - mfy_pos(i)) * dMfx_pos(i)
              dfI_dNy  =  (momnxx(i) - mfx_pos(i)) * dMfy_pos(i)
            
              !< Gradients f_II
              dfII_dMx  = -(momnyy(i) - mfy_neg(i))
              dfII_dMy  = -(momnxx(i) - mfx_neg(i))
              dfII_dMxy =  two * momnxy(i)
              dfII_dNx  =  (momnyy(i) - mfy_neg(i)) * dMfx_neg(i)
              dfII_dNy  =  (momnxx(i) - mfx_neg(i)) * dMfy_neg(i)
      
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
            
              !< Seconds membres
              b1 = fI  + dfI_dNx*(Dm11*depsxx(i)+Dm12*depsyy(i))               &
                       + dfI_dNy*(Dm21*depsxx(i)+Dm22*depsyy(i))               &
                       + dfI_dMx*(Db11(i)*depbxx(i)+Db12(i)*depbyy(i)+         &
                                                    Db13(i)*depbxy(i))         &
                       + dfI_dMy*(Db21(i)*depbxx(i)+Db22(i)*depbyy(i)+         &
                                                    Db23(i)*depbxy(i))         &
                       + dfI_dMxy*(Db31(i)*depbxx(i)+Db32(i)*depbyy(i)+        &
                                                     Db33(i)*depbxy(i))
            
              b2 = fII + dfII_dNx*(Dm11*depsxx(i)+Dm12*depsyy(i))              &
                       + dfII_dNy*(Dm21*depsxx(i)+Dm22*depsyy(i))              &
                       + dfII_dMx*(Db11(i)*depbxx(i)+Db12(i)*depbyy(i)+        &
                                                     Db13(i)*depbxy(i))        &
                       + dfII_dMy*(Db21(i)*depbxx(i)+Db22(i)*depbyy(i)+        &
                                                     Db23(i)*depbxy(i))        &
                       + dfII_dMxy*(Db31(i)*depbxx(i)+Db32(i)*depbyy(i)+       &
                                                      Db33(i)*depbxy(i))
            
              !< Résolution 2x2 par Cramer
              det2 = A11*A22 - A12*A21
              if (abs(det2) < em20) goto 999
              lam_I  = max(zero, (b1*A22 - b2*A12) / det2)
              lam_II = max(zero, (A11*b2 - A21*b1) / det2)
            
              !< Incréments plastiques combinés
              deps0p_x = lam_I*dfI_dNx  + lam_II*dfII_dNx
              deps0p_y = lam_I*dfI_dNy  + lam_II*dfII_dNy
              dkp_x    = lam_I*dfI_dMx  + lam_II*dfII_dMx
              dkp_y    = lam_I*dfI_dMy  + lam_II*dfII_dMy
              dkp_xy   = lam_I*dfI_dMxy + lam_II*dfII_dMxy
            
            end if

            !< =====================================================================
            !< Mise à jour des variables internes et des efforts
            !< =====================================================================
            
            !< Déformations plastiques cumulées
            ! eps0p_xx(i) = eps0p_xx(i) + deps0p_x
            ! eps0p_yy(i) = eps0p_yy(i) + deps0p_y
            ! kap_p_xx(i) = kap_p_xx(i) + dkp_x
            ! kap_p_yy(i) = kap_p_yy(i) + dkp_y
            ! kap_p_xy(i) = kap_p_xy(i) + dkp_xy
            
            !< Efforts de rappel (écrouissage cinématique, éq. 4.31)
            ! Nrx_old  = Cm11 * eps0p_xx(i) + Cm12 * eps0p_yy(i)
            ! Nry_old  = Cm21 * eps0p_xx(i) + Cm22 * eps0p_yy(i)
            ! Mrx_old  = Cb11 * kap_p_xx(i) + Cb12 * kap_p_yy(i)
            ! Mry_old  = Cb21 * kap_p_xx(i) + Cb22 * kap_p_yy(i)
            ! Mrxy_old = Cb33 * kap_p_xy(i)
            
            !< Efforts finaux (éq. 6.23)
            Nx  = Nx  - Dm11*deps0p_x - Dm12*deps0p_y
            Ny  = Ny  - Dm21*deps0p_x - Dm22*deps0p_y
            Mx  = Mx  - Db11(i)*dkp_x - Db12(i)*dkp_y - Db13(i)*dkp_xy
            My  = My  - Db21(i)*dkp_x - Db22(i)*dkp_y - Db23(i)*dkp_xy
            Mxy = Mxy - Db31(i)*dkp_x - Db32(i)*dkp_y - Db33(i)*dkp_xy

          enddo            
!
          !=====================================================================
        end subroutine sigeps136g
      end module sigeps136g_mod
