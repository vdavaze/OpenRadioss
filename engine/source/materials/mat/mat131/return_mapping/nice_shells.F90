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
      module nice_shells_mod
      contains
      subroutine nice_shells(                                                  &
        nel      ,matparam ,rho      ,nvartmp  ,vartmp   ,                     &
        depsxx   ,depsyy   ,depsxy   ,depsyz   ,depszx   ,                     &
        sigoxx   ,sigoyy   ,sigoxy   ,sigoyz   ,sigozx   ,                     &
        signxx   ,signyy   ,signxy   ,signyz   ,signzx   ,                     &
        soundsp  ,off      ,pla      ,dpla     ,seq      ,et       ,           &
        sigy     ,timestep ,epsd     ,temp     ,shf      ,thk      ,thkly    , &
        asrate   ,l_sigb   ,sigb     ,epsd_pg  ,nuvar    ,uvar     )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use constant_mod
        use precision_mod, only : WP
        use elasto_plastic_trial_stress_mod
        use elasto_plastic_eq_stress_mod
        use elasto_plastic_yield_stress_mod
        use elasto_plastic_kinematic_hardening_mod
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        integer,                       intent(in)    :: nel       !< Number of elements in the group
        type(matparam_struct_),        intent(in)    :: matparam  !< Material parameters data
        real(kind=WP), dimension(nel), intent(in)    :: rho       !< Density at current time
        integer,                       intent(in)    :: nvartmp   !< Number of variables used in tabulated variables
        integer, dimension(nel,nvartmp), intent(inout) :: vartmp  !< Temporary variables for tabulated hardening
        real(kind=WP), dimension(nel), intent(in)    :: depsxx    !< Strain increment xx
        real(kind=WP), dimension(nel), intent(in)    :: depsyy    !< Strain increment yy
        real(kind=WP), dimension(nel), intent(in)    :: depsxy    !< Strain increment xy
        real(kind=WP), dimension(nel), intent(in)    :: depsyz    !< Strain increment yz
        real(kind=WP), dimension(nel), intent(in)    :: depszx    !< Strain increment zx
        real(kind=WP), dimension(nel), intent(in)    :: sigoxx    !< Previous stress xx
        real(kind=WP), dimension(nel), intent(in)    :: sigoyy    !< Previous stress yy
        real(kind=WP), dimension(nel), intent(in)    :: sigoxy    !< Previous stress xy
        real(kind=WP), dimension(nel), intent(in)    :: sigoyz    !< Previous stress yz
        real(kind=WP), dimension(nel), intent(in)    :: sigozx    !< Previous stress zx
        real(kind=WP), dimension(nel), intent(inout) :: signxx    !< Current stress xx
        real(kind=WP), dimension(nel), intent(inout) :: signyy    !< Current stress yy
        real(kind=WP), dimension(nel), intent(inout) :: signxy    !< Current stress xy
        real(kind=WP), dimension(nel), intent(inout) :: signyz    !< Current stress yz
        real(kind=WP), dimension(nel), intent(inout) :: signzx    !< Current stress zx
        real(kind=WP), dimension(nel), intent(inout) :: soundsp   !< Current sound speed
        real(kind=WP), dimension(nel), intent(inout) :: off       !< Element failure flag
        real(kind=WP), dimension(nel), intent(inout) :: pla       !< Accumulated plastic strain
        real(kind=WP), dimension(nel), intent(inout) :: dpla      !< Plastic strain increment
        real(kind=WP), dimension(nel), intent(inout) :: seq       !< Equivalent stress
        real(kind=WP), dimension(nel), intent(inout) :: et        !< Hourglass stabilization variable
        real(kind=WP), dimension(nel), intent(inout) :: sigy      !< Current yield stress
        real(kind=WP), intent(in)                    :: timestep  !< Time step
        real(kind=WP), dimension(nel), intent(inout) :: epsd      !< Plastic strain rate
        real(kind=WP), dimension(nel), intent(inout) :: temp      !< Temperature
        real(kind=WP), dimension(nel), intent(in)    :: shf       !< Shear correction factor
        real(kind=WP), dimension(nel), intent(inout) :: thk       !< Current thickness
        real(kind=WP), dimension(nel), intent(in)    :: thkly     !< Integration point layer thickness
        real(kind=WP),                 intent(in)    :: asrate    !< Strain rate filtering weighting factor
        integer,                       intent(in)    :: l_sigb    !< Size of backstress array
        real(kind=WP),dimension(nel,l_sigb),intent(inout) :: sigb !< Backstress components for kinematic hardening
        real(kind=WP),dimension(nel),  intent(in)    :: epsd_pg   !< Global equivalent strain rate
        integer,                       intent(in)    :: nuvar     !< Number of user variables
        real(kind=WP),dimension(nel,nuvar), intent(inout) :: uvar  !< User variables
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: i,j,ii,nindx,indx(nel),vpflag,ikine
        real(kind=WP), dimension(nel,6,6) :: cstf,N
        real(kind=WP) :: dlam,dsigxx_dlam,dsigyy_dlam,dsigxy_dlam,dseq_dlam,   &
          dpla_dlam,dphi_dseq,dphi_dsigy,dphi_dlam,sig_dseqdsig,dsigy_dlam,    &
          dsigbxx_dlam,dsigbyy_dlam,dsigbzz_dlam,dsigbxy_dlam,chard,dphi
        real(kind=WP), dimension(nel) :: pla0,normxx,normyy,normzz,normxy,     &
          normyz,normzx,phi,young,dsigy_dpla,dtemp_dpla,s13,s23,depzz,         &
          sigbxx,sigbyy,sigbzz,sigbxy,sigy0,dsigy0_dpla,dtemp0_dpla,zeros,     &
          dsigxx,dsigyy,dsigxy,phi0,sig0xx,sig0yy,sig0zz,sig0xy,sig0yz,sig0zx
        real(kind=WP), dimension(nel,l_sigb) :: dsigb_dlam
        real(kind=WP), dimension(nel) :: signzz,sigozz,depszz,dezz
        integer, dimension(nel,nvartmp) :: ipos0
!
        integer, parameter :: eltype = 2               !< Element type (1 - Solids, 2 - Shells)
        integer, parameter :: iresp = 0                !< Response type (0 - standard)
        zeros(1:nel) = zero
!===============================================================================
!
        !=======================================================================
        !< - Initialisation of computation on time step
        !=======================================================================
        !< Viscoplastic formulation flag
        vpflag = matparam%iparam(10)
        !< Total strain-rate computation
        if (vpflag > 1) then
          epsd(1:nel) = asrate*epsd_pg(1:nel) + (one-asrate)*epsd(1:nel)
        !< Plastic strain rate recovering
        else
          epsd(1:nel) = uvar(1:nel,1)
        endif
        !< Kinematic hardening flag
        ikine = matparam%iparam(22)
        !< Mixed kinematic/isotropic hardening parameter
        chard = matparam%uparam(matparam%iparam(20) + 1)
        !< Initialisation of the hourglass control variable
        et(1:nel) = one
        !< Increment of cumulated plastic strain
        dpla(1:nel) = zero
        !< Derivative of temperature w.r.t. cumulated plastic strain
        dtemp_dpla(1:nel) = zero
        !< Save the initial cumulated plastic strain value
        pla0(1:nel) = pla(1:nel)
        !< Initialize out-of-plane plastic strain increment for shell elements
        depzz(1:nel) = zero
        !< Initialize dummy stress components for shell elements
        signzz(1:nel) = zero
        sigozz(1:nel) = zero
        depszz(1:nel) = zero
        !< Save initial values of the stress tensor components for shell elements
        sig0xx(1:nel) = sigoxx(1:nel)
        sig0yy(1:nel) = sigoyy(1:nel)
        sig0zz(1:nel) = sigozz(1:nel)
        sig0xy(1:nel) = sigoxy(1:nel)
        sig0yz(1:nel) = sigoyz(1:nel)
        sig0zx(1:nel) = sigozx(1:nel)
        !< Recover previous value of the yield function
        phi0(1:nel) = uvar(1:nel,2)
        !< Initial index array 
        nindx = nel
        do i = 1,nel
          indx(i) = i
        enddo
!
        !=======================================================================
        !< - Computation of the elastic trial stress tensor
        !=======================================================================
        call elasto_plastic_trial_stress(                                      &
          matparam ,nel      ,soundsp  ,cstf     ,young    ,rho      ,         &
          depsxx   ,depsyy   ,depszz   ,depsxy   ,depsyz   ,depszx   ,         &
          sigoxx   ,sigoyy   ,sigozz   ,sigoxy   ,sigoyz   ,sigozx   ,         &
          signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,         &
          eltype   ,shf      ,s13      ,s23      )
!
        !=======================================================================
        !< - Computation of the initial yield stress
        !=======================================================================
        call elasto_plastic_yield_stress(                                      &
          matparam ,nel      ,nindx    ,indx      ,sigy     ,pla      ,        &
          epsd     ,dsigy_dpla,nvartmp  ,vartmp   ,temp     ,dtemp_dpla)
!
        !=======================================================================
        !< - Backstress tensor computation for kinematic hardening models
        !=======================================================================
        if (ikine > 0) then
          do i=1,nel
            sigbxx(i) = zero
            sigbyy(i) = zero
            sigbzz(i) = zero
            sigbxy(i) = zero
            !< Compute the backstress tensor from all C-R kinematic hardenings
            do j = 1, l_sigb/6
              sigbxx(i) = sigbxx(i) + sigb(i,6*(j-1) + 1)
              sigbyy(i) = sigbyy(i) + sigb(i,6*(j-1) + 2)
              sigbzz(i) = sigbzz(i) + sigb(i,6*(j-1) + 3)
              sigbxy(i) = sigbxy(i) + sigb(i,6*(j-1) + 4)
            enddo
            !< Add the kinematic hardening contribution to stress tensors
            signxx(i) = signxx(i) - (sigbxx(i) - sigbzz(i))
            signyy(i) = signyy(i) - (sigbyy(i) - sigbzz(i))
            signxy(i) = signxy(i) - sigbxy(i)
            sig0xx(i) = sig0xx(i) - (sigbxx(i) - sigbzz(i))
            sig0yy(i) = sig0yy(i) - (sigbyy(i) - sigbzz(i))
            sig0xy(i) = sig0xy(i) - sigbxy(i)
          enddo
          !< Initial yield stress computation for kinematic hardening models
          zeros(1:nel) = zero
          ipos0(1:nel,1:nvartmp) = 0
          call elasto_plastic_yield_stress(                                    &
            matparam ,nel      ,nindx    ,indx      ,sigy0    ,zeros    ,      &
            epsd     ,dsigy0_dpla,nvartmp,ipos0     ,temp     ,dtemp0_dpla)
          !< Update of the yield stress for kinematic hardening models
          sigy(1:nel) = (one - chard)*sigy(1:nel) + chard*sigy0(1:nel)
        endif
!
        !=======================================================================
        !< - Computation of the trial equivalent stress and its 1st derivative
        !=======================================================================
        call elasto_plastic_eq_stress(                                         &
          matparam ,nel      ,nindx    ,indx     ,iresp    ,eltype   ,         &
          signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,         &
          normxx   ,normyy   ,normzz   ,normxy   ,normyz   ,normzx   ,         &
          N        ,.false.  ,seq      )
!
        !=======================================================================
        !< - Computation of the trial yield function and count yielding elements
        !=======================================================================
        nindx = 0
        do i = 1,nel
          indx(i) = 0
          phi(i)  = (seq(i)/sigy(i))**2 - one
        enddo
        do i = 1,nel
          if ((phi(i) >= zero).and.(off(i) == one)) then
            nindx = nindx + 1
            indx(nindx) = i
          endif
        enddo
!
        !=======================================================================
        !< - Return mapping 1 step using Next Increment Correct Error (N.I.C.E)
        !=======================================================================
        if (nindx > 0) then
!
          ! Note     : in this part, the purpose is to compute in one iteration
          ! a plastic multiplier allowing to update internal variables to satisfy
          ! the consistency condition. 
          ! Its expression is : dlam = - (phi0 + dphi)/dphi_dlam
          ! -> phi0 : old value of yield function (known)
          ! -> dphi : yield function prediction (to compute)
          ! -> dphi_dlambda : derivative of phi with respect to dlambda by taking
          !                   into account of internal variables kinetic : 
          !                plasticity, strain-rate ... (to be computed)
!
          !< Computation of equivalent stress of the previous stress tensor
          call elasto_plastic_eq_stress(                                       &
            matparam ,nel      ,nindx    ,indx     ,iresp    ,eltype   ,       &
            sig0xx   ,sig0yy   ,sig0zz   ,sig0xy   ,sig0yz   ,sig0zx   ,       &
            normxx   ,normyy   ,normzz   ,normxy   ,normyz   ,normzx   ,       &
            N        ,.false.  ,seq      )
!
          !< Loop over yielding elements
#include "vectorize.inc"
          do ii = 1,nindx
            i = indx(ii)
!
            !< 1 - Derivative of equivalent stress sigeq w.r.t lambda
            !< -----------------------------------------------------------------
!
            !<  a) Derivatives of stress tensor w.r.t lambda
            !<  ----------------------------------------------------------------
            dsigxx_dlam = -(cstf(i,1,1)*normxx(i) + cstf(i,1,2)*normyy(i))
            dsigyy_dlam = -(cstf(i,1,2)*normxx(i) + cstf(i,2,2)*normyy(i))
            dsigxy_dlam = -(cstf(i,4,4)*normxy(i))
!
            !<  b) Assembling derivative of eq. stress sigeq w.r.t lambda
            !<  ----------------------------------------------------------------
            dseq_dlam = normxx(i)*dsigxx_dlam + normyy(i)*dsigyy_dlam +        &
                        normxy(i)*dsigxy_dlam
!
            !< 2 - Derivative of yield stress ystrs w.r.t lambda
            !< -----------------------------------------------------------------
!
            !<  a) Derivative of eq. plastic strain w.r.t lambda
            !<  ----------------------------------------------------------------
            sig_dseqdsig = signxx(i)*normxx(i) +                               &
                           signyy(i)*normyy(i) +                               &
                           signxy(i)*normxy(i)
            dpla_dlam = sig_dseqdsig/max(sigy(i),em20)
!
            !<  b) Assembling derivative of ystrs w.r.t lambda
            !<  ----------------------------------------------------------------
            dsigy_dlam = (one - chard)*dsigy_dpla(i)*dpla_dlam
!
            !< 3 - Add kinematic hardening to the derivative of eq.stress 
            !   w.r.t lambda
            !<  ----------------------------------------------------------------              
            if (ikine > 0) then
              !< a - Derivative of backstress tensor w.r.t lambda
              !<  --------------------------------------------------------------
              call elasto_plastic_kinematic_hardening(                         &
                matparam ,1     ,l_sigb ,dsigb_dlam(i,1),dsigy_dpla(i),chard,  &
                normxx(i),normyy(i),normzz(i),normxy(i),normyz(i),normzx(i),   &
                dpla_dlam,sigb(i,1),eltype   ) 
              !< b - Assembling the backstress contribution to the derivative  
              !  of eq. stress w.r.t lambda
              !<  --------------------------------------------------------------
              dsigbxx_dlam = zero
              dsigbyy_dlam = zero
              dsigbzz_dlam = zero
              dsigbxy_dlam = zero
              do j = 1, l_sigb/6
                dsigbxx_dlam = dsigbxx_dlam + dsigb_dlam(i,6*(j-1) + 1)
                dsigbyy_dlam = dsigbyy_dlam + dsigb_dlam(i,6*(j-1) + 2)
                dsigbzz_dlam = dsigbzz_dlam + dsigb_dlam(i,6*(j-1) + 3)
                dsigbxy_dlam = dsigbxy_dlam + dsigb_dlam(i,6*(j-1) + 4)
              enddo
              dseq_dlam = dseq_dlam - normxx(i)*dsigbxx_dlam -                 &
                                      normyy(i)*dsigbyy_dlam +                 &
                            (normxx(i)+normyy(i))*dsigbzz_dlam -               &
                                      normxy(i)*dsigbxy_dlam
            endif
!
            !< 4 - Assembling the derivative of phi w.r.t lambda
            !< -----------------------------------------------------------------
!
            !<  a) Derivative of phi w.r.t eq. stress sigeq
            !<  ----------------------------------------------------------------
            dphi_dseq  =  two*seq(i)/(sigy(i)**2)
!
            !<  b) Derivative of phi w.r.t yield stress ystrs
            !<  ----------------------------------------------------------------
            dphi_dsigy = -two*(seq(i)**2)/(sigy(i)**3)
!
            !<  c) Derivative of phi w.r.t lambda
            !<  ----------------------------------------------------------------
            dphi_dlam = dphi_dseq*dseq_dlam + dphi_dsigy*dsigy_dlam
            dphi_dlam = sign(max(abs(dphi_dlam),em20),dphi_dlam)
!
            !< 5 - Computation of plastic multiplier and variables update
            !< -----------------------------------------------------------------
!
            !<  a) Computation of the plastic multiplier increment dlam
            !<  ----------------------------------------------------------------
            !< Computation of the trial stress increment
            dsigxx(i) = signxx(i) - sig0xx(i)
            dsigyy(i) = signyy(i) - sig0yy(i)
            dsigxy(i) = signxy(i) - sig0xy(i) 
            !< Computation of yield surface trial increment dphi       
            dphi = dphi_dseq * (                                               &
                   normxx(i) * dsigxx(i)                                       &
                 + normyy(i) * dsigyy(i)                                       &
                 + normxy(i) * dsigxy(i))
            !< Assembling plastic multiplier
            dlam = -(phi0(i) + dphi)/dphi_dlam
!
            !<  b) Stress tensor update
            !<  ----------------------------------------------------------------
            signxx(i) = signxx(i) + dsigxx_dlam*dlam
            signyy(i) = signyy(i) + dsigyy_dlam*dlam
            signxy(i) = signxy(i) + dsigxy_dlam*dlam
!
            !<  c) Update the plastic strain related variables
            !<  ----------------------------------------------------------------
            !< Equivalent plastic strain increment
            dpla(i) = max(dpla(i) + dpla_dlam*dlam,zero)
            !< Equivalent plastic strain
            pla(i)  = pla0(i) + dpla(i)
            !< Temperature
            temp(i) = temp(i) + dtemp_dpla(i)*dpla_dlam*dlam
            !< Out-of-plane plastic strain for shell elements
            depzz(i) = depzz(i) + dlam*normzz(i)
          enddo
!
          !<  d) Yield stress update
          !<  ----------------------------------------------------------------
          call elasto_plastic_yield_stress(                                    &
            matparam ,nel      ,nindx    ,indx      ,sigy     ,pla      ,      &
            epsd     ,dsigy_dpla,nvartmp ,vartmp    ,temp     ,dtemp_dpla)
!
          !<  e) Backstress tensor update
          !<  ----------------------------------------------------------------
          !< Update of the backstress tensor (if kinematic hardening)
          if (ikine > 0) then
            !< Loop over yielding elements
#include "vectorize.inc"
            do ii = 1,nindx
              i = indx(ii)
              ! -> Remove kinematic hardening contribution
              signxx(i) = signxx(i) + (sigbxx(i) - sigbzz(i))
              signyy(i) = signyy(i) + (sigbyy(i) - sigbzz(i))
              signxy(i) = signxy(i) + sigbxy(i)
              ! -> Add the evolution of backstress tensor
              sigbxx(i) = sigbxx(i) + dsigbxx_dlam*dlam
              sigbyy(i) = sigbyy(i) + dsigbyy_dlam*dlam
              sigbzz(i) = sigbzz(i) + dsigbzz_dlam*dlam
              sigbxy(i) = sigbxy(i) + dsigbxy_dlam*dlam
              ! -> Add the kinematic hardening contribution
              signxx(i) = signxx(i) - (sigbxx(i) - sigbzz(i))
              signyy(i) = signyy(i) - (sigbyy(i) - sigbzz(i))
              signxy(i) = signxy(i) - sigbxy(i)
              ! -> Update of the backstress components
              do j = 1, l_sigb/6
                sigb(i,6*(j-1) + 1) = sigb(i,6*(j-1) + 1) +                    &
                                     dsigb_dlam(i,6*(j-1) + 1)*dlam
                sigb(i,6*(j-1) + 2) = sigb(i,6*(j-1) + 2) +                    &
                                     dsigb_dlam(i,6*(j-1) + 2)*dlam
                sigb(i,6*(j-1) + 3) = sigb(i,6*(j-1) + 3) +                    &
                                     dsigb_dlam(i,6*(j-1) + 3)*dlam
                sigb(i,6*(j-1) + 4) = sigb(i,6*(j-1) + 4) +                    &
                                     dsigb_dlam(i,6*(j-1) + 4)*dlam
              enddo
              !< Update of the yield stress for kinematic hardening models
              sigy(i) = (one - chard)*sigy(i) + chard*sigy0(i)
            enddo
          endif
!
          !<  f) Equivalent stress update
          !<  ----------------------------------------------------------------
          call elasto_plastic_eq_stress(                                     &
            matparam ,nel      ,nindx    ,indx     ,iresp    ,eltype   ,     &
            signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,     &
            normxx   ,normyy   ,normzz   ,normxy   ,normyz   ,normzx   ,     &
            N        ,.false.  ,seq      )
!
          !< Loop over yielding elements
#include "vectorize.inc"
          do ii = 1, nindx
            i = indx(ii)
            !<  g) Yield function update
            !<  ----------------------------------------------------------------
            phi(i) = (seq(i)/sigy(i))**2 - one
            !< Update the hourglass stabilization variable
            et(i) = dsigy_dpla(i) / (dsigy_dpla(i) + young(i))
          enddo
!
        endif
!
        !=======================================================================
        !< - Update filtered plastic strain rate
        !=======================================================================
        if (vpflag == 1) then 
          do i = 1,nel
            epsd(i) = asrate*dpla(i)/max(timestep,em20) + (one - asrate)*uvar(i,1)
            uvar(i,1) = epsd(i)
          enddo
        endif
!
        !=======================================================================
        !< - Remove backstress contribution from stress tensor
        !=======================================================================        
        if (ikine > 0) then
          do i=1,nel
            signxx(i) = signxx(i) + (sigbxx(i) - sigbzz(i))
            signyy(i) = signyy(i) + (sigbyy(i) - sigbzz(i))
            signxy(i) = signxy(i) + sigbxy(i)
          enddo
        endif
!        
        !=======================================================================
        !< - Update remaining error after return mapping
        !=======================================================================
        do i = 1,nel
          if (dpla(i) > zero) then 
            uvar(i,2) = phi(i)
          else
            uvar(i,2) = zero
          endif
        enddo
!
        !=======================================================================
        !< Update thickness for shell elements
        !=======================================================================
        do i = 1,nel
          dezz(i) = - s13(i)*(signxx(i) - sigoxx(i))                           &
                    - s23(i)*(signyy(i) - sigoyy(i))                           &
                    + depzz(i)
          thk(i)  = thk(i) + dezz(i)*thkly(i)*off(i) 
        enddo
!
       end subroutine nice_shells
       end module nice_shells_mod
