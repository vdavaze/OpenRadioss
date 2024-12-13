!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2024 Altair Engineering Inc.
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
      module sigeps81_mod
        contains
! ======================================================================================================================
! \brief   Drucker-Prager with cap hardening material law /MAT/LAW81 (DPRAG_CAP)
! \details Material law based on Drucker-Prager pressure dependent plastic model with cap hardening. 
! ======================================================================================================================
        subroutine sigeps81 (                                                  &
           nel     ,nuvar   ,uvar     ,matparam ,nfunc    ,ifunc   ,           &
           npf     ,tf      ,snpc     ,stf      ,rho0     ,rho     ,           &
           volume  ,amu     ,defp     ,soundsp  ,viscmax  ,dt1     ,           &
           depsxx  ,depsyy  ,depszz   ,depsxy   ,depsyz   ,depszx  ,           &
           sigoxx  ,sigoyy  ,sigozz   ,sigoxy   ,sigoyz   ,sigozx  ,           &
           signxx  ,signyy  ,signzz   ,signxy   ,signyz   ,signzx  ,           &
           sigvxx  ,sigvyy  ,sigvzz   ,nvartmp  ,vartmp   ,seq     )
!-------------------------------------------------------------------------------
!   M o d u l e s
!-------------------------------------------------------------------------------
          use matparam_def_mod 
          use constant_mod      
!-------------------------------------------------------------------------------
!   I m p l i c i t   T y p e s
!-------------------------------------------------------------------------------
          implicit none 
#include  "my_real.inc"
!-------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-------------------------------------------------------------------------------
          integer, intent(in) :: nel                           !< number of elements in the group
          integer, intent(in) :: nuvar                         !< number of user variables
          my_real, dimension(nel,nuvar), intent(inout) :: uvar !< user variables
          type(matparam_struct_), intent(in) :: matparam       !< material parameters data
          integer, intent(in) :: nfunc                         !< number of functions
          integer, dimension(nfunc), intent(in) :: ifunc       !< function index
          integer, dimension(snpc), intent(in) :: npf          !< number of values for each function
          my_real, dimension(stf), intent(in) :: tf            !< function parameters
          integer, intent(in) :: snpc                          !< number of parameters for each function
          integer, intent(in) :: stf                           !< number of parameters
          my_real, dimension(nel), intent(in) :: rho0          !< initial density
          my_real, dimension(nel), intent(in) :: rho           !< current density
          my_real, dimension(nel), intent(in) :: volume        !< integration point associated volume
          my_real, dimension(nel), intent(in) :: amu           !< element volumetric strain
          my_real, dimension(nel,2), intent(inout) :: defp     !< element plastic strain
          my_real, dimension(nel), intent(inout) :: soundsp    !< sound speed
          my_real, dimension(nel), intent(inout) :: viscmax    !< maximum viscosity
          my_real, intent(in) :: dt1                           !< time step
          my_real, dimension(nel), intent(in)  :: depsxx       !< strain increment component xx
          my_real, dimension(nel), intent(in)  :: depsyy       !< strain increment component yy
          my_real, dimension(nel), intent(in)  :: depszz       !< strain increment component zz
          my_real, dimension(nel), intent(in)  :: depsxy       !< strain increment component xy
          my_real, dimension(nel), intent(in)  :: depsyz       !< strain increment component yz
          my_real, dimension(nel), intent(in)  :: depszx       !< strain increment component zx
          my_real, dimension(nel), intent(in)  :: sigoxx       !< old stress component xx
          my_real, dimension(nel), intent(in)  :: sigoyy       !< old stress component yy
          my_real, dimension(nel), intent(in)  :: sigozz       !< old stress component zz
          my_real, dimension(nel), intent(in)  :: sigoxy       !< old stress component xy
          my_real, dimension(nel), intent(in)  :: sigoyz       !< old stress component yz
          my_real, dimension(nel), intent(in)  :: sigozx       !< old stress component zx
          my_real, dimension(nel), intent(out) :: signxx       !< new stress component xx
          my_real, dimension(nel), intent(out) :: signyy       !< new stress component yy
          my_real, dimension(nel), intent(out) :: signzz       !< new stress component zz
          my_real, dimension(nel), intent(out) :: signxy       !< new stress component xy
          my_real, dimension(nel), intent(out) :: signyz       !< new stress component yz
          my_real, dimension(nel), intent(out) :: signzx       !< new stress component zx 
          my_real, dimension(nel), intent(out) :: sigvxx       !< viscous stress component xx
          my_real, dimension(nel), intent(out) :: sigvyy       !< viscous stress component yy
          my_real, dimension(nel), intent(out) :: sigvzz       !< viscous stress component zz
          integer, intent(in) :: nvartmp                       !< number of temporary variables
          integer, dimension(nel,nvartmp), intent(inout) :: vartmp !< temporary variables
          my_real, dimension(nel), intent(inout) :: seq        !< Von Mises equivalent stress
!-------------------------------------------------------------------------------
!  L o c a l   V a r i a b l e s
!-------------------------------------------------------------------------------
          integer :: i,ii,soft_flag,iter,nindx
          integer, dimension(nel) :: indx,ipos,iad,ilen
          integer, parameter :: niter = 3
          my_real :: kini,gini,tgphi,tgpsi,alpha,max_dilat,epspvol0,kwater,    &
            por0,sat0,u0,tolmu,viscfac,cini,capini
          my_real :: delta,depspd_dlam,depspv_dlam,df_dlam,dfdc,dfdp,dfdrc,    &
            dfdseq,dfdsig_dsigdlam,dfdsigxx,dfdsigyy,dfdsigzz,dfdsigxy,        &
            dfdsigyz,dfdsigzx,dgdp,dgdseq,dgdsigxx,dgdsigyy,dgdsigzz,          &
            dgdsigxy,dgdsigyz,dgdsigzx,dlam,dpdsigxx,dpdsigyy,dpdsigzz,        &
            dpdsigxy,dpdsigyz,dpdsigzx,drcdpa,drcdpb,dseqdsigxx,dseqdsigyy,    &
            dseqdsigzz,dseqdsigxy,dseqdsigyz,dseqdsigzx,dsigxxdlam,dsigyydlam, &
            dsigzzdlam,dsigxydlam,dsigyzdlam,dsigzxdlam,drcdp,ldav,trdep,      &
            trdgpds,fac
          my_real, dimension(nel) :: k,g,lame,g2,c,pa,pb,p0,dcdepsp,dpbdepsp,  &
            deri,epspd,epspv,depspd,depspv,dcdepspd,dpbdepspv,dpadepspv,f,     &
            p,sxx,syy,szz,sxy,syz,szx,pu,rc,a,dpxx,dpyy,dpzz,dpxy,dpyz,dpzx,   &
            muw,u,dudmu,por
!=============================================================================== 
!
          !=====================================================================
          ! - INITIALISATION OF COMPUTATION ON TIME STEP
          !=====================================================================
          !< Recovering integer model parameter
          soft_flag = matparam%iparam(1)  !< Softening flag
          !< Recovering real model paramter 
          kini      = matparam%bulk       !< Initial bulk modulus
          gini      = matparam%shear      !< Initial shear modulus
          tgphi     = matparam%uparam(1)  !< Friction angle
          tgpsi     = matparam%uparam(2)  !< Plastic flow potential angle
          cini      = matparam%uparam(3)  !< Initial material cohesion
          capini    = matparam%uparam(4)  !< Initial cap limit pressure
          alpha     = matparam%uparam(5)  !< Ratio Pa/Pb
          max_dilat = matparam%uparam(6)  !< Maximum dilatancy
          epspvol0  = matparam%uparam(5)  !< Initial volumetric plastic strain
          kwater    = matparam%uparam(6)  !< Pore water bulk modulus
          por0      = matparam%uparam(7)  !< Initial porosity
          sat0      = matparam%uparam(8)  !< Initial saturation
          u0        = matparam%uparam(9)  !< Initial pore water pressure
          tolmu     = matparam%uparam(10) !< Tolerance for capshift viscosity
          viscfac   = matparam%uparam(11) !< Viscosity factor
!
          !=====================================================================
          !< - Recovering user variables and state variables
          !=====================================================================
          epspd(1:nel)  = defp(1:nel,1) !< Deviatoric Equivalent Plastic Strain  
          epspv(1:nel)  = defp(1:nel,2) !< Volumetric Plastic Strain 
          depspv(1:nel) = zero
          depspd(1:nel) = zero
          dpxx(1:nel)   = zero
          dpyy(1:nel)   = zero
          dpzz(1:nel)   = zero
          dpxy(1:nel)   = zero
          dpyz(1:nel)   = zero
          dpzx(1:nel)   = zero
!
          !=====================================================================
          !< - Recovering elastic parameters
          !=====================================================================
          !<  Bulk modulus computation
          !  -> Interpolated with volumetric plastic strain
          if (ifunc(1) > 0) then
            ipos(1:nel) = vartmp(1:nel,1)
            iad(1:nel)  = npf(ifunc(1)) / 2 + 1
            ilen(1:nel) = npf(ifunc(1)+1) / 2 - iad(1:nel) - ipos(1:nel)
            call vinter2(tf,iad,ipos,ilen,nel,epspv,deri,k) 
            k(1:nel) = kini*k(1:nel)
            vartmp(1:nel,1) = ipos(1:nel)
          !  -> Constant value
          else
            k(1:nel) = kini
          endif
          !<  Shear modulus computation
          !   -> Interpolated with volumetric plastic strain
          if (ifunc(2) > 0) then
            ipos(1:nel) = vartmp(1:nel,2)
            iad(1:nel)  = npf(ifunc(2)) / 2 + 1
            ilen(1:nel) = npf(ifunc(2)+1) / 2 - iad(1:nel) - ipos(1:nel)
            call vinter2(tf,iad,ipos,ilen,nel,epspv,deri,g)
            g(1:nel) = gini*g(1:nel)
            vartmp(1:nel,2) = ipos(1:nel)
          !   -> Constant value
          else
            g(1:nel) = gini
          endif
          !< Two*shear modulus
          g2(1:nel) = two*g(1:nel)
          !< Lame coefficient
          lame(1:nel) = k(1:nel) - two_third*g(1:nel)
!
          !=====================================================================
          !< - Yield criterion parameter computation
          !=====================================================================
          !< Material cohesion
          !  -> Interpolated with deviatoric plastic strain
          if (ifunc(3) > 0) then 
            ipos(1:nel) = vartmp(1:nel,3)
            iad(1:nel)  = npf(ifunc(3)) / 2 + 1
            ilen(1:nel) = npf(ifunc(3)+1) / 2 - iad(1:nel) - ipos(1:nel)
            call vinter2(tf,iad,ipos,ilen,nel,epspd,dcdepspd,c)
            c(1:nel) = cini*c(1:nel)
            dcdepspd(1:nel) = cini*dcdepspd(1:nel)
            vartmp(1:nel,3) = ipos(1:nel)
          !  -> Constant value
          else
            c(1:nel) = cini
            dcdepspd(1:nel) = zero
          endif
          !< Cap limit pressure
          !  -> Interpolated with volumetric plastic strain
          if (ifunc(4) > 0) then
            ipos(1:nel) = vartmp(1:nel,4)
            iad(1:nel)  = npf(ifunc(4)) / 2 + 1
            ilen(1:nel) = npf(ifunc(4)+1) / 2 - iad(1:nel) - ipos(1:nel)
            call vinter2(tf,iad,ipos,ilen,nel,epspv,dpbdepspv,pb)
            pb(1:nel) = capini*pb(1:nel)
            dpbdepspv(1:nel) = capini*dpbdepspv(1:nel)
            vartmp(1:nel,4) = ipos(1:nel)
          !  -> Constant value
          else
            pb(1:nel) = capini
            dpbdepspv(1:nel) = zero
          endif
          !< Transition pressure yield surface to cap
          pa(1:nel) = alpha*pb(1:nel)  
          dpadepspv(1:nel) = alpha*dpbdepspv(1:nel)  
          !< Null criterion derivative pressure (dfdp = 0)
          do i = 1,nel
            delta = (pa(i)*tgphi + c(i))**2 + eight*((pb(i)-pa(i))**2)*(tgphi**2)
            if (tgphi > zero) then
              p0(i) = pa(i) + (-(pa(i)*tgphi+c(i)) + sqrt(delta))/(four*tgphi)
            else
              p0(i) = pa(i)
            endif
          enddo
! 
          !=====================================================================
          !< - Computation of trial stress tensor, Von Mises and pressure
          !=====================================================================
          do i=1,nel
            !< Trial Cauchy stress tensor
            ldav      = lame(i)*(depsxx(i) + depsyy(i) + depszz(i))
            signxx(i) = sigoxx(i) + g2(i)*depsxx(i) + ldav
            signyy(i) = sigoyy(i) + g2(i)*depsyy(i) + ldav
            signzz(i) = sigozz(i) + g2(i)*depszz(i) + ldav
            signxy(i) = sigoxy(i) +  g(i)*depsxy(i)
            signyz(i) = sigoyz(i) +  g(i)*depsyz(i)
            signzx(i) = sigozx(i) +  g(i)*depszx(i)
!    
            !< Trial devriatoric stress tensor
            p(i)   = -third*(signxx(i) + signyy(i) + signzz(i))
            sxx(i) = signxx(i) + p(i)
            syy(i) = signyy(i) + p(i)
            szz(i) = signzz(i) + p(i)
            sxy(i) = signxy(i)
            syz(i) = signyz(i)
            szx(i) = signzx(i)
!    
            !< Trial Von Mises stress
            seq(i) = three_half*(sxx(i)**2 + syy(i)**2 + szz(i)**2)            &
                   +      three*(sxy(i)**2 + syz(i)**2 + szx(i)**2)
            seq(i) = sqrt(seq(i))
          enddo
!
          !=====================================================================
          !< - Porosity computation (if activated)
          !=====================================================================
          if (sat0 > zero) then
            do i = 1,nel
              por(i) = one - (one - por0)*exp(epspv(i) - epspvol0)
              fac    = max(em03,por(i)/por0)
              muw(i) = sat0/fac*amu(i)
              if (muw(i) >= tolmu) then
                u(i) = kwater*muw(i)
                dudmu(i) = kwater
              elseif (muw(i) > -tolmu) then
                u(i) = (kwater/(four*tolmu))*(muw(i)+tolmu)**2
                dudmu(i) = (kwater/(two*tolmu))*(muw(i)+tolmu)
              else
                u(i) = zero
                dudmu(i) = zero
              endif
            enddo
          else
            muw(1:nel)   = -one
            u(1:nel)     = zero
            dudmu(1:nel) = zero
          endif 
!    
          !=====================================================================
          ! - Computation of yield function and check element behavior
          !=====================================================================
          nindx = 0
          indx(1:nel) = 0
          do i=1,nel 
            !< Taking into account pore water pressure         
            if (p(i) < p0(i)) then
              pu(i) = p(i)
            elseif (p(i) <= p0(i) + u(i)) then
              pu(i) = p0(i)
            else
              pu(i) = p(i) - u(i)
            endif
            !< Regular return mapping treatment
            if (p(i) > -c(i)/tgphi) then
              !< Transition surface to cap hardening factor
              if (p(i)<=pa(i)) then
                rc(i) = one
              else
                rc(i) = one - ((p(i)-pa(i))/(pb(i)-pa(i)))**2
                rc(i) = sqrt(max(rc(i),zero))
              endif
              !< Pressure dependent factor
              a(i) = max(zero,p(i)*tgphi + c(i))
              f(i) = seq(i) - rc(i)*a(i)
              if (f(i) > zero) then 
                nindx = nindx + 1
                indx(nindx) = i
              endif
            !< Tri-traction treatment (apex of the yield surface)
            else
              depspv(i) = -(p(i) + (c(i)/tgphi))/k(i)
              if (soft_flag == 1) then
                epspv(i) = epspv(i) + max(depspv(i),zero)
              else
                epspv(i) = epspv(i) + depspv(i)
              endif
              p(i) = -c(i)/tgphi
              seq(i) = zero
              signxx(i) = -p(i)
              signyy(i) = -p(i)
              signzz(i) = -p(i)
              signxy(i) = zero
              signyz(i) = zero
              signzx(i) = zero
            endif
          enddo
!
          !=====================================================================
          ! - PLASTIC CORRECTION WITH CUTTING PLANE (NEWTON-ITERATION) METHOD
          !=====================================================================
          if (nindx > 0) then
!
            !< Loop over the iterations   
            do iter = 1, niter
              !< Loop over yielding elements
              do ii = 1, nindx
                i = indx(ii)
!       
                ! Note: in this part, the purpose is to compute for each iteration
                ! a plastic multiplier allowing to update internal variables to 
                ! satisfy the consistency condition using the cutting plane method
                ! within an iterative procedure.
                ! Its expression at each iteration is : dlambda = - f/df_dlambda
                ! -> f       : current value of yield function (known)
                ! -> df_dlam : derivative of f with respect to dlambda by taking
                !              into account of internal variables kinetic : 
                !              plasticity, damage ... (to be computed)
!      
                !< 1 - Derivative of yield criterion w.r.t plastic multiplier
                !      Contribution of the stress tensor
                !---------------------------------------------------------------
                !< Derivative of Von Mises stress w.r.t stress tensor
                dseqdsigxx = three_half*sxx(i)/max(seq(i),em20)
                dseqdsigyy = three_half*syy(i)/max(seq(i),em20)
                dseqdsigzz = three_half*szz(i)/max(seq(i),em20)
                dseqdsigxy =      three*sxy(i)/max(seq(i),em20)
                dseqdsigyz =      three*syz(i)/max(seq(i),em20)
                dseqdsigzx =      three*szx(i)/max(seq(i),em20)
!
                !< Derivative of hydrostatic pressure w.r.t stress tensor
                dpdsigxx = -third
                dpdsigyy = -third
                dpdsigzz = -third
                dpdsigxy =   zero
                dpdsigyz =   zero
                dpdsigzx =   zero
!
                !< Derivative of yield criterion with respect to pressure and 
                !  Von Mises stress
                dfdp = -rc(i)*tgphi
                if ((p(i) > pa(i)).and.(rc(i) > zero)) then 
                  drcdp = -(p(i) - pa(i))/(rc(i)*(pb(i) - pa(i))**2)
                  dfdp  = dfdp - drcdp*a(i)
                endif
                dfdseq = one
!
                !< Derivative of plastic potential with respect to pressure and 
                !  Von Mises stress
                if (p(i) <= pa(i)) then
                  dgdp = -tgpsi
                elseif (p(i) <= p0(i)) then
                  dgdp = -tgpsi*(one - ((p(i) - pa(i))/(p0(i) - pa(i))))
                elseif (p(i) > p0(i)) then 
                  dgdp = dfdp
                endif  
                dgdseq = one
!
                !< Check if maximum dilatancy is reached 
                !  (maximum dilatancy always negative)
                if (rho(i) <= (one + max_dilat)*rho0(i)) then
                  dgdp = max(zero,dgdp)
                  dfdp = max(zero,dfdp)
                endif
!
                !< Assembling derivative of yield criterion w.r.t stress tensor
                dfdsigxx = dfdseq*dseqdsigxx + dfdp*dpdsigxx
                dfdsigyy = dfdseq*dseqdsigyy + dfdp*dpdsigyy
                dfdsigzz = dfdseq*dseqdsigzz + dfdp*dpdsigzz
                dfdsigxy = dfdseq*dseqdsigxy + dfdp*dpdsigxy
                dfdsigyz = dfdseq*dseqdsigyz + dfdp*dpdsigyz
                dfdsigzx = dfdseq*dseqdsigzx + dfdp*dpdsigzx
!
                !< Assembling derivative of plastic potential w.r.t stress tensor
                dgdsigxx = dgdseq*dseqdsigxx + dgdp*dpdsigxx
                dgdsigyy = dgdseq*dseqdsigyy + dgdp*dpdsigyy
                dgdsigzz = dgdseq*dseqdsigzz + dgdp*dpdsigzz
                dgdsigxy = dgdseq*dseqdsigxy + dgdp*dpdsigxy
                dgdsigyz = dgdseq*dseqdsigyz + dgdp*dpdsigyz
                dgdsigzx = dgdseq*dseqdsigzx + dgdp*dpdsigzx
!
                !< Derivative of stress tensor w.r.t plastic multiplier
                trdgpds    = dgdsigxx + dgdsigyy + dgdsigzz  
                dsigxxdlam = -(dgdsigxx*g2(i) + lame(i)*trdgpds)
                dsigyydlam = -(dgdsigyy*g2(i) + lame(i)*trdgpds)
                dsigzzdlam = -(dgdsigzz*g2(i) + lame(i)*trdgpds) 
                dsigxydlam = -  dgdsigxy*g(i)
                dsigyzdlam = -  dgdsigyz*g(i)
                dsigzxdlam = -  dgdsigzx*g(i)   
!
                !< Contribution of the stress tensor to the derivative of 
                !  the yield criterion with respect to the plastic multiplier
                dfdsig_dsigdlam = dfdsigxx*dsigxxdlam + dfdsigyy*dsigyydlam +  &  
                                  dfdsigzz*dsigzzdlam + dfdsigxy*dsigxydlam +  &
                                  dfdsigyz*dsigyzdlam + dfdsigzx*dsigzxdlam
!
                !< 2 - Derivative of yield criterion w.r.t pressures Pa and Pb
                !---------------------------------------------------------------
                dfdrc = -one
                if ((p(i) > pa(i)).and.(rc(i) > zero)) then
                  drcdpb =  ((p(i) - pa(i))**2)/(rc(i)*(pb(i) - pa(i))**3)     
                  drcdpa = -((p(i) - pa(i))*(p(i) - pb(i)))/                   &
                            (rc(i)*(pb(i) - pa(i))**3)
                else
                  drcdpa = zero
                  drcdpb = zero
                endif
!
                !< 3 - Derivative of yield criterion w.r.t material cohesion
                !---------------------------------------------------------------
                dfdc = -rc(i)
!
                !< 4 - Derivative of deviatoric and volumetric plastic strain 
                !      w.r.t plastic multiplier
                !---------------------------------------------------------------
                depspd_dlam = (signxx(i)*dgdsigxx + signyy(i)*dgdsigyy +       &
                               signzz(i)*dgdsigzz + signxy(i)*dgdsigxy +       &
                               signyz(i)*dgdsigyz + signzx(i)*dgdsigzx)        &
                               /max(seq(i),em20)
                depspv_dlam = dgdsigxx + dgdsigyy + dgdsigzz
!
                !< 5 - Derivative of yield criterion w.r.t plastic multiplier
                !--------------------------------------------------------------- 
                df_dlam =  dfdsig_dsigdlam +                                   &
                           dfdrc*drcdpa*dpadepspv(i)*depspv_dlam +             &
                           dfdrc*drcdpb*dpbdepspv(i)*depspv_dlam +             &
                           dfdc*dcdepspd(i)*depspd_dlam
                df_dlam = sign(max(abs(df_dlam),em20),df_dlam)
!
                !< 6 - Computation of plastic multiplier
                !---------------------------------------------------------------             
                dlam = -f(i)/df_dlam
!
                !< 7 - Update plastic strain related variables
                !--------------------------------------------------------------- 
                !< Plastic strain tensor increment (on the current iteration)
                dpxx(i) = dlam * dgdsigxx
                dpyy(i) = dlam * dgdsigyy
                dpzz(i) = dlam * dgdsigzz
                dpxy(i) = dlam * dgdsigxy
                dpyz(i) = dlam * dgdsigyz
                dpzx(i) = dlam * dgdsigzx
                !< Volumetric plastic strain update
                depspv(i) = depspv(i) + depspv_dlam*dlam
                if (soft_flag == 1) then
                  epspv(i) = epspv(i) + max(depspv(i),zero)
                else
                  epspv(i) = epspv(i) + depspv(i)
                endif
                !< Deviatoric plastic strain update
                depspd(i) = max(depspd(i) + depspd_dlam*dlam,zero)
                epspd(i)  = epspd(i)  + depspd(i)
! 
                !< 8 - Update stress tensor, pressure and Von Mises stress
                !---------------------------------------------------------------
                trdep     = dpxx(i) + dpyy(i) + dpzz(i)
                signxx(i) = signxx(i) - (g2(i)*dpxx(i) + lame(i)*trdep)
                signyy(i) = signyy(i) - (g2(i)*dpyy(i) + lame(i)*trdep)
                signzz(i) = signzz(i) - (g2(i)*dpzz(i) + lame(i)*trdep)
                signxy(i) = signxy(i) -   g(i)*dpxy(i)
                signyz(i) = signyz(i) -   g(i)*dpyz(i)
                signzx(i) = signzx(i) -   g(i)*dpzx(i)
!
                !< New deviatoric stress tensor
                p(i)   = -(signxx(i) + signyy(i) + signzz(i))/three
                sxx(i) = signxx(i) + p(i)
                syy(i) = signyy(i) + p(i)
                szz(i) = signzz(i) + p(i)
                sxy(i) = signxy(i)
                syz(i) = signyz(i)
                szx(i) = signzx(i)
!
                !< New Von Mises stress
                seq(i) = three_half*(sxx(i)**2 + syy(i)**2 + szz(i)**2)        &
                       +      three*(sxy(i)**2 + syz(i)**2 + szx(i)**2)
                seq(i) = sqrt(seq(i))
!
              enddo 
!
              !< 9 - Update yield function value
              !-----------------------------------------------------------------
              !< Update material cohesion
              if (ifunc(3) > 0) then 
                ipos(1:nel) = vartmp(1:nel,3)
                iad(1:nel)  = npf(ifunc(3)) / 2 + 1
                ilen(1:nel) = npf(ifunc(3)+1) / 2 - iad(1:nel) - ipos(1:nel)
                call vinter2(tf,iad,ipos,ilen,nel,epspd,dcdepspd,c)
                c(1:nel) = cini*c(1:nel)
                dcdepspd(1:nel) = cini*dcdepspd(1:nel)
                vartmp(1:nel,3) = ipos(1:nel)
              endif
!
              !< Update cap limit pressure
              if (ifunc(4) > 0) then
                ipos(1:nel) = vartmp(1:nel,4)
                iad(1:nel)  = npf(ifunc(4)) / 2 + 1
                ilen(1:nel) = npf(ifunc(4)+1) / 2 - iad(1:nel) - ipos(1:nel)
                call vinter2(tf,iad,ipos,ilen,nel,epspv,dpbdepspv,pb)
                pb(1:nel) = capini*pb(1:nel)
                dpbdepspv(1:nel) = capini*dpbdepspv(1:nel)
                vartmp(1:nel,4) = ipos(1:nel)
                pa(1:nel) = alpha*pb(1:nel)
                dpadepspv(1:nel) = alpha*dpbdepspv(1:nel)  
              endif  
!     
              ! if (ifunc(1) > 0) k(i) = kini*finter(ifunc(1),epspv(i),npf,tf,deri)
              ! if (ifunc(2) > 0) g(i) = gini*finter(ifunc(2),epspv(i),npf,tf,deri)
!
              !< Loop over yielding elements
              do ii = 1,nindx
                i = indx(ii)
!
                !< Update null criterion derivative pressure (dfdp = 0)
                delta = (pa(i)*tgphi+c(i))**2+eight*((pb(i)-pa(i))**2)*(tgphi**2)
                if (tgphi > zero) then
                  p0(i) = pa(i)+(-(pa(i)*tgphi+c(i))+sqrt(delta))/(four*tgphi)
                else
                  p0(i) = pa(i)
                endif
!
                !< Something about porosity
                ! if (sat0 > zero) then
                !   por(i) = one - (one-por0)*exp(epspv(i)-epspvol0)
                !   fac = max(em03,por(i)/por0)
                !   muw(i) = sat0/fac*amu(i)
                !   if (muw(i) >= tolmu) then
                !     u(i) = kwater*muw(i)
                !   elseif (muw(i) > -tolmu) then
                !     u(i)=kwater/four/tolmu*(muw(i)+tolmu)**2
                !   else
                !     u(i)=zero
                !   endif
                !   if (muw(i) >= -tolmu) dudmu(i)=max(dudmu(i),kwater*sat0/fac)
                ! endif
!
                !< Transition yield surface to cap hardening factor
                if (p(i) <= pa(i)) then
                  rc(i) = one
                else
                  rc(i) = one - ((p(i)-pa(i))/(pb(i)-pa(i)))**2
                  rc(i) = sqrt(max(rc(i),zero))
                endif
                !< Pressure dependent factor
                a(i) = max(zero,p(i)*tgphi + c(i))
                !< New yield function value
                f(i) = seq(i) - rc(i)*a(i)
              enddo
            enddo      
          endif
          !=====================================================================
          ! - END OF PLASTIC CORRECTION WITH CUTTING PLANE METHOD
          !=====================================================================
!
          !< Update porosity variable
          if (sat0 > zero) then
            do i=1,nel
              uvar(i,5) = u(i)   ! cap shift
              uvar(i,6) = por(i)
              uvar(i,5) = muw(i) + one
              !< pore pressure is calculated here
              if (muw(i) >  zero) then
                u(i) = kwater*muw(i)
              else
                u(i) = zero
              endif
            enddo
            do i=1,nel
              viscmax(i) = zero
              !< fp_poro adding viscosity close to saturation
              if (muw(i) > -tolmu) then 
                viscmax(i) = viscfac*(sqrt(kwater*rho(i))*volume(i)**third)
                u(i) = u(i) - viscmax(i)*(depsxx(i)+depsyy(i)+depszz(i))/dt1
              endif 
              !< fp_poro the pore pressure is stored in the viscous stress
              !< for practical reasons including compatibility with ale
              sigvxx(i) = -u(i)
              sigvyy(i) = -u(i)
              sigvzz(i) = -u(i)
              uvar(i,3) = u(i)
            enddo
          endif
!
          !< Update user variables and compute the sound speed
          do i = 1,nel
            uvar(i,1)  = epspd(i)        !< Deviatoric equivalent plastic strain
            uvar(i,2)  = epspv(i)        !< Volumetric plastic strain
            uvar(i,3)  = c(i)            !< Material cohesion
            uvar(i,4)  = pb(i)           !< Cap limit pressure 
            uvar(i,5)  = pa(i)           !< Transition pressure yield surface to cap
            uvar(i,6)  = p0(i)           !< Null criterion derivative pressure (dfdp = 0)
            uvar(i,7)  = u(i)            !< Pore water pressure
            uvar(i,8)  = por(i)          !< Porosity
            uvar(i,9)  = muw(i) + one    !< Saturation
            uvar(i,10) = u(i)            !< Cap shift
            defp(1:nel,1) = epspd(1:nel) !< Deviatoric Equivalent Plastic Strain  
            defp(1:nel,2) = epspv(1:nel) !< Volumetric Plastic Strain  
            !< Sound speed in the material
            soundsp(i) = sqrt((k(i) + four_over_3*g(i) + dudmu(i))/rho0(i))
          enddo 
!
        end subroutine sigeps81
      end module sigeps81_mod