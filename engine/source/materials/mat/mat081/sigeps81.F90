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
! \brief Johnson-Holmquist 1 material law /MAT/LAW126
! \details Material law based on Johnson-Holmquist version 1 theory. Dedicated to concrete application. 
! ======================================================================================================================
        subroutine sigeps81 (
     1     nel    ,nuvar  ,uvar    ,matparam,nfunc   ,ifunc  ,
     2     ngl    ,npf    ,tf      ,time    ,rho0    ,rho    ,
     3     volume ,amu    ,defp    ,soundsp ,viscmax ,
     4     epspxx ,epspyy ,epspzz  ,epspxy  ,epspyz  ,epspzx ,
     5     depsxx ,depsyy ,depszz  ,depsxy  ,depsyz  ,depszx ,
     6     sigoxx ,sigoyy ,sigozz  ,sigoxy  ,sigoyz  ,sigozx ,
     7     signxx ,signyy ,signzz  ,signxy  ,signyz  ,signzx ,
     8     sigvxx ,sigvyy ,sigvzz  )
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
#include  "units_c.inc"
!-------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-------------------------------------------------------------------------------
          integer, intent(in) :: nel                           !< number of elements in the group
          integer, intent(in) :: nuvar                         !< number of user variables
          my_real, dimension(nel,nuvar), intent(inout) :: uvar !< user variables
          type(matparam_struct_), intent(in) :: matparam       !< material parameters data
          integer, intent(in) :: nfunc                         !< number of functions
          integer, dimension(nfunc), intent(in) :: ifunc       !< function index
          integer, dimension(nel), intent(in) :: ngl           !< element user IDs index table
          integer, dimension(*), intent(in) :: npf(*)          !< number of values for each function
          my_real, dimension(*), intent(in) :: tf(*)           !< function parameters
          my_real, intent(in) :: time                          !< current time
          my_real, dimension(nel), intent(in) :: rho0          !< initial density
          my_real, dimension(nel), intent(in) :: rho           !< current density
          my_real, dimension(nel), intent(in) :: volume        !< element volume
          my_real, dimension(nel), intent(in) :: amu           !< element volumetric strain
          my_real, dimension(nel,2), intent(inout) :: defp     !< element plastic strain
          my_real, dimension(nel), intent(out) :: soundsp      !< sound speed
          my_real, dimension(nel), intent(out) :: viscmax      !< maximum viscosity
          my_real, dimension(nel), intent(in)  :: epspxx       !< strain rate component xx
          my_real, dimension(nel), intent(in)  :: epspyy       !< strain rate component yy
          my_real, dimension(nel), intent(in)  :: epspzz       !< strain rate component zz
          my_real, dimension(nel), intent(in)  :: epspxy       !< strain rate component xy
          my_real, dimension(nel), intent(in)  :: epspyz       !< strain rate component yz
          my_real, dimension(nel), intent(in)  :: epspzx       !< strain rate component zx
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
!-------------------------------------------------------------------------------
!  L o c a l   V a r i a b l e s
!-------------------------------------------------------------------------------
          integer :: i,nit,soft_flag,iter
          integer, parameter :: niter = 3 
          my_real, parameter :: small = 1.0d-10
          my_real, parameter :: tol   = 1.0d-20
          my_real :: alpha,tgb,tgp,max_dilat,kini,gini,yldini,capini,delta,    &
            drcdp,dfdp,dfdpb,dfdc,dgdp,hh,dlambda,dwpx2,qx2,depspv,depslv,     &
            dfds1,dfds2,dfds3,dfds4,dfds5,dfds6,depsp1,depsp2,depsp3,depsp4,
            depsp5,depsp6,depsl1,depsl2,depsl3,depsl4,depsl5,depsl6,
            sy1,sy2,sy3,sy4,sy5,sy6,qy,fy,py,uy,deri,fac,
            epspvol0,rhow0,kwater,por0,viscfac,sat0,u0,tolmu,pu,
            drcdpb,aa,rc,rp,pmpa,f1,f2,x,x1,x2,dx    
          my_real, dimension(nel) :: k,g,g2,c,pb,depsvol,dav,d1,d2,d3,         &
            ds1,ds2,ds3,ds4,ds5,ds6,dp,s1,s2,s3,s4,s5,s6,p,                    &
            dcdepsp,dpbdepsp,pa,pbmpa,pbmpa2,s1n,s2n,s3n,s4n,s5n,s6n,          &
            p0,pn,qn,q2,q,f,de1,de2,de3,de4,de5,de6,por,muw,u,du,dudmu,        &
            po,uo,so1,so2,so3,so4,so5,so6,fo
          my_real ,dimension(nel) :: epspd,epspv
!=============================================================================== 
!
      !=========================================================================
      ! - INITIALISATION OF COMPUTATION ON TIME STEP
      !=========================================================================
      !< Recovering integer model parameter
      soft_flag = matparam%iparam(1)   !< Softening flag
      !< Recovering real model paramter 
      kini      = matparam%bulk        !< Initial bulk modulus
      gini      = matparam%shear       !< Initial shear modulus
      tgphi     = matparam%uparam(1)   !< Friction angle
      tgpsi     = matparam%uparam(2)   !< Plastic flow potential angle
      alpha     = matparam%uparam(3)   !< Ratio Pa/Pb
      max_dilat = matparam%uparam(4)   !< Maximum dilatancy
      epspvol0  = matparam%uparam(5)   !< Initial volumetric plastic strain
      kwater    = matparam%uparam(6)   !< Pore water bulk modulus
      por0      = matparam%uparam(7)   !< Initial porosity
      sat0      = matparam%uparam(8)   !< Initial saturation
      u0        = matparam%uparam(9)   !< Initial pore water pressure
      tolmu     = matparam%uparam(10)  !< Tolerance for capshift viscosity
      viscfac   = matparam%uparam(11)  !< Viscosity factor
      yldini    = matparam%uparam(12)  !< Material cohesion
      capini    = matparam%uparam(13)  !< Initial cap limit pressure
!
      !========================================================================
      !< - Recovering user variables and state variables
      !========================================================================
      epspd(1:nel) = defp(1:nel,1)   !< Deviatoric Equivalent Plastic Strain  
      epspv(1:nel) = defp(1:nel,2)   !< Volumetric Plastic Strain  
!
      !========================================================================
      !< - Recovering elastic parameters
      !========================================================================
      !<  Bulk modulus computation
      !  -> Interpolated with volumetric plastic strain
      if (ifunc(1) > 0) then
        call vinter2(tf,iad,ipos,ilen,nel,epspv,deri,k)
        k(1:nel) = kini*k(1:nel)
      !  -> Constant
      else
        k(1:nel) = kini
      endif
      !<  Shear modulus computation
      !   -> Interpolated with volumetric plastic strain
      if (ifunc(2) > 0) then
        call vinter2(tf,iad,ipos,ilen,nel,epspv,deri,g)
        g(1:nel) = gini*g(1:nel)
      !   -> Constant
      else
        g(1:nel) = gini
      endif
      !< Two*shear modulus
      g2(1:nel) = two*g(1:nel)
      !< Lame coefficient
      lam(1:nel) = k(1:nel) - two*g(1:nel)/three
!
      !========================================================================
      !< - Recovering yield criterion parameters
      !========================================================================
      !< Material cohesion
      !  -> Interpolated with deviatoric plastic strain
      if (ifunc(3) > 0) then 
        call vinter2(tf,iad,ipos,ilen,nel,epspd,dcdepsp,c)
        c(1:nel) = yldini*c(1:nel)
      !  -> Constant
      else
        c(1:nel) = yldini
      endif
      !< Cap limit pressure
      !  -> Interpolated with volumetric plastic strain
      if (ifunc(4) > 0) then
        call vinter2(tf,iad,ipos,ilen,nel,epspv,dpbdepsp,pb)
        pb(1:nel) = capini*pb(1:nel)
      !  -> Constant
      else
        pb(1:nel) = capini
      endif
      !< Transition pressures yield surface to cap
      pa(1:nel) = alpha*pb(1:nel)  
      do i = 1,nel
        delta = (pa(i)*tgb + c(i))**2 + eight*((pb(i)-pa(i))**2)*tgb**2
        if (tgb > small) then
          p0(i) = pa(i) + (-(pa(i)*tgb+c(i)) + sqrt(delta))/four/tgb
        else
          p0(i) = pa(i)
        endif
      enddo
!      
      !< User variables initilization     
      if ((isigi == 0).and.(time == zero)) then  
        do i=1,nel
          uvar(i,3) = c(i)
          uvar(i,4) = pb(i)
          uvar(i,5) = u0
          uvar(i,6) = por0
          uvar(i,7) = sat0
          muw(i)    = sat0-one
          if (muw(i) >= tolmu) then
            uvar(i,6) = kwater*muw(i)
          elseif (muw(i) > -tolmu) then
            uvar(i,6) = kwater/four/tolmu*(muw(i)+tolmu)**2
          else
            uvar(i,6) = zero
          endif
        enddo
      endif
! 
      !========================================================================
      !< - Computation of trial stress tensor, Von Mises and pressure
      !========================================================================
      do i=1,nel
!
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
        sxx(i) = signxx(i) - (signxx(i) + signyy(i) + signzz(i))/three
        syy(i) = signyy(i) - (signxx(i) + signyy(i) + signzz(i))/three
        szz(i) = signzz(i) - (signxx(i) + signyy(i) + signzz(i))/three
        sxy(i) = signxy(i)
        syz(i) = signyz(i)
        szx(i) = signzx(i)
!
        !< Trial Von Mises stress
        seq(i) = three_half*(sxx(i)**2 + syy(i)**2 + szz(i)**2)                &
               +      three*(sxy(i)**2 + syz(i)**2 + szx(i)**2)
        seq(i) = sqrt(seq(i))
      enddo
!
      !< Porosity computation
      if (sat0 > zero) then
        do i=1,nel
          por(i) = one - (one-por0)*exp(epspv(i) - epspvol0)
          fac    = max(em03,por(i)/por0)
          muw(i) = sat0/fac*amu(i)
          if (muw(i) >= tolmu) then
            u(i) = kwater*muw(i)
          elseif (muw(i) > -tolmu) then
            u(i) = kwater/four/tolmu*(muw(i)+tolmu)**2
          else
            u(i) = zero
          endif
          du(i) = u(i) - uvar(i,6)
          dudmu(i) = zero
          if (muw(i) >= -tolmu) dudmu(i) = kwater*sat0/fac
        enddo
      else
        muw(1:nel)   = -one
        u(1:nel)     = zero
        du(1:nel)    = zero
        dudmu(1:nel) = zero
      endif  
!
      !========================================================================
      ! - Computation of yield function and check element behavior
      !========================================================================
      nindx = 0
      indx(1:nel) = 0
      do i=1,nel          
        ! if (p(i) < p0(i)) then
        !   pu(i) = p(i)
        ! elseif (p(i) <= p0(i) + u(i)) then
        !   pu(i) = p0(i)
        ! else
        !   pu(i) = p(i) - u(i)
        ! endif
        dp(i) = p(i) - pa(i)
        if (dp(i) > zero) then
          rc(i) = one - min(one, dp(i)**2/((pb(i)-pa(i))**2))
          rc(i) = sqrt(rc(i))
        else
          rc(i) = one
        endif
        a(i) = max(zero, p(i)*tgb + c(i))
        f(i) = seq(i) - rc(i)*a(i)
        if (f(i) > zero) then 
          nindx = nindx + 1
          indx(nindx) = i
        enddo
      enddo
!
      !========================================================================
      ! - PLASTIC CORRECTION WITH CUTTING PLANE (NEWTON-ITERATION) METHOD
      !======================================================================== 
      if (nindx > 0) then
!
        !< Loop over the iterations   
        do iter = 1, niter
          !< Loop over yielding elements
          do ii = 1, nindx
            i = indx(ii)
!       
            ! Note: in this part, the purpose is to compute for each iteration
            ! a plastic multiplier allowing to update internal variables to satisfy
            ! the consistency condition using the cutting plane method within an
            ! iterative procedure.
            ! Its expression at each iteration is : dlambda = - f/df_dlambda
            ! -> f          : current value of yield function (known)
            ! -> df_dlam    : derivative of f with respect to dlambda by taking
            !                 into account of internal variables kinetic : 
            !                 plasticity, damage ... (to compute)
!      
            ! 1 - Derivative of yield criterion w.r.t plastic multiplier
            !     Contribution of the stress tensor
            !------------------------------------------------------------------
            !< Derivative of Von Mises stress w.r.t stress tensor
            dseqdsigxx = three_half*sxx(i)/max(seq(i),em20)
            dseqdsigyy = three_half*syy(i)/max(seq(i),em20)
            dseqdsigzz = three_half*szz(i)/max(seq(i),em20)
            dseqdsigxy = three*sxy(i)/max(seq(i),em20)
            dseqdsigyz = three*syz(i)/max(seq(i),em20)
            dseqdsigzx = three*szx(i)/max(seq(i),em20)
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
            if (p(i) <= pa(i)) then
              dfdp  = -tgphi
            else
              drcdp = (pa(i) - p(i))/(rc(i)*(pb(i) - pa(i))**2)
              dfdp  = -(drcdp*a(i) + rc(i)*tgphi)
            endif
            dfdseq = one
!
            !< Derivative of plastic potential with respect to pressure and 
            !  Von Mises stress
            if (p(i) <= pa(i)) then
              dgdp = -tgpsi
            elseif (p(i) <= p0(i)) then
              dgdp = -tgpsi*(one - ((p(i) - pa(i))/(p0(i) - pa(i))))
            elseif (p(i) > po(i)) then 
              dgdp = dfdp
            endif  
            dgdseq = one
!
            !< Assembling derivative of yield criterion w.r.t stress tensor
            dfdsigxx = dfdseq*dseqdsigxx + dfdp*dpdsigxx
            dfdsigyy = dfdseq*dseqdsigyy + dfdp*dpdsigyy
            dfdsigzz = dfdseq*dseqdsigzz + dfdp*dpdsigzz
            dfdsigxy = dfdseq*dseqdsigxy + dfdp*dpdsigxy
            dfdsigyz = dfdseq*dseqdsigyz + dfdp*dpdsigyz
            dfdsigzx = dfdseq*dseqdsigzx + dfdp*dpdsigzx
!
            !<Assembling derivatives of plastic potential w.r.t stress tensor
            dgdsigxx = dgdseq*dseqdsigxx + dgdp*dpdsigxx
            dgdsigyy = dgdseq*dseqdsigyy + dgdp*dpdsigyy
            dgdsigzz = dgdseq*dseqdsigzz + dgdp*dpdsigzz
            dgdsigxy = dgdseq*dseqdsigxy + dgdp*dpdsigxy
            dgdsigyz = dgdseq*dseqdsigyz + dgdp*dpdsigyz
            dgdsigzx = dgdseq*dseqdsigzx + dgdp*dpdsigzx
!
            !< Derivative of stress tensor w.r.t plastic multiplier
            trdgpds    = dgdsigxx + dgdsigyy + dgdsigzz  
            dsigxxdlam = -(dgdsigxx*g2(i) + lam(i)*trdgpds)
            dsigyydlam = -(dgdsigyy*g2(i) + lam(i)*trdgpds)
            dsigzzdlam = -(dgdsigzz*g2(i) + lam(i)*trdgpds) 
            dsigxydlam =  - dgdsigxy*g(i)
            dsigyzdlam =  - dgdsigyz*g(i)
            dsigzxdlam =  - dgdsigzx*g(i)   
!
            !< Contribution of the stress tensor to the derivative of 
            !< the yield criterion with respect to the plastic multiplier
            dfdsig_dsigdlam = dfdsigxx*dsigxxdlam + dfdsigyy*dsigyydlam +      &  
                              dfdsigzz*dsigzzdlam + dfdsigxy*dsigxydlam +      &
                              dfdsigyz*dsigyzdlam + dfdsigzx*dsigzxdlam
!
            ! 2 - Derivative of yield criterion w.r.t pressures Pa and Pb
            !-------------------------------------------------------------------
            dfdrc = -one
            if (p(i) <= pa(i)) then
              drcdpb = zero
              drcdpa = zero
            else
              drcdpb = ((p(i) - pa(i))**2)/(rc(i)*(pb(i) - pa(i))**3)
              drcdpa = -((p(i)-pa(i))*(p(i)-pb(i)))/(rc(i)*(pb(i)-pa(i))**3)
            endif
!
            ! 3 - Derivative of yield criterion w.r.t material cohesion
            !-------------------------------------------------------------------
            dfdc = -rc(i)
!
            ! 4 - Derivative of deviatoric and volumetric plastic strain w.r.t
            !     plastic multiplier
            !-------------------------------------------------------------------
            depspd_dlam = (signxx(i)*dgdsigxx(i) + signyy(i)*dgdsigyy(i) +     &
                           signzz(i)*dgdsigzz(i) + signxy(i)*dgdsigxy(i) +     &
                           signyz(i)*dgdsigyz(i) + signzx(i)*dgdsigzx(i))      &
                           /max(seq(i),em20)
            depspv_dlam = dgdsixx + dgdsiyy + dgdsizz
!
            ! 4 - Derivative of yield criterion w.r.t plastic multiplier
            !------------------------------------------------------------------- 
            df_dlam = -dfdsig_dsigdlam +                                       &
                       dfdrc*(drcdpa*alpha+drcdpb)*dpbdepspv(i)*depspv_dlam +  &
                       dfdc*dcdepspd(i)*depspd_dlam
            df_dlam = sign(max(abs(df_dlam),em20),df_dlam)
!
            ! 5 - Computation of plastic multiplier
            !-------------------------------------------------------------------               
            dlam = -f(i)/df_dlam
!
            ! 6 - Update plastic strain related variables
            !------------------------------------------------------------------- 
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
            depspd(i) = depspv(i) + depspd_dlam*dlam
            epspd(i)  = epspd(i)  + depspd(i)
! 
            ! 6 - Update stress tensor and yield function
            !-------------------------------------------------------------------
            trdep     = dpxx(i) + dpyy(i) + dpzz(i)
            signxx(i) = signxx(i) - (g2(i)*dpxx(i) + lam(i)*trdep)
            signyy(i) = signyy(i) - (g2(i)*dpyy(i) + lam(i)*trdep)
            signzz(i) = signzz(i) - (g2(i)*dpzz(i) + lam(i)*trdep)
            signxy(i) = signxy(i) - ( g(i)*dpxy(i))
            signyz(i) = signyz(i) - ( g(i)*dpyz(i))
            signzx(i) = signzx(i) - ( g(i)*dpzx(i))
!
            !< New deviatoric stress tensor
            sxx(i) = signxx(i) - (signxx(i) + signyy(i) + signzz(i))/three
            syy(i) = signyy(i) - (signxx(i) + signyy(i) + signzz(i))/three
            szz(i) = signzz(i) - (signxx(i) + signyy(i) + signzz(i))/three
            sxy(i) = signxy(i)
            syz(i) = signyz(i)
            szx(i) = signzx(i)
!
            !< New Von Mises stress
            seq(i) = three_half*(sxx(i)**2 + syy(i)**2 + szz(i)**2)            &
                   +      three*(sxy(i)**2 + syz(i)**2 + szx(i)**2)
            seq(i) = sqrt(seq(i))
! 
          enddo 
!
          !< Update material cohesion
          if (ifunc(3) > 0) then 
            ilen = 
            ipos = 
            iad  = 
            call vinter2(tf,iad,ipos,ilen,nel,epspd,dcdepsp,c)
            c(1:nel) = yldini*c(1:nel)
            dcdepspd(1:nel) = yldini*dcdepspd(1:nel)
          endif
!
          !< Cap limit pressure
          if (ifunc(4) > 0) then
            ilen 
            ipos 
            iad 
            call vinter2(tf,iad,ipos,ilen,nel,epspv,dpbdepsp,pb)
            pb(1:nel) = capini*pb(1:nel)
            pa(1:nel) =  alpha*pb(1:nel)
            dpbdepsp(1:nel) = capini*dpbdepsp(1:nel)
          endif  
!     
          ! if (ifunc(1) > 0) k(i) = kini*finter(ifunc(1),epspv(i),npf,tf,deri)
          ! if (ifunc(2) > 0) g(i) = gini*finter(ifunc(2),epspv(i),npf,tf,deri)
!
          !< Loop over yielding elements
          do ii = 1, nindx
            i = indx(ii)
!
            !< Pressure
            delta = (pa(i)*tgphi+c(i))**2+eight*((pb(i)-pa(i))**2)*tgphi**2
            if (tgphi > small) then
              p0(i) = pa(i)+(-(pa(i)*tgphi+c(i))+sqrt(delta))/four/tgphi
            else
              p0(i) = pa(i)
            endif

            if (sat0 > zero) then
              por(i) = one - (one-por0)*exp(epspv(i)-epspvol0)
              fac = max(em03,por(i)/por0)
              muw(i) = sat0/fac*amu(i)
              if (muw(i) >= tolmu) then
                u(i) = kwater*muw(i)
              elseif (muw(i) > -tolmu) then
                u(i)=kwater/four/tolmu*(muw(i)+tolmu)**2
              else
                u(i)=zero
              endif
              if (muw(i) >= -tolmu) dudmu(i)=max(dudmu(i),kwater*sat0/fac)
            endif

            fy = fcrit(pn(i),u(i),tgb,c(i),pa(i),p0(i),pbmpa2(i),qn(i))

            if (fy > small) then
              if (tgb*pn(i) + c(i) <= zero) then
                pn(i)  =-c(i)/tgb
                s1n(i) = zero
                s2n(i) = zero
                s3n(i) = zero
                s4n(i) = zero
                s5n(i) = zero
                s6n(i) = zero
                write(6,*)'tri-traction failure'
              else 
                if (pn(i) < p0(i)) then
                  pu = pn(i)
                else if (pn(i) <= p0(i)+u(i)) then
                  pu = p0(i)
                else
                  pu = pn(i) - u(i)
                endif
                if (pn(i) > pa(i)) then
                  rc = sqrt(max(zero,one - ((pu-pa(i))**2/pbmpa2(i))))
                else
                  rc = one
                endif
                if (rc > tol ) then        ! changed by marian at 2017.12.06
                  if (qn(i) > small) then
                    x = rc*(pu*tgb+c(i))/qn(i)
                    if (x < one-em02 .or. x > one) then
                      write(7,*)'reproj q',x,fy,pn(i),pu,qn(i)
                      write(6,*)'reproj q',x,fy,pn(i),pu,qn(i)
                    endif
                    s1n(i) = x*s1n(i)
                    s2n(i) = x*s2n(i)
                    s3n(i) = x*s3n(i)
                    s4n(i) = x*s4n(i)
                    s5n(i) = x*s5n(i)
                    s6n(i) = x*s6n(i)
                  endif
                else
                  s1n(i) = zero
                  s2n(i) = zero 
                  s3n(i) = zero
                  s4n(i) = zero 
                  s5n(i) = zero 
                  s6n(i) = zero
                  pn(i)  = pb(i)+u(i)                    
                endif               
              endif
            endif
            f(i) = fcrit(p(i),u(i),tgb,c(i),pa(i),p0(i),pbmpa2(i),q(i))
          enddo
        enddo      
      endif
      !========================================================================
      ! - END OF PLASTIC CORRECTION WITH CUTTING PLANE METHOD
      !======================================================================== 
!
      !< Update user variables
      do i = 1,nel
        uvar(i,1)  = epspd(i)
        uvar(i,2)  = epspv(i)
        uvar(i,3)  = c(i)
        uvar(i,4)  = pb(i)
        soundsp(i) = sqrt((k(i) + four_over_3*g(i) + dudmu(i))/rho(i))
      enddo 
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
            u(i) = u(i) - viscmax(i)*(epspxx(i)+epspyy(i)+epspzz(i))
          endif 
          !< fp_poro the pore pressure is stored in the viscous stress
          !< for practical reasons including compatibility with ale
          sigvxx(i) = -u(i)
          sigvyy(i) = -u(i)
          sigvzz(i) = -u(i)
          uvar(i,3) = u(i)
        enddo
      endif
    end subroutine sigeps81
  end module sigeps81_mod

