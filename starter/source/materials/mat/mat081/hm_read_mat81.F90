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
      module hm_read_mat81_mod
      contains
! ======================================================================================================================
! \brief Reading material parameters of /MAT/LAW81
! \details Reading material parameters of /MAT/LAW81
! ======================================================================================================================
      subroutine hm_read_mat81(matparam ,nuvar    ,ifunc    ,maxfunc  ,        &
                               nfunc    ,parmat   ,mat_id   ,titr     ,        &
                               unitab   ,lsubmodel,mtag     ,iout     )
!-------------------------------------------------------------------------------
!   M o d u l e s
!-------------------------------------------------------------------------------
      use unitab_mod
      use message_mod
      use submodel_mod
      use matparam_def_mod    
      use elbuftag_mod      
      use constant_mod   
!-------------------------------------------------------------------------------
!   I m p l i c i t   T y p e s
!-------------------------------------------------------------------------------
      implicit none 
#include  "my_real.inc"
!-------------------------------------------------------------------------------
!   D u m m y   A r g u m e n t s
!-------------------------------------------------------------------------------
      type(matparam_struct_) ,intent(inout)      :: matparam  !< matparam data structure
      integer, intent(inout)                     :: nuvar     !< number of material law user variables
      integer, dimension(maxfunc), intent(inout) :: ifunc     !< ids. of functions associated to the material
      integer, intent(in)                        :: maxfunc   !< maximal size of ifunc
      integer, intent(inout)                     :: nfunc     !< number of functions associated to the material
      my_real, dimension(100),intent(inout)      :: parmat    !< material parameter global table 1
      integer, intent(in)                        :: mat_id    !< material law user ID 
      character(len=nchartitle),intent(in)       :: titr      !< material law user title
      type(unit_type_),intent(in)                :: unitab    !< units table
      type(submodel_data),dimension(nsubmod),intent(in) :: lsubmodel !< submodel data structure
      type(mlaw_tag_), intent(inout)             :: mtag      !< material tag for internal variables in element buffer
      integer, intent(in)                        :: iout      !< output file number
!------------------------------------------------------------------------------
!   L o c a l   V a r i a b l e s
!------------------------------------------------------------------------------
      my_real :: kini,gini,yldini,capini,beta,psi,alpha,max_dilat,epsvini,     &
        kwater,por0,sat0,u0,muw0,viscfac,tol,rho0,rhor,fac_unit
      integer :: soft_flag,ilaw
      logical :: is_available,is_encrypted
!-------------------------------------------------------------------------------
!     S o u r c e 
!-------------------------------------------------------------------------------    
      is_encrypted = .false.
      is_available = .false.
!
      ilaw = 81
!-------------------------------------------------------------------------------
      call hm_option_is_encrypted(is_encrypted)
!-------------------------------------------------------------------------------
!Card1
      call hm_get_floatv('MAT_RHO'  ,rho0     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('REFER_RHO',rhor     ,is_available, lsubmodel, unitab)
!Card2
      call hm_get_floatv('K0'       ,kini     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_G0'   ,gini     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_COH0' ,yldini   ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_PB0'  ,capini   ,is_available, lsubmodel, unitab)
!Card3
      call hm_get_floatv('MAT_BETA' ,beta     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('PSI'      ,psi      ,is_available, lsubmodel, unitab)
!Card4
      call hm_get_floatv('MAT_ALPHA',alpha    ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_EPS'  ,max_dilat,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_SRP'  ,epsvini  ,is_available, lsubmodel, unitab)
!Card5
      call hm_get_intv  ('FUN_A1'   ,ifunc(1) ,is_available, lsubmodel)
      call hm_get_intv  ('FUN_A2'   ,ifunc(2) ,is_available, lsubmodel)
      call hm_get_intv  ('FUN_A3'   ,ifunc(3) ,is_available, lsubmodel)
      call hm_get_intv  ('FUN_A4'   ,ifunc(4) ,is_available, lsubmodel)
      call hm_get_intv  ('IFLAG'    ,soft_flag,is_available, lsubmodel)
!Card6
      call hm_get_floatv('MAT_KW'   ,kwater   ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_POR0' ,por0     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_SAT0' ,sat0     ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_MUE0' ,u0       ,is_available, lsubmodel, unitab)
!Card7
      call hm_get_floatv('MAT_TOL'  ,tol      ,is_available, lsubmodel, unitab)
      call hm_get_floatv('MAT_VIS'  ,viscfac  ,is_available, lsubmodel, unitab)
!      
!-------------------------------------------------------------------------------
!     Default values
!-------------------------------------------------------------------------------
      if (kwater == zero) kwater = one
      if (por0 == zero) por0 = zero
      if (sat0 == zero) sat0 = zero
      if (u0 == zero) u0 = zero
      if (tol == zero) tol = em04
      if (viscfac == zero) viscfac = half
      if (max_dilat == zero) max_dilat = -infinity
      max_dilat = -abs(max_dilat)
      if (alpha == zero) alpha = half
!-------------------------------------------------------------------------------
!     Data checking
!-------------------------------------------------------------------------------
      if (kini<=zero) then
        call ancmsg(msgid=1012,                                                &
                    msgtype=msgerror,                                          &
                    anmode=aninfo_blind_1,                                     &
                    i1=mat_id,                                                 &
                    c1=titr)
      endif
      if (gini<=zero) then
        call ancmsg(msgid=1013,                                                &
                    msgtype=msgerror,                                          &
                    anmode=aninfo_blind_1,                                     &
                    i1=mat_id,                                                 &
                    c1=titr)
      endif
      if (ifunc(1) == 0) then
        call ancmsg(msgid=1014,                                                &
                    msgtype=msgwarning,                                        &
                    anmode=aninfo,                                             &
                    i1=mat_id,                                                 &
                    c1=titr)
      endif
      if (ifunc(2) == 0) then
        call ancmsg(msgid=1015,                                                &
                    msgtype=msgwarning,                                        &
                    anmode=aninfo,                                             &
                    i1=mat_id,                                                 &
                    c1=titr)
      endif
      if (ifunc(3) == 0) then
        if (yldini == zero) then 
          call ancmsg(msgid=1713,                                              &
                      msgtype=msgerror,                                        &
                      anmode=aninfo_blind_1,                                   &
                      i1=mat_id,                                               &
                      c1=titr)         
        else
          call ancmsg(msgid=1016,                                              &
                      msgtype=msgwarning,                                      &
                      anmode=aninfo,                                           &
                      i1=mat_id,                                               &
                      c1=titr)
        endif
      endif
      if (ifunc(4) == 0) then
        if (capini == zero) then 
          call ancmsg(msgid=1714,                                             &
                      msgtype=msgerror,                                       &
                      anmode=aninfo_blind_1,                                  &
                      i1=mat_id,                                              &
                      c1=titr)         
        else
          call ancmsg(msgid=1017,                                             &
                      msgtype=msgwarning,                                     &
                      anmode=aninfo,                                          &
                      i1=mat_id,                                              &
                      c1=titr)        
        endif
      endif
      ! Default value for scale factors of c and Pb
      if (yldini == zero) then 
        call hm_get_floatv_dim('mat_coh0',fac_unit,is_available,lsubmodel,    &
                               unitab    )
        yldini = one*fac_unit
      endif
      if (capini == zero) then 
        call hm_get_floatv_dim('mat_pb0',fac_unit,is_available,lsubmodel,     &
                              unitab    )
        capini = one*fac_unit
      endif
      if (sat0 /= zero) then
        if (kwater <= zero) then
          call ancmsg(msgid=1085,                                             &
                      msgtype=msgerror,                                       &
                      anmode=aninfo_blind_1,                                  &
                      i1=mat_id,                                              &
                      c1=titr)
        endif
        muw0 = 0
        if (por0 == zero) then
          call ancmsg(msgid=1086,                                             &
                      msgtype=msgerror,                                       &
                      anmode=aninfo_blind_1,                                  &
                      i1=mat_id,                                              &
                      c1=titr)
        elseif (u0 > zero) then
          muw0 = u0/kwater
          sat0 = one + muw0
        else
          muw0 = sat0 - one
          if (muw0 > zero) then
            u0 = kwater*muw0
          else 
            u0 = zero
          endif
        endif
      else
        muw0 = -one
      endif
!
!-------------------------------------------------------------------------------
!     Filling buffer tables
!------------------------------------------------------------------------------- 
      !< Number of integer material parameters
      matparam%niparam = 1
      !< Number of real material parameters
      matparam%nuparam = 15
      !< Number of user variables 
      nuvar = 6
      !< Number of functions
      nfunc = 4
!          
      !< Allocation of material parameters tables
      allocate (matparam%iparam(matparam%niparam))
      allocate (matparam%uparam(matparam%nuparam))
!     
      !< Integer material parameter
      matparam%iparam(1)  = soft_flag
!
      !< Real material parameters
      matparam%uparam(1)  = tand(beta)
      matparam%uparam(2)  = tand(psi)
      matparam%uparam(3)  = alpha
      matparam%uparam(4)  = max_dilat
      matparam%uparam(5)  = epsvini
      matparam%uparam(6)  = kwater
      matparam%uparam(7)  = por0
      matparam%uparam(8)  = muw0 + one
      matparam%uparam(9)  = u0
      matparam%uparam(10) = tol
      matparam%uparam(11) = viscfac
      matparam%uparam(12) = yldini
      matparam%uparam(13) = capini
!
      !< Elastic parameters
      matparam%bulk  = kini
      matparam%shear = gini
      matparam%young = nine*kini*gini/(three*kini+gini)
      matparam%nu    = (three*kini-two*gini)/(six*kini+two*gini)
!
      !< PARMAT table
      parmat(1)  = kini + four_over_3*gini
      parmat(16) = 2
      parmat(17) = two*gini/(kini + four_over_3*gini) 
!
      !< PM table
      if (rhor == zero) rhor = rho0
      matparam%rho  = rhor
      matparam%rho0 = rho0
!
      !< MTAG variable activation
      mtag%g_seq = 1
      mtag%l_seq = 1
      mtag%g_pla = 2
      mtag%l_pla = 2
!
      !< Properties compatibility  
      call init_mat_keyword(matparam,"SOLID_ISOTROPIC") 
      call init_mat_keyword(matparam,"SPH")      
! 
      !< MATPARAM keywords   
      call init_mat_keyword(matparam,"HOOK")
      call init_mat_keyword(matparam,"COMPRESSIBLE")
      call init_mat_keyword(matparam,"INCREMENTAL" )
      call init_mat_keyword(matparam,"LARGE_STRAIN")
      call init_mat_keyword(matparam,"HYDRO_EOS") 
      call init_mat_keyword(matparam,"ISOTROPIC") 
!
!-------------------------------------------------------------------------------
!     Parameters printout
!-------------------------------------------------------------------------------
      write(iout,1000) trim(titr),mat_id,ilaw
      write(iout,1100)
      if (is_encrypted) then
        write(iout,'(5x,a,//)')'CONFIDENTIAL DATA'
      else
        write(iout,1200) rho0
        write(iout,1300) kini,gini,yldini,capini
        write(iout,1400) beta,psi
        write(iout,1500) alpha,max_dilat,epsvini
        write(iout,1600) ifunc(1),ifunc(2),ifunc(3),ifunc(4),soft_flag
        write(iout,1700) kwater,por0,sat0,u0
        write(iout,1800) tol,viscfac
      endif     
!-------------------------------------------------------------------------------
 1000 format(/                                                                 &
       5x,a,/,                                                                 &
       5x,'MATERIAL NUMBER. . . . . . . . . . . . . . .=',i10/,                &
       5x,'MATERIAL LAW . . . . . . . . . . . . . . . .=',i10/)
 1100 format(                                                                  &
       5x,'---------------------------------------------',/,                   &
       5x,'  MATERIAL MODEL: DRUCKER PRAGER WITH CAP    ',/,                   &
       5x,'---------------------------------------------',/)
 1200 FORMAT(                                                                  &
       5x,'INITIAL DENSITY  . . . . . . . . . . . . . .=',1pg20.13/)  
 1300 FORMAT(                                                                  &
       5x,'INITIAL BULK MODULUS (K0) . . . . . . . . . =',1pg20.13/,           &
       5x,'INITIAL SHEAR MODULUS (G0). . . . . . . . . =',1pg20.13/,           &
       5x,'INITIAL MATERIAL COHESION (C0). . . . . . . =',1pg20.13/,           &
       5x,'INITIAL CAP LIMIT PRESSURE (PB0). . . . . . =',1pg20.13/)
 1400 FORMAT(                                                                  &
       5x,'FRICTION ANGLE (PHI). . . . . . . . . . . . =',1pg20.13/,           &
       5x,'PLASTIC FLOW ANGLE (PSI). . . . . . . . . . =',1pg20.13/)           
 1500 FORMAT(                                                                  &
       5x,'ALPHA: RATIO PA/PB. . . . . . . . . . . . . =',1pg20.13/,           &
       5x,'MAXIMUM DILATANCY (EPS_MAX) . . . . . . . . =',1pg20.13/,           &
       5x,'INITIAL PLASTIC VOLUMETRIC STRAIN (EPSPV0). =',1pg20.13/)
 1600 FORMAT(                                                                  &
       5x,'BULK MODULUS K SCALE FUNCTION . . . . . . . =',i10/,                &
       5x,'SHEAR MODULUS G SCALE FUNCTION  . . . . . . =',i10/,                &
       5x,'MATERIAL COHESION C SCALE FUNCTION  . . . . =',i10/,                &
       5x,'CAP LIMIT PRESSURE PB SCALE FUNCTION  . . . =',i10/,                &
       5x,'FLAG: =0 SOFTENING ALLOWED, =1 NO SOFTENING =',i10/)
 1700 FORMAT(                                                                  &
       5x,'BULK MODULUS OF WATER . . . . . . . . . . . =',1pg20.13/            &
       5x,'INITIAL POROSITY POR0 . . . . . . . . . . . =',1pg20.13/            &
       5x,'INITIAL WATER SATURATION SAT0 . . . . . . . =',1pg20.13/            &
       5x,'INITIAL PORE PRESSURE U0  . . . . . . . . . =',1pg20.13/)
 1800 FORMAT(                                                                  &
     & 5x,'TOLERANCE FOR THE CRITERION SHIFT . . . . . =',1pg20.13/            &
     & 5x,'VISCOSITY FACTOR  . . . . . . . . . . . . . =',1pg20.13/)
!-------------------------------------------------------------------------------
      end subroutine hm_read_mat81
      end module hm_read_mat81_mod

