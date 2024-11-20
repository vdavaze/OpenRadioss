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
      module hm_read_mat129_mod
      contains
! ======================================================================================================================
! \brief Reading material parameters of /MAT/LAW129
! \details Reading material parameters of /MAT/LAW129
! ======================================================================================================================
      subroutine hm_read_mat129(                                               &       
                   nuvar    ,npropm   ,iout     ,mtag     ,parmat   ,unitab   ,&
                   pm       ,lsubmodel,israte   ,mat_id   ,titr     ,matparam ,&
                   nvartmp  )                  
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use unitab_mod
      use message_mod
      use submodel_mod
      use matparam_def_mod    
      use elbuftag_mod      
      use constant_mod    
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
      implicit none 
#include  "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)                          :: mat_id,npropm,iout
      integer, intent(inout)                       :: nuvar,nvartmp
      integer, intent(inout)                       :: israte
      type(mlaw_tag_), intent(inout)               :: mtag
      my_real, dimension(100),intent(inout)        :: parmat
      type (unit_type_),intent(in)                 :: unitab 
      my_real, dimension(npropm) ,intent(inout)    :: pm   
      type(submodel_data), dimension(*),intent(in) :: lsubmodel  
      character(len=nchartitle),intent(in)         :: titr 
      type(matparam_struct_) ,intent(inout)        :: matparam  
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: tab_lch00,tab_lch45,tab_lch90,tab_lchbi,tab_lchsh,ilaw
      integer :: tab_lcr00,tab_lcr45,tab_lcr90,tab_lcrbi,tab_lcrsh
      my_real :: hosf,mexp,nu,rho0,young,a11,a12,bulk,asrate
      logical :: is_available,is_encrypted
!=======================================================================
      is_encrypted = .false.
      is_available = .false.
      ilaw = 129
!------------------------------------------
      call hm_option_is_encrypted(is_encrypted)
!------------------------------------------
! - Density
      call hm_get_floatv('Rho'      ,rho0     ,is_available,lsubmodel, unitab)
! - Elastic properties flags and yield criterion parameters
      call hm_get_floatv('E'        ,young    ,is_available,lsubmodel, unitab)
      call hm_get_floatv('Nu'       ,nu       ,is_available,lsubmodel, unitab)
      call hm_get_floatv('LSD_LCHSH',hosf     ,is_available,lsubmodel, unitab)
      call hm_get_floatv('LSD_MAT_M',mexp     ,is_available,lsubmodel, unitab)
! - Uniaxial loading curves
      call hm_get_intv  ('LSD_LCH00',tab_lch00,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCH45',tab_lch45,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCH90',tab_lch90,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCHBI',tab_lchbi,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCHSH',tab_lchsh,is_available,lsubmodel)      
! - Lankford coefficients curves
      call hm_get_intv  ('LSD_LCR00',tab_lcr00,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCR45',tab_lcr45,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCR90',tab_lcr90,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCRBI',tab_lcrbi,is_available,lsubmodel)      
      call hm_get_intv  ('LSD_LCRSH',tab_lcrsh,is_available,lsubmodel) 
! 
!------------------------------------------------------------------------------
!     Default values and check parameters
!------------------------------------------------------------------------------
      ! Elastic properties
      bulk = young/(three*(one-two*nu))
      a11  = young/(one-nu*nu)
      a12  = nu*a11
      ! Default exponent
      if (mexp == zero) mexp = six
      ! Default strain rate filtering frequency
      if (israte == 0) then
        asrate = zero
      else
        asrate = 10000.0d0*unitab%FAC_T_WORK
      endif
      if (tab_lch00 == 0) then
        write(*,*) 'Error: LSD_LCH00 is not defined'
        stop
      endif
!------------------------------------------------------------------------------
!     Filling buffer tables
!------------------------------------------------------------------------------ 
      ! Number of integer material parameters
      matparam%niparam = 0
      ! Number of real material parameters
      matparam%nuparam = 6
      ! Number of table material parameters
      matparam%ntable = 10 
      ! Number of user variables 
      nuvar = 0
      ! Number of temporary variables
      nvartmp = 0
!          
      ! Allocation of material parameters tables
      allocate(matparam%iparam(matparam%niparam))
      allocate(matparam%uparam(matparam%nuparam))
      allocate(matparam%table(matparam%ntable))
!    
      ! Real material parameters
      matparam%uparam(1) = young
      matparam%uparam(2) = nu
      matparam%uparam(3) = a11
      matparam%uparam(4) = a12
      matparam%uparam(5) = hosf
      matparam%uparam(6) = mexp
! 
      ! Table material parameters
      matparam%table( 1)%notable = tab_lch00  
      matparam%table( 2)%notable = tab_lch45
      matparam%table( 3)%notable = tab_lch90
      matparam%table( 4)%notable = tab_lchbi
      matparam%table( 5)%notable = tab_lchsh
      matparam%table( 6)%notable = tab_lcr00
      matparam%table( 7)%notable = tab_lcr45
      matparam%table( 8)%notable = tab_lcr90
      matparam%table( 9)%notable = tab_lcrbi
      matparam%table(10)%notable = tab_lcrsh
!
      ! PARMAT table
      parmat(1) = bulk
      parmat(2) = young
      parmat(3) = nu
      parmat(4) = israte
      parmat(5) = asrate
!
      ! PM table
      pm(1)  = rho0
      pm(89) = rho0
!
      ! MTAG variable activation
      mtag%g_epsd = 1
      mtag%l_epsd = 1
      mtag%g_pla  = 1
      mtag%l_pla  = 1
!
      ! Properties compatibility  
      call init_mat_keyword(matparam,"SHELL_ORTHOTROPIC") 
! 
      ! Properties compatibility  
      call init_mat_keyword(matparam ,"ELASTO-PLASTIC")
      call init_mat_keyword(matparam ,"INCREMENTAL" )
      call init_mat_keyword(matparam ,"LARGE_STRAIN")
      call init_mat_keyword(matparam ,"ORTHOTROPIC") 
!
!--------------------------
!     Parameters printout
!--------------------------
      write(iout, 900) trim(titr),mat_id,129
      write(iout,1000)
      if (is_encrypted) then
        write(iout,'(5x,a,//)') 'confidential data'
      else
        write(iout,1100) rho0
        write(iout,1200) young,nu 
        write(iout,1300) hosf,mexp,tab_lch00,tab_lch45,tab_lch90,     &
          tab_lchbi,tab_lchsh,tab_lcr00,tab_lcr45,tab_lcr90,tab_lcrbi,&
          tab_lcrsh
      endif
  900 format(/                                                        &
        5X,A,/,                                                       &
        5X,'MATERIAL NUMBER. . . . . . . . . . . . . . .=',I10/,      &
        5X,'MATERIAL LAW . . . . . . . . . . . . . . . .=',I10/)
 1000 format(/                                                        &
        5X,'-------------------------------------------',/            &
        5X,'MATERIAL MODEL: EXTENDED 3 PARAMETER BARLAT',/,           &
        5X,'-------------------------------------------',/)
 1100 format(/                                                        &
        5X,'INITIAL DENSITY. . . . . . . . . . . . . . .=',1PG20.13/)  
 1200 format(/                                                        &
        5X,'ELASTIC PARAMETERS:                          ',/,         &
        5X,'-------------------                          ',/,         &
        5X,'YOUNG MODULUS E (E). . . . . . . . . . . . .=',1PG20.13/, &
        5X,"POISSON'S RATIO (NU) . . . . . . . . . . . .=",1PG20.13/)
 1300 format(/                                                        &
        5X,'YIELD CRITERION PARAMETERS:                  ',/,         &
        5X,'---------------------------                  ',/,         &
        5X,'HOSFORD PARAMETER (HOSF) . . . . . . . . . .=',1PG20.13/, &
        5X,'BARLAT EXPONENT (M). . . . . . . . . . . . .=',1PG20.13/, &
        5X,'                                             ',/,         &
        5X,'DIR-0  UNIAXIAL LOADING TABLE ID (LCH00) . .=',I10/,      &
        5X,'DIR-45 UNIAXIAL LOADING TABLE ID (LCH45) . .=',I10/,      &
        5X,'DIR-90 UNIAXIAL LOADING TABLE ID (LCH90) . .=',I10/,      &
        5X,'BIAXIAL LOADING TABLE ID (LCHBI) . . . . . .=',I10/,      &
        5X,'SHEAR LOADING TABLE ID (LCHSH) . . . . . . .=',I10/,      &
        5X,'                                             ',/,         &
        5X,'DIR-0  LANKFORD COEF TABLE ID (LCR00). . . .=',I10/,      &
        5X,'DIR-45 LANKFORD COEF TABLE ID (LCR45). . . .=',I10/,      &
        5X,'DIR-90 LANKFORD COEF TABLE ID (LCR90). . . .=',I10/,      &
        5X,'BIAXIAL LANKFORD COEF TABLE ID (LCRBI) . . .=',I10/,      &
        5X,'SHEAR LANKFORD COEF TABLE ID (LCRSH) . . . .=',I10/,      &
        5X,'                                             ',/)
  end subroutine hm_read_mat129
!
end module hm_read_mat129_mod






