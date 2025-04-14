!Copyright>        OpenRadioss
!Copyright>        Copyright (C) 1986-2025 Altair Engineering Inc.
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
      module hm_read_mat88_mod
      contains
      subroutine hm_read_mat88(                                                &
                   matparam ,nvartmp  ,parmat   ,unitab   ,mat_id   ,titr     ,&
                   mtag     ,lsubmodel,iout     ,nuvar    ,ilaw     ,ntable   ,&
                   table    ,imatvis  ,israte   ,maxfunc  )
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
        use unitab_mod
        use submodel_mod
        use matparam_def_mod    
        use elbuftag_mod      
        use constant_mod    
        use hm_option_read_mod 
        use table_mod
        use message_mod
        use func_table_copy_mod
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
        implicit none 
#include  "my_real.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        type(matparam_struct_),intent(inout) :: matparam !< Material parameters data structure
        integer, intent(inout) :: nvartmp                !< Number of temporary variables
        my_real, intent(inout) :: parmat(100)            !< Material parameters local table
        type(unit_type_),intent(in) :: unitab             !< Unit conversion table
        integer, intent(in) :: mat_id                    !< Material ID
        character(len=nchartitle), intent(in) :: titr    !< Material title
        type(mlaw_tag_),intent(inout) :: mtag            !< Material tag data structure
        type(submodel_data), dimension(nsubmod), intent(in) :: lsubmodel !< Submodel data structure
        integer, intent(in) :: iout                      !< Output unit
        integer, intent(inout) :: nuvar                  !< Number of user variables
        integer, intent(in) :: ilaw                      !< Material law number
        integer, intent(in) :: ntable                    !< Number of tables
        type(ttable),dimension(ntable),intent(in) :: table !< Table data structure
        integer, intent(inout) :: imatvis                !< Material viscosity flag
        integer, intent(inout) :: israte                 !< Strain rate filtering flag
        integer, intent(in) :: maxfunc                   !< Maximum number of functions
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        my_real ::                                                             &
          k,nu,g,rate(maxfunc+1),hys,rho0,rhor,bulk,fcut,yfac(maxfunc+1),      &
          yfac_unl,shape,gs,e,x1scale,x2scale                            
        integer ::                                                             &
          i,ii,nl,ifunc(maxfunc+1),ifunc_unload,itens,iunl_for,iadd
        logical is_available,is_encrypted
        my_real, parameter :: zep495 = 0.495d0        
!-----------------------------------------------
!     S o u r c e 
!-----------------------------------------------
        is_encrypted = .false.
        is_available = .false.
        imatvis = 1 
        iadd = 0
!-------------------------------------------------------------------------------
        call hm_option_is_encrypted(is_encrypted)
!-------------------------------------------------------------------------------
        !< Density
        call hm_get_floatv('MAT_RHO'      ,rho0  ,is_available, lsubmodel, unitab)
        call hm_get_floatv('Refer_Rho'    ,rhor  ,is_available, lsubmodel, unitab)
!-------------------------------------------------------------------------------
        !< 1st line of material card
        call hm_get_floatv('LAW88_Nu'     ,nu    ,is_available, lsubmodel, unitab)
        call hm_get_floatv('LAW88_K'      ,bulk  ,is_available, lsubmodel, unitab)
        call hm_get_floatv('LAW88_Fcut'   ,fcut  ,is_available, lsubmodel, unitab)
        call hm_get_intv  ('LAW88_Fsmooth',israte,is_available, lsubmodel)      
        call hm_get_intv  ('LAW88_NL'     ,nl    ,is_available, lsubmodel)
        ! -> Check if no loading curve is defined
        if (nl == 0) then
          call ancmsg(msgid = 866,                                             &                                                                                                                                     
                      msgtype = msgerror,                                      &
                      anmode = aninfo_blind,                                   &
                      i1 = mat_id,                                             &
                      c1 = titr)
        endif 
!-------------------------------------------------------------------------------
        !< 2nd line of material card
        call hm_get_intv  ('LAW88_fct_IDunL',ifunc_unload,is_available, lsubmodel)
        call hm_get_floatv('LAW88_FscaleunL',yfac_unl    ,is_available, lsubmodel, unitab)
        if (yfac_unl == zero) then
          call hm_get_floatv_dim('LAW88_FscaleunL',yfac_unl,is_available,lsubmodel,unitab)
        endif
        call hm_get_floatv('LAW88_Hys'      ,hys         ,is_available, lsubmodel, unitab)
        call hm_get_floatv('LAW88_Shape'    ,shape       ,is_available, lsubmodel, unitab)
        call hm_get_intv  ('LAW88_Tension'  ,itens       ,is_available, lsubmodel)
!-------------------------------------------------------------------------------
        !< 3rd line of material card
        ! -> Read the loading curves 
        do i = 1,nl
          call hm_get_int_array_index  ('LAW88_arr1',ifunc(i),i,is_available, lsubmodel)
          call hm_get_float_array_index('LAW88_arr2',yfac(i) ,i,is_available, lsubmodel, unitab)
          if (yfac(i) == zero) then
            call hm_get_float_array_index_dim('LAW88_arr2',yfac(i),i,is_available,lsubmodel, unitab)
          endif
          call hm_get_float_array_index('LAW88_arr3',rate(i) ,i,is_available, lsubmodel, unitab)
          if (rate(i) == zero) then
            call hm_get_float_array_index_dim('LAW88_arr3',rate(i),i,is_available,lsubmodel, unitab)
          endif
        enddo      
!
        !< Check values and set default values
        do i = 2,nl
          if (rate(i) < rate(i-1)) then
            call ancmsg(msgid=478,                                             &
                        msgtype=msgerror,                                      &
                        anmode=aninfo_blind_1,                                 &
                        i1=mat_id,                                             &
                        c1=titr)
            exit
          endif
        enddo
        if (shape == zero) shape = one
        if (hys   == zero) hys   = one
        if (nu    == zero) nu    = zep495
        gs = three_half*bulk*(one - two*nu)/(one + nu)
        e  = two*gs*(one + nu)
        if (gs <= 0) then
          call ancmsg(msgid = 828,                                             &
                      msgtype = msgerror,                                      &
                      anmode = anstop,                                         &
                      i1 = mat_id,                                             &
                      c1 = titr)
        endif
        if (fcut == zero .and. nl > 1) then
          fcut = ep03*unitab%fac_t_work
          israte = 1  
        endif
!-------------------------------------------------------------------------------
        !< Filling buffer tables
!------------------------------------------------------------------------------- 
        !< Number of integer material parameters
        matparam%niparam = 2
        !< Number of real material parameters
        matparam%nuparam = 2
        !< Number of user variables 
        nuvar = 32 
        !< Number of tables and temporary variables
        matparam%ntable = 2
        nvartmp = 2
!          
        !< Allocation of material parameters tables
        allocate(matparam%iparam(matparam%niparam))
        allocate(matparam%uparam(matparam%nuparam))
        allocate(matparam%table (matparam%ntable ))
!     
        !< Integer material parameter
        matparam%iparam(1)  = itens
        matparam%iparam(2)  = iunl_for
!    
        !< Real material parameters
        matparam%young      = e
        matparam%nu         = nu
        matparam%shear      = gs
        matparam%bulk       = bulk
        matparam%uparam(1)  = hys
        matparam%uparam(2)  = shape
!    
        !< Transform series of functions into material table
        ! -> Loading curves
        x1scale = one 
        x2scale = one
        call func_table_copy(matparam%table(1) ,titr     ,mat_id   ,           &
                 nl       ,ifunc(1:nl),rate(1:nl),x1scale,x2scale  ,           &
                 yfac(1:nl),ntable    ,table     ,ierr   )      
        ! -> Unloading curves
        iunl_for = 0
        if (ifunc_unload > 0) then 
          iunl_for = 1
          yfac(nl+1) = yfac_unl
          rate(nl+1) = zero
          ifunc(nl+1) = ifunc_unload
          call func_table_copy(matparam%table(2) ,titr     ,mat_id   ,         &
                   1        ,ifunc(nl+1),rate(nl+1),x1scale,x2scale  ,         &
                   yfac(nl+1),ntable    ,table    ,ierr     )   
        elseif (hys /= zero) then
          if (nl == 1) then 
            iunl_for = 2 ! based on the energy 
          else
            iunl_for = 3 ! based on the energy
          endif
          hys = abs(hys)
        else
          if (nl == 1) then
            iunl_for = 0  ! no unloading curve, 
          else
            iunl_for = 1
            yfac(nl+1) = yfac_unl
            rate(nl+1) = zero
            ifunc(nl+1) = ifunc_unload
            call func_table_copy(matparam%table(2) ,titr     ,mat_id   ,       &
                     1        ,ifunc(nl+1),rate(nl+1),x1scale,x2scale  ,       &
                     yfac(nl+1),ntable    ,table    ,ierr    )    
          endif            
        endif    
! 
        !< Material density
        if (rhor == zero) rhor = rho0
        matparam%rho0 = rho0
        matparam%rho  = rhor
!      
        !< PARMAT table
        parmat(1)  = two*gs 
        parmat(2)  = e
        parmat(3)  = nu
        parmat(4)  = israte
        parmat(5)  = fcut
        parmat(16) = 2
        parmat(17) = two*gs/(bulk + four_over_3*gs)
!
        ! MTAG variable activation
        mtag%l_epsd = 1
        mtag%g_epsd = 1
! 
        !< Material model keywords
        call init_mat_keyword(matparam,"COMPRESSIBLE")
        call init_mat_keyword(matparam,"TOTAL")
        call init_mat_keyword(matparam,"HOOK")
!
        !< Properties compatibility
        call init_mat_keyword(matparam,"SOLID_ISOTROPIC")
        call init_mat_keyword(matparam,"SHELL_ISOTROPIC")
!
!-------------------------------------------------------------------------------
        !< Printing out the material data
!-------------------------------------------------------------------------------
        write(iout,1010) trim(titr),mat_id,ilaw
        write(iout,1000)     
        if (is_encrypted) then
          write(iout,'(5X,A,//)')'CONFIDENTIAL DATA'
        else     
          write(iout,1020) rho0 
          write(iout,1100) nu,bulk,itens,nl-iadd
          write(iout,1200)(ifunc(i),yfac(i),rate(i),i=1,nl)     
          if(iunl_for == 1) then
           write(iout,1300) ifunc(nl+1),yfac_unl 
          elseif(iunl_for == 2 .or. iunl_for == 3) then
           write(iout,1400) hys, shape
          endif
          write(iout,1500) itens
        endif     
!-----------------
 1000 format(/                                                                 &
       5X,'-------------------------------------------',/,                     &
       5X,'     MATERIAL MODEL :TABULATED OGDEN       ',/,                     &
       5X,'-------------------------------------------',/)                     
 1010 format(/                                                                 &
       5X,A,/,                                                                 &
       5X,'MATERIAL NUMBER. . . . . . . . . . . . . .=',I10/,                  &
       5X,'MATERIAL LAW . . . . . . . . . . . . . . .=',I10/) 
 1020 format(/                                                                 &
       5X,'INITIAL DENSITY. . . . . . . . . . . . . .=',1PG20.13/) 
 1100 format(/                                                                 &
       5X,'POISSON RATIO. . . . . . . . . .  . . . . =',1PG20.13/,             &
       5X,'BULK MODULUS. . . . . . . . . . . . . . . =',1PG20.13/,             &
       5X,'STRAIN RATE EFFECT FLAG  . .. . . . . . . =',I10/,                  &
       5X,'NUMBER OF LOADING  FUNCTION . . .. . . . .=',I10/)
 1200 format(/                                                                 &
       5X,'LOADING STRESS-STRAIN FUNCTION NUMBER. . .=',I10/,                  &
       5X,'STRESS SCALE FACTOR. . . . . . . . . . . .=',1PG20.13/,             &
       5X,'STRAIN RATE . . . . . . . . . . . . . . . =',1PG20.13/)   
 1300 format(/                                                                 &
       5X,'UNLOADING STRESS-STRAIN FUNCTION NUMBER. .=',I10/,                  &
       5X,'STRESS SCALE FACTOR. . . . . . . . . . . .=',1PG20.13/)    
 1400 format(/                                                                 &
       5X,'HYSTERETIC UNLOADING FACTOR. . . . .  . . =',1PG20.13/,             &
       5X,'SHAPE UNLOADING FACTOR. . . . . . . . . . =',1PG20.13/) 
 1500 format(/                                                                 &
       5X,'ITENSION : PARAMETER FOR UNLOADING . . . .=',I10/)
!-----------------
      end subroutine hm_read_mat88
      end module hm_read_mat88_mod 
