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
      module hm_read_mat136_mod
      contains
        subroutine hm_read_mat136(                                             &
          matparam ,nuvar    ,nfunc    ,parmat   ,unitab   ,mat_id   ,titr   , &
          mtag     ,nvartmp  ,lsubmodel,iout     )
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
          use unitab_mod
          use submodel_mod
          use matparam_def_mod
          use elbuftag_mod
          use constant_mod
          use mat_table_copy_mod
          use hm_option_read_mod
          use message_mod
          use precision_mod, only: WP
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
          implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          type(matparam_struct_) ,intent(inout) :: matparam  !< Material parameters structure
          integer,                intent(inout) :: nuvar     !< Number of user variables
          integer,                intent(inout) :: nfunc     !< Number of functions
          real(kind=WP), dimension(100),intent(inout) :: parmat !< Material parameter local array
          type (unit_type_),      intent(in)    :: unitab    !< Units table
          integer,                intent(in)    :: mat_id    !< Material identification number
          character(len=nchartitle),intent(in)  :: titr      !< Material title
          type(mlaw_tag_),        intent(inout) :: mtag      !< Material tags structure
          integer,                intent(inout) :: nvartmp   !< Number of temporary variables
          type(submodel_data), dimension(nsubmod), intent(in) :: lsubmodel !< Submodel data structure
          integer, intent(in)                   :: iout      !< Output file number
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
          real(kind=WP) :: ec, nuc, rho0, bulk, shear, lambda_m, mu_m,         &
            lambda_b, mu_b, tshear, f_c, f_t, kfiss1, kfiss2, dmax,            &
            qp1, qp2, gamma, cn1x, cn1y, cm1x, cm1y, cm1xy, cn2x, cn2y, cm2x,  &
            cm2y, cm2xy
          integer :: i, nlayer 
          real(kind=WP), dimension(:), allocatable :: omega_x, omega_y, rho_x, &
            rho_y, sigy
          integer :: ilaw
          logical :: is_available,is_encrypted
!-----------------------------------------------
!   S o u r c e   L i n e s
!-----------------------------------------------
          is_encrypted = .false.
          is_available = .false.
          ilaw = 136
          !----------------------------------------------------------------------------------
          call hm_option_is_encrypted(is_encrypted)
          !----------------------------------------------------------------------------------
          !< Density
          call hm_get_floatv('MAT_RHO'          ,rho0   ,is_available, lsubmodel, unitab)
          !----------------------------------------------------------------------------------
          !< 1st line of material card
          call hm_get_floatv('MAT_EC'           ,ec     ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_NUC'          ,nuc    ,is_available, lsubmodel, unitab)
          !----------------------------------------------------------------------------------
          !< 2nd line of material card
          call hm_get_floatv('MAT_FT'           ,f_t    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_FC'           ,f_c    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_GAMMA'        ,gamma  ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_QP1'          ,qp1    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_QP2'          ,qp2    ,is_available, lsubmodel, unitab)
          !----------------------------------------------------------------------------------
          !< 3rd line of material card
          call hm_get_intv  ('NLAYER_REINF'     ,nlayer ,is_available, lsubmodel)
          !----------------------------------------------------------------------------------
          !< 4th line of material card
          allocate(sigy(nlayer),omega_x(nlayer),omega_y(nlayer),rho_x(nlayer),rho_y(nlayer))
          do i = 1, nlayer
            call hm_get_float_array_index('MAT_SIGY'   ,sigy(i)   ,i,is_available,lsubmodel,unitab)
            call hm_get_float_array_index('MAT_OMEGA_X',omega_x(i),i,is_available,lsubmodel,unitab)
            call hm_get_float_array_index('MAT_OMEGA_Y',omega_y(i),i,is_available,lsubmodel,unitab)
            call hm_get_float_array_index('MAT_RHO_X'  ,rho_x(i)  ,i,is_available,lsubmodel,unitab)
            call hm_get_float_array_index('MAT_RHO_Y'  ,rho_y(i)  ,i,is_available,lsubmodel,unitab)
          enddo
          !----------------------------------------------------------------------------------
          !< 5th line of material card
          call hm_get_float_array_index('MAT_CN1X'   ,cn1x    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CN1Y'   ,cn1y    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM1X'   ,cm1x    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM1Y'   ,cm1y    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM1XY'  ,cm1xy   ,1,is_available,lsubmodel,unitab)
          !-----------------------------------------------------------------------------------
          !< 6th line of material card
          call hm_get_float_array_index('MAT_CN2X'   ,cn2x    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CN2Y'   ,cn2y    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM2X'   ,cm2x    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM2Y'   ,cm2y    ,1,is_available,lsubmodel,unitab)
          call hm_get_float_array_index('MAT_CM2XY'  ,cm2xy   ,1,is_available,lsubmodel,unitab)
!
          !----------------------------------------------------------------------------------
          !< Parameters default values
          !----------------------------------------------------------------------------------
!
          !----------------------------------------------------------------------------------
          !< Elastic constants
          !----------------------------------------------------------------------------------
          !< Membrane elastic parameters
          lambda_m = ec*nuc/(one - nuc*nuc)
          mu_m     = ec/(two*(one + nuc))
          !< Bending elastic parameters
          lambda_b = ec*nuc/(12.0d0*(one - nuc*nuc))
          mu_b     = ec/(24.0d0*(one + nuc))
          !< Shear modulus
          shear    = ec/(two*(one + nuc))
          !< Bulk modulus
          bulk     = ec/(three*(one - two*nuc))
!
          !----------------------------------------------------------------------------------
          !< Damage constants
          !----------------------------------------------------------------------------------
          dmax = (one - qp2/qp1)/(qp2/qp1 - gamma)
!
          !----------------------------------------------------------------------------------
          !< Filling buffer tables
          !----------------------------------------------------------------------------------
          !< Number of integer material parameters
          matparam%niparam = 1
          !< Number of real material parameters
          matparam%nuparam = 8 + 5*nlayer
          !< Number of user variables
          nuvar = 3
          !< Number of functions
          nfunc = 0
          !< Number of tables and temporary variables
          matparam%ntable = 0
          nvartmp = 0
!
          !< Allocation of material parameters tables
          allocate(matparam%iparam(matparam%niparam))
          allocate(matparam%uparam(matparam%nuparam))
          allocate(matparam%table (matparam%ntable ))
!
          !< Integer material parameters
          matparam%iparam(1)  = nlayer
!
          !< Real material parameters
          matparam%young      = ec
          matparam%nu         = nuc
          matparam%shear      = shear
          matparam%bulk       = bulk
          matparam%uparam(1)  = lambda_m
          matparam%uparam(2)  = mu_m
          matparam%uparam(3)  = lambda_b
          matparam%uparam(4)  = mu_b
          matparam%uparam(5)  = f_t
          matparam%uparam(6)  = f_c
          matparam%uparam(7)  = gamma
          matparam%uparam(8)  = dmax
          do i = 1, nlayer
            matparam%uparam(8 + 5*(i-1) + 1) = sigy(i)
            matparam%uparam(8 + 5*(i-1) + 2) = omega_x(i)
            matparam%uparam(8 + 5*(i-1) + 3) = omega_y(i)
            matparam%uparam(8 + 5*(i-1) + 4) = rho_x(i)
            matparam%uparam(8 + 5*(i-1) + 5) = rho_y(i)
          enddo
!
          !< PARMAT table
          parmat(1)  = bulk
          parmat(2)  = ec
          parmat(3)  = nuc
          parmat(16) = 2 
          parmat(17) = two*shear/(bulk+four_over_3*shear)
!
          !< Initial and reference density
          matparam%rho0 = rho0
          matparam%rho  = rho0
!
          !< MTAG variable activation
          mtag%g_epsd = 1
          mtag%l_epsd = 1
          mtag%g_pla  = 1
          mtag%l_pla  = 1
          mtag%g_dmg  = 3
          mtag%l_dmg  = 3
          mtag%l_sigb = 5
!
          ! Number of output modes 
          matparam%nmod = 2
          allocate(matparam%mode(matparam%nmod))
          matparam%mode(1) = "Positive bending damage"
          matparam%mode(2) = "Negative bending damage"
!
          !< Properties compatibility
          call init_mat_keyword(matparam,"SHELL_ISOTROPIC")
!
          !< Material model keywords
          call init_mat_keyword(matparam ,"INCREMENTAL"   )
          call init_mat_keyword(matparam ,"LARGE_STRAIN"  )
          call init_mat_keyword(matparam ,"HOOK"          )
          call init_mat_keyword(matparam ,"ISOTROPIC"     )
!
          !----------------------------------------------------------------------------------
          !< Listing output
          !----------------------------------------------------------------------------------
          write(iout,1001) trim(titr),mat_id,ilaw
          write(iout,1000)
          if (is_encrypted) then
            write(iout,'(5X,A,//)') 'CONFIDENTIAL DATA'
          else
            write(iout,1002) rho0
            write(iout,1003) ec,nuc
            write(iout,1004) f_t,f_c,gamma,qp1,qp2,dmax
            write(iout,1005) 
            do i = 1,nlayer
              write(iout,1006) i,sigy(i),omega_x(i),omega_y(i),rho_x(i),rho_y(i)
            enddo
            write(iout,1007)
          endif
!
          !< Array deallocation
          deallocate(sigy,omega_x,omega_y,rho_x,rho_y)
!
          !----------------------------------------------------------------------------------
          !< Output formats
          !----------------------------------------------------------------------------------
1000      format(/                                                                 &
            5X,"========================================================",/        &
            5X,"       MATERIAL MODEL: GLOBAL REINFORCED CONCRETE       ",/,       &
            5X,"========================================================",/)
1001      format(/                                                                 &
            5X,A,/,                                                                &
            5X,"MATERIAL NUMBER . . . . . . . . . . . . . . . . . . . .=",I10/,    &
            5X,"MATERIAL LAW. . . . . . . . . . . . . . . . . . . . . .=",I10/)
1002      format(/                                                                 &
            5X,"INITIAL DENSITY . . . . . . . . . . . . . . . . . . . .=",1PG20.13/)
1003      format(/                                                                 &
            5X,"ELASTIC PARAMETERS:                                     ",/,       &
            5X,"-------------------                                     ",/,       &
            5X,"YOUNG MODULUS (E) . . . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"POISSON RATIO (NU). . . . . . . . . . . . . . . . . . .=",1PG20.13/)
1004      format(/                                                                 &
            5X,"BENDING DAMAGE PARAMETERS:                              ",/,       &
            5X,"--------------------------                              ",/,       &
            5X,"TENSILE STRENGTH (FT) . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"COMPRESSIVE STRENGTH (FC) . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"DAMAGE SLOPE RATIO BEFORE/AFTER CRACKING (GAMMA). . . .=",1PG20.13/&
            5X,"SLOPE QUOTIENT FOR POSITIVE BENDING (QP1) . . . . . . .=",1PG20.13/&
            5X,"SLOPE QUOTIENT FOR NEGATIVE BENDING (QP2) . . . . . . .=",1PG20.13/&
            5X,"MAXIMUM DAMAGE COMPUTED . . . . . . . . . . . . . . . .=",1PG20.13/)
1005      format(/                                                                 &
            5X,"STEEL REINFORCEMENT PLASTICITY PARAMETERS:              ",/,       &
            5X,"------------------------------------------              ",/)
1006      format(/                                                                 &
            5X,"REINFORCEMENT LAYER NO # ",I3/                                     &
            5X,"STEEL YIELD STRESS (SIGY) . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"AREA OF THE REINFORCEMENT IN X DIRECTION (OMEGA_X). . .=",1PG20.13/&
            5X,"AREA OF THE REINFORCEMENT IN Y DIRECTION (OMEGA_Y). . .=",1PG20.13/&
            5X,"POSITION IN THICKNESS OF REINF. IN X DIRECTION (RHO_X).=",1PG20.13/&
            5X,"POSITION IN THICKNESS OF REINF. IN Y DIRECTION (RHO_Y).=",1PG20.13/)
1007 format(/                                                                  &
            5X,"========================================================",/)
!
        end subroutine hm_read_mat136
      end module hm_read_mat136_mod
