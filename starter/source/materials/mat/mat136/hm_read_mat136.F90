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
          mtag     ,nvartmp  ,lsubmodel,ntable   ,table    ,iout     )
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
          use table_mod
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
          integer, intent(in)                   :: ntable    !< Number of tables
          type(ttable),dimension(ntable),intent(in) :: table !< Tables data structure
          integer, intent(in)                   :: iout      !< Output file number
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
          real(kind=WP) :: e, nu, rho0, bulk, shear, lambda_m, mu_m, lambda_b, &
             mu_b, tshear, f_c, f_t, kfiss1, kfiss2, dmax, qp1, qp2, gamma
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
          call hm_get_floatv('MAT_E'            ,e      ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_NU'           ,nu     ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_FT'           ,f_t    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_FC'           ,f_c    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_GAMMA'        ,gamma  ,is_available, lsubmodel, unitab)
          !----------------------------------------------------------------------------------
          !< 2nd line of material card
          call hm_get_floatv('MAT_QP1'          ,qp1    ,is_available, lsubmodel, unitab)
          call hm_get_floatv('MAT_QP2'          ,qp2    ,is_available, lsubmodel, unitab)
          !----------------------------------------------------------------------------------
!
          !----------------------------------------------------------------------------------
          !< Parameters default values
          !----------------------------------------------------------------------------------
!
          !----------------------------------------------------------------------------------
          !< Elastic constants
          !----------------------------------------------------------------------------------
          !< Membrane elastic parameters
          lambda_m = e*nu/(one - nu*nu)
          mu_m     = e/(two*(one + nu))
          !< Bending elastic parameters
          lambda_b = e*nu/(12.0d0*(one - nu*nu))
          mu_b     = e/(24.0d0*(one + nu))
          !< Shear modulus
          shear    = e/(two*(one + nu))
          !< Bulk modulus
          bulk     = e/(three*(one - two*nu))
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
          matparam%nuparam = 8
          !< Number of user variables
          nuvar = 5
          !< Number of functions
          nfunc = 0
          !< Number of tables and temporary variables
          matparam%ntable = 0
          nvartmp = 0
          if (1 == 2) write(*,*) table(1)%ndim
!
          !< Allocation of material parameters tables
          allocate(matparam%iparam(matparam%niparam))
          allocate(matparam%uparam(matparam%nuparam))
          allocate(matparam%table (matparam%ntable ))
!
          !< Integer material parameters
          matparam%iparam(1)  = 1
!
          !< Real material parameters
          matparam%young     = e
          matparam%nu        = nu
          matparam%shear     = shear
          matparam%bulk      = bulk
          matparam%uparam(1) = lambda_m
          matparam%uparam(2) = mu_m
          matparam%uparam(3) = lambda_b
          matparam%uparam(4) = mu_b
          matparam%uparam(5) = f_t
          matparam%uparam(6) = f_c
          matparam%uparam(7) = gamma
          matparam%uparam(8) = dmax
!
          !< PARMAT table
          parmat(1)  = bulk
          parmat(2)  = e
          parmat(3)  = nu
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
            write(iout,1003) e,nu,f_t,f_c,gamma,qp1,qp2,dmax
          endif
!
          !----------------------------------------------------------------------------------
          !< Output formats
          !----------------------------------------------------------------------------------
1000      format(/                                                                 &
            5X,"-------------------------------------------------------",/         &
            5X,"       MATERIAL MODEL: GLOBAL REINFORCED CONCRETE      ",/,        &
            5X,"-------------------------------------------------------",/)
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
            5X,"POISSON RATIO (NU). . . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"TENSILE STRENGTH (FT) . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"COMPRESSIVE STRENGTH (FC) . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"GAMMA . . . . . . . . . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"SLOPE 1 . . . . . . . . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"SLOPE 2 . . . . . . . . . . . . . . . . . . . . . . . .=",1PG20.13/&
            5X,"MAXIMUM DAMAGE (DMAX) . . . . . . . . . . . . . . . . .=",1PG20.13/)
!
        end subroutine hm_read_mat136
      end module hm_read_mat136_mod
