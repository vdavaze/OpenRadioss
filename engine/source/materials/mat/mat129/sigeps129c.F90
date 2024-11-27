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
      module sigeps129c_mod
        contains
  ! ======================================================================================================================
  ! \brief   
  ! \details 
  ! ======================================================================================================================
         subroutine sigeps129c(                                                &
                      nel      ,matparam  ,rho      ,thk      ,thkly    ,      &
                      depsxx   ,depsyy    ,depsxy   ,depsyz   ,depszx   ,      &   
                      sigoxx   ,sigoyy    ,sigoxy   ,sigozx   ,sigoyz   ,      &
                      signxx   ,signyy    ,signxy   ,signzx   ,signyz   ,      &
                      off      ,sigy      ,etse     ,ssp      ,gs       ,      &
                      pla      ,dpla      ,epsp     ,vartmp   ,nvartmp  ,      &
                      table    ,ntable    )
!---------------------------------------------- -
!   M o d u l e s
!-----------------------------------------------
          use matparam_def_mod 
          use constant_mod      
          use table_mod
!-----------------------------------------------
!   I m p l i c i t   T y p e s
!-----------------------------------------------
          implicit none 
#include  "my_real.inc"
#include  "units_c.inc"
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
          integer, intent(in)                    :: nel      !< number of elements in the group
          type(matparam_struct_), intent(in)     :: matparam !< material parameters data
          my_real, dimension(nel), intent(in)    :: rho      !< material density
          my_real, dimension(nel), intent(inout) :: thk      !< element thikness
          my_real, dimension(nel), intent(in)    :: thkly    !< integration point related thikness
          my_real, dimension(nel), intent(in)    :: depsxx   !< strain increment xx
          my_real, dimension(nel), intent(in)    :: depsyy   !< strain increment yy
          my_real, dimension(nel), intent(in)    :: depsxy   !< strain increment xy
          my_real, dimension(nel), intent(in)    :: depsyz   !< strain increment yz
          my_real, dimension(nel), intent(in)    :: depszx   !< strain increment zx
          my_real, dimension(nel), intent(in)    :: sigoxx   !< old stress xx
          my_real, dimension(nel), intent(in)    :: sigoyy   !< old stress yy
          my_real, dimension(nel), intent(in)    :: sigoxy   !< old stress xy
          my_real, dimension(nel), intent(in)    :: sigozx   !< old stress zx
          my_real, dimension(nel), intent(in)    :: sigoyz   !< old stress yz
          my_real, dimension(nel), intent(inout) :: signxx   !< new stress xx
          my_real, dimension(nel), intent(inout) :: signyy   !< new stress yy
          my_real, dimension(nel), intent(inout) :: signxy   !< new stress xy
          my_real, dimension(nel), intent(inout) :: signzx   !< new stress zx
          my_real, dimension(nel), intent(inout) :: signyz   !< new stress yz
          my_real, dimension(nel), intent(inout) :: off      !< element deletion flag
          my_real, dimension(nel), intent(inout) :: sigy     !< yield stress
          my_real, dimension(nel), intent(inout) :: etse     !< ratio of rigidity
          my_real, dimension(nel), intent(inout) :: ssp      !< sound speed
          my_real, dimension(nel), intent(in)    :: gs       !< shear modulus with correction factor
          my_real, dimension(nel), intent(inout) :: pla      !< plastic strain
          my_real, dimension(nel), intent(inout) :: dpla     !< plastic strain increment
          my_real, dimension(nel), intent(inout) :: epsp     !< equiv. strain rate
          integer, dimension(nel,nvartmp), intent(inout) :: vartmp !< temporary variables
          integer, intent(in)                    :: nvartmp  !< number of temporary variables
          type(ttable), dimension(ntable), intent(in) :: table !< table data structure
          integer, intent(in)                    :: ntable   !< number of tables
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
          integer :: i
          integer, dimension(nel,2) :: ipos
          integer :: tab_lch00,tab_lch45,tab_lch90,tab_lchbi,tab_lchsh
          integer :: tab_lcr00,tab_lcr45,tab_lcr90,tab_lcrbi,tab_lcrsh
          my_real :: young,nu,a11,a12,g,hosf,mexp
          my_real, dimension(nel) :: dezz,hardp,hardr,sigy_00,sigy_45,sigy_90
          my_real, dimension(nel,2) :: xvec
!=============================================================================
! 
          !========================================================================
          ! - Initialisation of computation on time step
          !========================================================================
          !< Recovering real model parameters
          young = matparam%uparam(1)
          nu    = matparam%uparam(2)
          a11   = matparam%uparam(3)
          a12   = matparam%uparam(4)
          g     = matparam%uparam(5)
          hosf  = matparam%uparam(6)
          mexp  = matparam%uparam(7)
          !< Recovering table ids
          tab_lch00 = matparam%table( 1)%notable
          tab_lch45 = matparam%table( 2)%notable
          tab_lch90 = matparam%table( 3)%notable
          tab_lchbi = matparam%table( 4)%notable
          tab_lchsh = matparam%table( 5)%notable
          tab_lcr00 = matparam%table( 6)%notable
          tab_lcr45 = matparam%table( 7)%notable
          tab_lcr90 = matparam%table( 8)%notable
          tab_lcrbi = matparam%table( 9)%notable
          tab_lcrsh = matparam%table(10)%notable
          write(*,*) "param = ", tab_lch00,tab_lch45,tab_lch90
!
          !< Computation of yielding parameters
          xvec(1:nel,1) = pla(1:nel)
          xvec(1:nel,2) = epsp(1:nel)
          !< Yield stress in direction 00 
          ipos(1:nel,1) = vartmp(1:nel,1)
          ipos(1:nel,2) = 1
          call table2d_vinterp_log(table(tab_lch00),0,nel,nel,ipos,xvec,sigy_00,hardp,hardr)
          vartmp(1:nel,1) = ipos(1:nel,1)
          !< Yield stress in direction 45
          ipos(1:nel,1) = vartmp(1:nel,2)
          ipos(1:nel,2) = 1
          call table2d_vinterp_log(table(tab_lch45),0,nel,nel,ipos,xvec,sigy_45,hardp,hardr)
          vartmp(1:nel,2) = ipos(1:nel,1)
          !< Yield stress in direction 90
          ipos(1:nel,1) = vartmp(1:nel,3)
          ipos(1:nel,2) = 1
          call table2d_vinterp_log(table(tab_lch90),0,nel,nel,ipos,xvec,sigy_90,hardp,hardr)
          vartmp(1:nel,3) = ipos(1:nel,1) 
!
          !========================================================================
          ! - Computation of elastic deviatoric stresses and equivalent stress
          !========================================================================
          do i = 1,nel  
            signxx(i) = sigoxx(i) + a11*depsxx(i) + a12*depsyy(i)
            signyy(i) = sigoyy(i) + a11*depsyy(i) + a12*depsxx(i)
            signxy(i) = sigoxy(i) + depsxy(i)*g
            signyz(i) = sigoyz(i) + depsyz(i)*gs(i)
            signzx(i) = sigozx(i) + depszx(i)*gs(i)     
          enddo     
!
          !========================================================================
          ! - Update stress tensor and sound speed
          !========================================================================
          do i=1,nel
            ! Sound speed
            ssp(i)  = sqrt(a11/rho(i))
            ! Coefficient for hourglass 
            etse(i) = one
            ! Update thickness
            thk(i)  = thk(i) + dezz(i)*thkly(i)*off(i)  
            if (1 == 2) write(*,*) epsp(i), sigy(i), pla(i), dpla(i)
          enddo   
!-------------------------------------------------------------------------------------------
         end subroutine sigeps129c
      end module sigeps129c_mod 
