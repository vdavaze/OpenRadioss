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
     module law81_upd_mod
       contains
! \brief Update material law 81 to take into account tabulated stiffness
       subroutine law81_upd(                                                   &
         matparam,nfunc   ,ifunc   ,npc     ,pld     ,pm      ,npropm  )
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Modules
! ----------------------------------------------------------------------------------------------------------------------
         use constant_mod
         use matparam_def_mod
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Implicit none
! ----------------------------------------------------------------------------------------------------------------------
         implicit none
#include "my_real.inc"
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Arguments
! ----------------------------------------------------------------------------------------------------------------------
         type(matparam_struct_), intent(inout)     :: matparam
         integer, intent(in)                       :: nfunc
         integer, dimension(nfunc), intent(in)     :: ifunc
         integer, dimension(*), intent(in)         :: npc
         my_real, dimension(*), intent(in)         :: pld
         my_real, dimension(npropm), intent(inout) :: pm
         integer, intent(in)                       :: npropm
         my_real :: finter
         external finter
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Local variables
! ----------------------------------------------------------------------------------------------------------------------
         my_real :: deri,kini,kscale,gini,gscale
! ----------------------------------------------------------------------------------------------------------------------
!                                                   Body
! ----------------------------------------------------------------------------------------------------------------------
!
         !< Bulk modulus function (if exists)
         if (ifunc(1) > 0) then 
           !< Initial bulk modulus (or scale factor)
           kscale = matparam%uparam(1)
           !< Update accounting for the tabulated bulk modulus function
           kini = kscale*finter(ifunc(1),zero,npc,pld,deri)
           !< Save the new value of the bulk modulus
           matparam%bulk = kini
           !< Update PM table
           pm(32)  = matparam%bulk
           pm(100) = matparam%bulk
           pm(107) = two*pm(32)
         endif
!
         !< Shear modulus function (if exists)
         if (ifunc(2) > 0) then
           !< Initial shear modulus (or scale factor)
           gscale = matparam%uparam(2)
           !< Update accounting for the tabulated bulk modulus function
           gini = gscale*finter(ifunc(2),zero,npc,pld,deri)
           !< Save the new value of the bulk modulus
           matparam%shear = gini
           !< Update PM table
           pm(22) = matparam%shear
         endif
!
         !< Update elastic parameters in the material parameters structures
         if ((ifunc(1) > 0).or.(ifunc(2) > 0)) then
           kini = matparam%bulk
           gini = matparam%shear
           matparam%young = nine*kini*gini/(three*kini+gini)
           matparam%nu = (three*kini-two*gini)/(six*kini+two*gini)
           pm(20) = matparam%young
           pm(21) = matparam%nu
           pm(24) = matparam%young/(one - (matparam%nu)**2)
         endif
!
       end subroutine law81_upd
     end module law81_upd_mod

