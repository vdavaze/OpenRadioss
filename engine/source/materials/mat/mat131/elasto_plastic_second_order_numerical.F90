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
      module elasto_plastic_second_order_numerical_mod
      contains
      subroutine elasto_plastic_second_order_numerical(                        &
        matparam ,nel      ,eltype   ,icrit    ,                               &
        signxx   , signyy  ,signzz   ,signxy   ,signyz   ,signzx   ,           &
        N2       )
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use precision_mod, only : WP
        use constant_mod
        use yield_criterion_barlat1989_mod
        use yield_criterion_barlat2000_mod
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        type(matparam_struct_),        intent(in)    :: matparam
        integer,                       intent(in)    :: nel, eltype, icrit
        real(kind=WP), dimension(nel), intent(in)    :: signxx, signyy, signzz
        real(kind=WP), dimension(nel), intent(in)    :: signxy, signyz, signzx
        real(kind=WP), dimension(nel,6,6), intent(out) :: N2
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        real(kind=WP), dimension(nel) :: sig_p, sig_m, norm
        real(kind=WP), dimension(nel) :: seq_p,  seq_m
        real(kind=WP), dimension(nel) :: sxx_p, syy_p, szz_p, sxy_p, syz_p, szx_p
        real(kind=WP), dimension(nel) :: sxx_m, syy_m, szz_m, sxy_m, syz_m, szx_m
        real(kind=WP), dimension(nel) :: nxx_p, nyy_p, nzz_p, nxy_p, nyz_p, nzx_p
        real(kind=WP), dimension(nel) :: nxx_m, nyy_m, nzz_m, nxy_m, nyz_m, nzx_m
        real(kind=WP) :: h
        integer :: j, jj, njdir
        integer, dimension(6) :: jlist
!===============================================================================
!
        !< Initialisation of the output array
        N2(1:nel,1:6,1:6) = zero
!
        !< Computation of the pertubation step size
        h = 1.0d-8
        ! -> For solid elements 
        if (eltype == 1) then
          norm(1:nel) = signxx(1:nel)**2 + signyy(1:nel)**2                    &
                      + signzz(1:nel)**2 + 2*signxy(1:nel)**2                  &
                      + 2*signyz(1:nel)**2 + 2*signzx(1:nel)**2
        ! -> For shell elements
        else
          norm(1:nel) = signxx(1:nel)**2 + signyy(1:nel)**2                    &
                      + 2*signxy(1:nel)**2
        endif
        norm(1:nel) = sqrt(norm(1:nel))
        norm(1:nel) = max(norm(1:nel),one)
!
        !< List of the perturbed directions : 
        ! -> 6 for solids
        if (eltype == 1) then
          njdir = 6
          jlist = [1, 2, 3, 4, 5, 6]
        ! -> 3 for shells
        else
          njdir = 3
          jlist = [1, 2, 4, 0, 0, 0] 
        endif
!
        !< Loop over the perturbed directions
        do jj = 1, njdir
          j = jlist(jj)
!
          !< Perturbed stress states : +h and -h in the j-th direction
          sxx_p(1:nel) = signxx(1:nel)
          syy_p(1:nel) = signyy(1:nel)
          szz_p(1:nel) = signzz(1:nel)
          sxy_p(1:nel) = signxy(1:nel)
          syz_p(1:nel) = signyz(1:nel)
          szx_p(1:nel) = signzx(1:nel)
!
          sxx_m(1:nel) = signxx(1:nel)
          syy_m(1:nel) = signyy(1:nel)
          szz_m(1:nel) = signzz(1:nel)
          sxy_m(1:nel) = signxy(1:nel)
          syz_m(1:nel) = signyz(1:nel)
          szx_m(1:nel) = signzx(1:nel)
!
          select case (j)
            case (1)
              sxx_p(1:nel) = signxx(1:nel) + h*norm(1:nel)
              sxx_m(1:nel) = signxx(1:nel) - h*norm(1:nel)
            case (2)
              syy_p(1:nel) = signyy(1:nel) + h*norm(1:nel)
              syy_m(1:nel) = signyy(1:nel) - h*norm(1:nel)
            case (3)
              szz_p(1:nel) = signzz(1:nel) + h*norm(1:nel)
              szz_m(1:nel) = signzz(1:nel) - h*norm(1:nel)
            case (4)
              sxy_p(1:nel) = signxy(1:nel) + h*norm(1:nel)
              sxy_m(1:nel) = signxy(1:nel) - h*norm(1:nel)
            case (5)
              syz_p(1:nel) = signyz(1:nel) + h*norm(1:nel)
              syz_m(1:nel) = signyz(1:nel) - h*norm(1:nel)
            case (6)
              szx_p(1:nel) = signzx(1:nel) + h*norm(1:nel)
              szx_m(1:nel) = signzx(1:nel) - h*norm(1:nel)
          end select
!
          !< Evaluation of the yield criterion and its gradient for the 
          ! perturbed stress states
          select case(icrit)
            !-------------------------------------------------------------------
            !< Barlat 89 yield criterion
            !-------------------------------------------------------------------
            case(4)
              call yield_criterion_barlat1989(                                 &          
                matparam ,nel      ,seq_p    ,sxx_p    ,syy_p    ,sxy_p    ,   &
                nxx_p    ,nyy_p    ,nzz_p    ,nxy_p    ,nyz_p    ,nzx_p    ) 
              call yield_criterion_barlat1989(                                 &          
                matparam ,nel      ,seq_m    ,sxx_m    ,syy_m    ,sxy_m    ,   &
                nxx_m    ,nyy_m    ,nzz_m    ,nxy_m    ,nyz_m    ,nzx_m    ) 
            !-------------------------------------------------------------------
            !< Barlat 2000 yield criterion
            !-------------------------------------------------------------------
            case(5)
              call yield_criterion_barlat2000(                                 &          
                matparam ,nel      ,seq_p    ,sxx_p    ,syy_p    ,sxy_p    ,   &
                nxx_p    ,nyy_p    ,nzz_p    ,nxy_p    ,nyz_p    ,nzx_p    )
              call yield_criterion_barlat2000(                                 &          
                matparam ,nel      ,seq_m    ,sxx_m    ,syy_m    ,sxy_m    ,   &
                nxx_m    ,nyy_m    ,nzz_m    ,nxy_m    ,nyz_m    ,nzx_m    ) 
          end select
!
          !< Finite differences : N2(:,i,j) = (norm_i(+h) - norm_i(-h)) / 2h
          N2(1:nel,1,j) = (nxx_p(1:nel) - nxx_m(1:nel)) / (two*h*norm(1:nel))
          N2(1:nel,2,j) = (nyy_p(1:nel) - nyy_m(1:nel)) / (two*h*norm(1:nel))
          N2(1:nel,4,j) = (nxy_p(1:nel) - nxy_m(1:nel)) / (two*h*norm(1:nel))
          if (eltype == 1) then 
            N2(1:nel,3,j) = (nzz_p(1:nel) - nzz_m(1:nel)) / (two*h*norm(1:nel))
            N2(1:nel,5,j) = (nyz_p(1:nel) - nyz_m(1:nel)) / (two*h*norm(1:nel))
            N2(1:nel,6,j) = (nzx_p(1:nel) - nzx_m(1:nel)) / (two*h*norm(1:nel))
          endif
!
        enddo
!
      end subroutine elasto_plastic_second_order_numerical
      end module elasto_plastic_second_order_numerical_mod
