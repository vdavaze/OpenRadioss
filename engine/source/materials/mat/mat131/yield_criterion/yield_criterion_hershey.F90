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
      module yield_criterion_hershey_mod
      contains
      subroutine yield_criterion_hershey(                                      &
        matparam ,nel      ,nindx    ,indx     ,seq      ,iresp    ,           &
        signxx   ,signyy   ,signzz   ,signxy   ,signyz   ,signzx   ,           &
        normxx   ,normyy   ,normzz   ,normxy   ,normyz   ,normzx   ,           &
        eltype   )  
!----------------------------------------------------------------
!   M o d u l e s
!----------------------------------------------------------------
        use matparam_def_mod
        use constant_mod
        use mvsiz_mod
        use precision_mod, only : WP
        use yield_criterion_vonmises_mod
!----------------------------------------------------------------
!   I m p l i c i t   T y p e s
!----------------------------------------------------------------
        implicit none
!----------------------------------------------------------------
!  I n p u t   A r g u m e n t s
!----------------------------------------------------------------
        type(matparam_struct_),        intent(in)    :: matparam !< Material parameters data
        integer,                       intent(in)    :: nel      !< Number of elements in the group
        integer,                       intent(in)    :: nindx    !< Number of elements to consider in the computation (for partial updates)
        integer,       dimension(nel), intent(in)    :: indx     !< Indices of the elements to consider in the computation (for partial updates)
        real(kind=WP), dimension(nel), intent(inout) :: seq      !< Equivalent stress
        integer,                       intent(in)    :: iresp    !< Precision flag
        real(kind=WP), dimension(nel), intent(in)    :: signxx   !< Current stress xx
        real(kind=WP), dimension(nel), intent(in)    :: signyy   !< Current stress yy
        real(kind=WP), dimension(nel), intent(in)    :: signzz   !< Current stress zz
        real(kind=WP), dimension(nel), intent(in)    :: signxy   !< Current stress xy
        real(kind=WP), dimension(nel), intent(in)    :: signyz   !< Current stress yz
        real(kind=WP), dimension(nel), intent(in)    :: signzx   !< Current stress zx
        real(kind=WP), dimension(nel), intent(inout) :: normxx   !< 1st derivative of equivalent stress wrt stress xx
        real(kind=WP), dimension(nel), intent(inout) :: normyy   !< 1st derivative of equivalent stress wrt stress yy
        real(kind=WP), dimension(nel), intent(inout) :: normzz   !< 1st derivative of equivalent stress wrt stress zz
        real(kind=WP), dimension(nel), intent(inout) :: normxy   !< 1st derivative of equivalent stress wrt stress xy
        real(kind=WP), dimension(nel), intent(inout) :: normyz   !< 1st derivative of equivalent stress wrt stress yz
        real(kind=WP), dimension(nel), intent(inout) :: normzx   !< 1st derivative of equivalent stress wrt stress zx
        integer,                       intent(in)    :: eltype   !< Element type
!----------------------------------------------------------------
!  L o c a l  V a r i a b l e s
!----------------------------------------------------------------
        integer :: offset,i,j,k,ii
        real(kind=WP), dimension(mvsiz,6) :: strs
        real(kind=WP), dimension(mvsiz,3) :: pr_strs
        real(kind=WP), dimension(mvsiz,3,3) :: dir
        real(kind=WP) :: pres,sxx,syy,szz,sxy,syz,szx,nexp
        real(kind=WP) :: norm(nel),dsigeq_dsig1,dsigeq_dsig2,dsigeq_dsig3
        real(kind=WP) :: dsig1_dsig(6),dsig2_dsig(6),dsig3_dsig(6)
        real(kind=WP) :: dsig1_dsigxx,dsig1_dsigyy,dsig1_dsigxy
        real(kind=WP) :: dsig2_dsigxx,dsig2_dsigyy,dsig2_dsigxy
        real(kind=WP) :: center, rootv
        real(kind=WP) :: d2sigeq_dsig12,d2sigeq_dsig1sig2,d2sigeq_dsig1sig3,   &
                         d2sigeq_dsig22,d2sigeq_dsig2sig3,d2sigeq_dsig32,      &
                         d2sigeq_dsig2sig1,d2sigeq_dsig3sig1,d2sigeq_dsig3sig2
!===============================================================================
!
        !=======================================================================
        !< - Hershey yield criterion and its derivatives
        !=======================================================================
        offset = matparam%iparam(2)
        !< Hershey exponent
        nexp = matparam%uparam(offset + 1)
        !< Solid element
        if (eltype == 1) then 
#include "vectorize.inc"
          do ii = 1,nindx
            i = indx(ii)
            !< Normalization of the stress tensor
            pres = -(signxx(i) + signyy(i) + signzz(i))/three
            sxx = signxx(i) + pres
            syy = signyy(i) + pres
            szz = signzz(i) + pres
            sxy = signxy(i)
            syz = signyz(i)
            szx = signzx(i)
            norm(i) = half * (sxx**2 + syy**2 + syy**2) +                      &
                              sxy**2 + syz**2 + szx**2
            norm(i) = max(sqrt(three*norm(i)),one)
            strs(ii,1) = signxx(i)/norm(i)
            strs(ii,2) = signyy(i)/norm(i)
            strs(ii,3) = signzz(i)/norm(i)
            strs(ii,4) = signxy(i)/norm(i)
            strs(ii,5) = signyz(i)/norm(i)
            strs(ii,6) = signzx(i)/norm(i)
          enddo
          !< Compute principal strains and directions
          if (iresp == 1) then
            call valpvecdp_v(strs ,pr_strs ,dir ,nindx)
          else
            call valpvec_v(strs ,pr_strs ,dir ,nindx)
          endif
          do ii = 1,nindx
            i = indx(ii)
            !< Equivalent stress
            seq(i) = half * ((abs(pr_strs(ii,1) - pr_strs(ii,2)))**nexp +      &
                             (abs(pr_strs(ii,2) - pr_strs(ii,3)))**nexp +      &
                             (abs(pr_strs(ii,3) - pr_strs(ii,1)))**nexp )
            if (seq(i) > zero) then
              seq(i) = exp((one/nexp)*log(seq(i)))
            else
              seq(i) = zero
            endif
            !< Derivatives of eq. stress
            normxx(i) = zero
            normyy(i) = zero
            normzz(i) = zero
            normxy(i) = zero
            normyz(i) = zero
            normzx(i) = zero
            if (seq(i) > zero) then 
              dsigeq_dsig1 = half*(seq(i)**(1-nexp))*(                         &
                     ((abs(pr_strs(ii,1) - pr_strs(ii,2)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,1) - pr_strs(ii,2)) -                    &
                     ((abs(pr_strs(ii,3) - pr_strs(ii,1)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,3) - pr_strs(ii,1)))
              dsigeq_dsig2 = half*(seq(i)**(1-nexp))*(                         &
                     ((abs(pr_strs(ii,2) - pr_strs(ii,3)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,2) - pr_strs(ii,3)) -                    &
                     ((abs(pr_strs(ii,1) - pr_strs(ii,2)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,1) - pr_strs(ii,2)))
              dsigeq_dsig3 = half*(seq(i)**(1-nexp))*(                         &
                     ((abs(pr_strs(ii,3) - pr_strs(ii,1)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,3) - pr_strs(ii,1)) -                    &
                     ((abs(pr_strs(ii,2) - pr_strs(ii,3)))**(nexp-1))*         &
                  sign(one,pr_strs(ii,2) - pr_strs(ii,3)))
              !< Derivatives of principal stresses w.r.t. stress tensor
              dsig1_dsig(1) =     dir(ii,1,1)*dir(ii,1,1)
              dsig2_dsig(1) =     dir(ii,1,2)*dir(ii,1,2)
              dsig3_dsig(1) =     dir(ii,1,3)*dir(ii,1,3)
              dsig1_dsig(2) =     dir(ii,2,1)*dir(ii,2,1)
              dsig2_dsig(2) =     dir(ii,2,2)*dir(ii,2,2)
              dsig3_dsig(2) =     dir(ii,2,3)*dir(ii,2,3)
              dsig1_dsig(3) =     dir(ii,3,1)*dir(ii,3,1)
              dsig2_dsig(3) =     dir(ii,3,2)*dir(ii,3,2)
              dsig3_dsig(3) =     dir(ii,3,3)*dir(ii,3,3)
              dsig1_dsig(4) = two*dir(ii,1,1)*dir(ii,2,1)
              dsig2_dsig(4) = two*dir(ii,1,2)*dir(ii,2,2)
              dsig3_dsig(4) = two*dir(ii,1,3)*dir(ii,2,3)
              dsig1_dsig(5) = two*dir(ii,2,1)*dir(ii,3,1)
              dsig2_dsig(5) = two*dir(ii,2,2)*dir(ii,3,2)
              dsig3_dsig(5) = two*dir(ii,2,3)*dir(ii,3,3)
              dsig1_dsig(6) = two*dir(ii,1,1)*dir(ii,3,1)
              dsig2_dsig(6) = two*dir(ii,1,2)*dir(ii,3,2)
              dsig3_dsig(6) = two*dir(ii,1,3)*dir(ii,3,3)
              !< Assembly of the derivative of the eq. stress w.r.t. stress tensor
              normxx(i) = dsigeq_dsig1*dsig1_dsig(1) +                         &
                          dsigeq_dsig2*dsig2_dsig(1) +                         &
                          dsigeq_dsig3*dsig3_dsig(1)
              normyy(i) = dsigeq_dsig1*dsig1_dsig(2) +                         &
                          dsigeq_dsig2*dsig2_dsig(2) +                         &
                          dsigeq_dsig3*dsig3_dsig(2)     
              normzz(i) = dsigeq_dsig1*dsig1_dsig(3) +                         &
                          dsigeq_dsig2*dsig2_dsig(3) +                         &
                          dsigeq_dsig3*dsig3_dsig(3)      
              normxy(i) = dsigeq_dsig1*dsig1_dsig(4) +                         &
                          dsigeq_dsig2*dsig2_dsig(4) +                         &
                          dsigeq_dsig3*dsig3_dsig(4)          
              normyz(i) = dsigeq_dsig1*dsig1_dsig(5) +                         &
                          dsigeq_dsig2*dsig2_dsig(5) +                         &
                          dsigeq_dsig3*dsig3_dsig(5)          
              normzx(i) = dsigeq_dsig1*dsig1_dsig(6) +                         &
                          dsigeq_dsig2*dsig2_dsig(6) +                         &
                          dsigeq_dsig3*dsig3_dsig(6)  
            endif
            !< Remove normalization of stress tensor and its derivative
            seq(i) = seq(i)*norm(i)
          enddo
        !< Shell element
        elseif (eltype == 2) then 
#include "vectorize.inc"
          do ii = 1,nindx
            i = indx(ii)
            !< Normalization of the stress tensor
            norm(i) = signxx(i)**2 + signyy(i)**2 - signxx(i)*signyy(i)        &
                                                 + three*(signxy(i)**2)
            norm(i) = max(sqrt(norm(i)),one)
            strs(i,1) = signxx(i)/norm(i)
            strs(i,2) = signyy(i)/norm(i)
            strs(i,4) = signxy(i)/norm(i)
            !< Principal stresses under plane stress condition
            center = strs(i,1) + strs(i,2)
            rootv  = sqrt((strs(i,1)-strs(i,2))**2 + (two*strs(i,4))**2)                          
            pr_strs(i,1) = half*(center + rootv)
            pr_strs(i,2) = half*(center - rootv)
            rootv  = max(rootv,em20)
            !< Equivalent stress
            seq(i) = half*((abs( pr_strs(i,1) - pr_strs(i,2) ))**nexp +        &
                           (abs( pr_strs(i,2)))**nexp +                        &
                           (abs(-pr_strs(i,1)))**nexp)
            if (seq(i) > zero) then
              seq(i) = exp((one/nexp)*log(seq(i)))
            else
              seq(i) = zero
            endif
            !< Derivatives of eq. stress
            if (seq(i) > zero) then 
              dsigeq_dsig1 = half*(seq(i)**(1-nexp))*(                         &
                    ((abs(pr_strs(i,1) - pr_strs(i,2)))**(nexp-1))*            &
                  sign(one,pr_strs(i,1) - pr_strs(i,2)) -                      &
                    ((abs(-pr_strs(i,1)))**(nexp-1))*sign(one,-pr_strs(i,1)))
              dsigeq_dsig2 = half*(seq(i)**(1-nexp))*(                         &
                    ((abs(pr_strs(i,2)))**(nexp-1))*sign(one,pr_strs(i,2)) -   &
                    ((abs(pr_strs(i,1) - pr_strs(i,2)))**(nexp-1))*            &
                  sign(one,pr_strs(i,1) - pr_strs(i,2)))
              dsig1_dsigxx = half*(one + (strs(i,1)-strs(i,2))/rootv)
              dsig1_dsigyy = half*(one - (strs(i,1)-strs(i,2))/rootv)
              dsig1_dsigxy = two*strs(i,4)/rootv
              dsig2_dsigxx = half*(one - (strs(i,1)-strs(i,2))/rootv)
              dsig2_dsigyy = half*(one + (strs(i,1)-strs(i,2))/rootv)
              dsig2_dsigxy = -two*strs(i,4)/rootv
              normxx(i) = dsigeq_dsig1*dsig1_dsigxx +                          &
                          dsigeq_dsig2*dsig2_dsigxx
              normyy(i) = dsigeq_dsig1*dsig1_dsigyy +                          &
                          dsigeq_dsig2*dsig2_dsigyy
              normzz(i) = -normxx(i) - normyy(i)
              normxy(i) = dsigeq_dsig1*dsig1_dsigxy +                          &
                          dsigeq_dsig2*dsig2_dsigxy
            else 
              normxx(i) = zero
              normyy(i) = zero
              normzz(i) = zero
              normxy(i) = zero
              normyz(i) = zero
              normzx(i) = zero
            endif
            !< Remove normalization of stress tensor and its derivative
            seq(i) = seq(i)*norm(i)
          enddo
        endif  
!
      end subroutine yield_criterion_hershey
      end module yield_criterion_hershey_mod
