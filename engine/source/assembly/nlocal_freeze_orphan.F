Copyright>        OpenRadioss
Copyright>        Copyright (C) 2026 Siemens
Copyright>
Copyright>        This program is free software: you can redistribute it and/or modify
Copyright>        it under the terms of the GNU Affero General Public License as published by
Copyright>        the Free Software Foundation, either version 3 of the License, or
Copyright>        (at your option) any later version.
Copyright>
Copyright>        This program is distributed in the hope that it will be useful,
Copyright>        but WITHOUT ANY WARRANTY; without even the implied warranty of
Copyright>        MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
Copyright>        GNU Affero General Public License for more details.
Copyright>
Copyright>        You should have received a copy of the GNU Affero General Public License
Copyright>        along with this program.  If not, see <https://www.gnu.org/licenses/>.
Copyright>
Copyright>
Copyright>        Commercial Alternative: Simcenter Radioss Software
Copyright>
Copyright>        As an alternative to this open-source version, Siemens also offers Simcenter(TM) Radioss(R)
Copyright>        software under a commercial license.  Contact Siemens to discuss further if the
Copyright>        commercial version may interest you: 
Copyright>        https://www.siemens.com/en-us/products/simcenter/mechanical-simulation/radioss/.
!||====================================================================
!||    nlocal_freeze_orphan   ../engine/source/assembly/nlocal_freeze_orphan.F
!||--- called by ------------------------------------------------------
!||    resol                  ../engine/source/engine/resol.F
!||--- uses       -----------------------------------------------------
!||    constant_mod           ../common_source/modules/constant_mod.F
!||    nlocal_reg_mod         ../common_source/modules/nlocal_reg_mod.F
!||====================================================================
      SUBROUTINE NLOCAL_FREEZE_ORPHAN(NLOC_DMG, NODFT   ,NODLT)
C-----------------------------------------------
C   D e s c r i p t i o n
C-----------------------------------------------
C   Option (A) - zero-flux (homogeneous Neumann) crack boundary.
C
C   Freezes "orphan" non-local degrees of freedom: d.o.fs whose ALL
C   supporting elements have eroded. Such a d.o.f no longer receives any
C   non-local internal force (with the absorbing dashpot removed, an
C   eroded element contributes exactly zero), so it becomes a spurious
C   zero-energy mode that would drift indefinitely. It is detected from
C   the assembled non-local nodal stiffness STIFNL, which stays at the
C   EM20 seed floor when no live element supports the d.o.f. The drift is
C   arrested by zeroing the non-local velocity (VNL/VNL_OLD); since the
C   assembled force is zero, NLOCAL_ACC/VEL/INCR then hold UNL frozen at
C   its last valid value.
C
C   This routine MUST be called after STIFNL is fully assembled AND (for
C   MPI) frontier-exchanged, i.e. right after NLOCAL_ACC. It is the unique
C   consumer of the post-exchange STIFNL, so it also resets STIFNL to the
C   EM20 seed for the next cycle's element accumulation (this reset used
C   to live at the end of NLOCAL_DTNODA).
C-----------------------------------------------
C   M o d u l e s
C-----------------------------------------------
      USE NLOCAL_REG_MOD
      USE CONSTANT_MOD, ONLY : ZERO, EM20
C-----------------------------------------------
C   I m p l i c i t   T y p e s
C-----------------------------------------------
#include      "implicit_f.inc"
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER NODFT, NODLT
      TYPE(NLOCAL_STR_), TARGET :: NLOC_DMG
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER N,I,NPOS,NDDL
      my_real, POINTER, DIMENSION(:) ::
     .  STIFN,VEL,VEL0
      ! Stiffness threshold (~1000x the EM20 non-local stiffness floor)
      ! below which a d.o.f is considered orphaned. A live element adds
      ! O(VOL/NTN) >> this; an all-eroded d.o.f only carries the EM20 seed
      ! plus a few EM20 dead-element contributions (< tol for < ~1000 dead
      ! elements sharing a node). For extremely fine SI/meter-unit meshes
      ! this floor can be lowered if needed.
      my_real, PARAMETER :: STIF_ORPHAN_TOL = 1.0D-17
C=======================================================================
      IF (NLOC_DMG%L_NLOC <= 0) RETURN
C
      ! Recovering assembled nodal stiffness and non-local velocities
      STIFN => NLOC_DMG%STIFNL(1:NLOC_DMG%L_NLOC,1)
      VEL   => NLOC_DMG%VNL(1:NLOC_DMG%L_NLOC)
      VEL0  => NLOC_DMG%VNL_OLD(1:NLOC_DMG%L_NLOC)
C
      ! Loop over the nodes of the subset
      DO N = NODFT,NODLT
        ! Position of the first additional d.o.f
        NPOS = NLOC_DMG%POSI(N)
        ! Number of additional degrees of freedom
        NDDL = NLOC_DMG%POSI(N+1) - NPOS
        ! Freeze orphan d.o.fs and reset the stiffness seed
        DO I = NPOS,NPOS+NDDL-1
          IF (STIFN(I) <= STIF_ORPHAN_TOL) THEN
            VEL(I)  = ZERO
            VEL0(I) = ZERO
          ENDIF
          ! Re-seed the EM20 stiffness floor for next cycle's accumulation
          STIFN(I) = EM20
        ENDDO
      ENDDO
C
      RETURN
      END
