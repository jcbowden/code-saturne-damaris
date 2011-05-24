!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine rijthe &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ivar   , isou   , ipp    ,                                     &
   ia     ,                                                       &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , grarox , graroy , graroz , smbr   ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! TERMES DE GRAVITE
!   POUR Rij et EPSILON
! VAR  = R11 R22 R33 R12 R13 R23 EP

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! ia(*)            ! ia ! --- ! main integer work array                        !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! grarox,y,z       ! tr ! <-- ! tableau de travail pour grad rom               !
!  (ncelet)        !    !     !                                                !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ivar   , isou   , ipp

integer          ia(*)

double precision rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision grarox(ncelet), graroy(ncelet), graroz(ncelet)
double precision smbr(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          iel
integer          ir11ip, ir22ip, ir33ip, ir12ip, ir13ip, ir23ip
integer          ieiph

double precision uns3, const, kseps
double precision prdtur, r1t, r2t, r3t
double precision g11, g22, g33, g12, g13, g23, gkks3
double precision g11p, g22p, g33p
double precision phit11, phit22, phit33, phit12, phit13, phit23
double precision aa, bb


!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ir11ip = ir11
ir22ip = ir22
ir33ip = ir33
ir12ip = ir12
ir13ip = ir13
ir23ip = ir23
ieiph  = iep

if(iscalt.gt.0.and.nscal.ge.iscalt) then
  prdtur = sigmas(iscalt)
else
  prdtur = 1.d0
endif

const = -1.5d0*cmu/prdtur
uns3  = 1.d0/3.d0

!===============================================================================
! 2. TERMES POUR RIJ :
!      ROM*VOLUME*dRij/dt =
!                     ... + (Gij - CRIJ3*(Gij-Delta ij Gkk/3))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!===============================================================================


if     (ivar.eq.ir11ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit11 = -crij3*(g11-gkks3)

    smbr(iel) = smbr(iel) + (g11+phit11)*volume(iel)

  enddo

elseif (ivar.eq.ir22ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit22 = -crij3*(g22-gkks3)

    smbr(iel) = smbr(iel) + (g22+phit22)*volume(iel)

  enddo

elseif (ivar.eq.ir33ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g11 = const*kseps*2.d0*(r1t*gx       )
    g22 = const*kseps*2.d0*(r2t*gy       )
    g33 = const*kseps*2.d0*(r3t*gz       )
    gkks3 = uns3*(g11+g22+g33)

    phit33 = -crij3*(g33-gkks3)

    smbr(iel) = smbr(iel) + (g33+phit33)*volume(iel)

  enddo

elseif (ivar.eq.ir12ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g12 = const*kseps*     (r1t*gy+r2t*gx)

    phit12 = -crij3* g12

    smbr(iel) = smbr(iel) + (g12+phit12)*volume(iel)

  enddo

elseif (ivar.eq.ir13ip) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g13 = const*kseps*     (r1t*gz+r3t*gx)

    phit13 = -crij3* g13

    smbr(iel) = smbr(iel) + (g13+phit13)*volume(iel)

  enddo

elseif (ivar.eq.ir23ip) then

  do iel = 1, ncel

    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    kseps = (rtpa(iel,ir11ip)+rtpa(iel,ir22ip)+rtpa(iel,ir33ip))  &
           /(2.d0*rtpa(iel,ieiph))

    g23 = const*kseps*(r2t*gz+r3t*gy)

    phit23 = -crij3* g23

    smbr(iel) = smbr(iel) + (g23+phit23)*volume(iel)

  enddo

!===============================================================================
! 3. TERMES POUR EPSILON :
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*(EPS/K)*MAX(0,(Gkk/2))*VOLUME
!            Avec Gij = -(1.5 CMU/PRDTUR) (K/EPS) (Rit Gj + Rjt Gi)
!                 Rit = Rik dROM/dxk (somme sur k)
!            On simplifie (EPS/K) en notant
!                GijP = -(1.5 CMU/PRDTUR)         (Rit Gj + Rjt Gi)
!      ROM*VOLUME*dEps/dt =
!                     ... + CEPS1*        MAX(0,(GkkP/2))*VOLUME
!===============================================================================


elseif (ivar.eq.ieiph ) then

  do iel = 1, ncel

    r1t = rtpa(iel,ir11ip)*grarox(iel)                            &
        + rtpa(iel,ir12ip)*graroy(iel)                            &
        + rtpa(iel,ir13ip)*graroz(iel)
    r2t = rtpa(iel,ir12ip)*grarox(iel)                            &
        + rtpa(iel,ir22ip)*graroy(iel)                            &
        + rtpa(iel,ir23ip)*graroz(iel)
    r3t = rtpa(iel,ir13ip)*grarox(iel)                            &
        + rtpa(iel,ir23ip)*graroy(iel)                            &
        + rtpa(iel,ir33ip)*graroz(iel)

    g11p = const*      2.d0*(r1t*gx)
    g22p = const*      2.d0*(r2t*gy)
    g33p = const*      2.d0*(r3t*gz)

    aa = 0.d0
    bb = 0.5d0*(g11p+g22p+g33p)
    smbr(iel) = smbr(iel) + ce1*max(aa,bb)*volume(iel)

  enddo

endif

return

end subroutine
