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

subroutine lageqp &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , propce , propfa , propfb ,                            &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt ,                                              &
   fmala  , fmalb  ,                                              &
   ul     , vl     , wl     , alphal , phia   , phi    ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!          RESOLUTION D'UNE EQUATION DE POISSON

!            div[ALPHA grad(PHI)] = div(ALPHA <Up>)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! smbrs(ncelet     ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! fmala(nfac)      ! tr ! --- ! flux de masse au faces internes                !
! fmalb(nfabor)    ! tr ! --- ! flux de masse au faces de bord                 !
! ul,vl,wl         ! tr ! <-- ! vitesse lagrangien                             !
! (ncelet)         !    !     !                                                !
! alphal           ! tr ! <-- ! taux de presence                               !
! (ncelet)         !    !     !                                                !
! phi , phia       ! tr ! --> ! terme de correction en n et n-1                !
! (ncelet)         !    !     !                                                !
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
use cstnum
use pointe
use lagpar
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ia(*)

double precision ul(ncelet), vl(ncelet), wl(ncelet)
double precision phia(ncelet), phi(ncelet), alphal(ncelet)
double precision fmala(nfac) , fmalb(nfabor)
double precision dt(ncelet)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)
double precision ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra , ifinia , ifinra
integer          idtva0, ivar
integer          ifac, iel
integer          ipp
integer          nswrgp, imligp, iwarnp , iescap
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp
integer          imgrp, ncymxp, nitmfp
integer          icoefax,icoefay,icoefaz,icefap
integer          icoefbx,icoefby,icoefbz,icefbp

double precision epsrgp, climgp, extrap, blencp, epsilp, epsrsp
double precision relaxp, thetap

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, w2, w3

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))

idebia = idbia0
idebra = idbra0

CHAINE = 'Correction pression'
write(nfecra,1000) chaine(1:19)

!===============================================================================
! 2. TERMES SOURCES
!===============================================================================

! --> Initialisation

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo
do iel = 1, ncel
  phi(iel)  = 0.d0
  phia(iel) = 0.d0
enddo

!     "VITESSE" DE DIFFUSION FACE

  call viscfa                                                     &
  !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   alphal ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

! CALCUL  de div(Alpha Up) avant correction

do iel = 1, ncel
  w1(iel) = -ul(iel)*alphal(iel)
  w2(iel) = -vl(iel)*alphal(iel)
  w3(iel) = -wl(iel)*alphal(iel)
enddo

! --> Calcul du gradient de W1
!     ========================

!       On alloue localement 6 tableaux de NFABOR pour le calcul
!         de COEFA et COEFB de W1,W2,W3

icoefax = idebra
icoefbx = icoefax + nfabor
icoefay = icoefbx + nfabor
icoefby = icoefay + nfabor
icoefaz = icoefby + nfabor
icoefbz = icoefaz + nfabor
ifinra  = icoefbz + nfabor
call rasize ('lageqp',ifinra)
!==========

do ifac = 1, nfabor
  iel = ifabor(ifac)

  ra(icoefax+ifac-1) = w1(iel)
  ra(icoefbx+ifac-1) = zero

  ra(icoefay+ifac-1) = w2(iel)
  ra(icoefby+ifac-1) = zero

  ra(icoefaz+ifac-1) = w3(iel)
  ra(icoefbz+ifac-1) = zero

enddo

call diverv                                                       &
!==========
 ( idebia , ifinra ,                                              &
   nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     ,                                                       &
   smbrs  , w1     , w2     , w3     ,                            &
   ra(icoefax) , ra(icoefay) , ra(icoefaz) ,                      &
   ra(icoefbx) , ra(icoefby) , ra(icoefbz) ,                      &
   ra     )

!      On libere la place dans RA

ifinra = idebra

! --> Conditions aux limites sur PHI
!     ==============================

!       On alloue localement 2 tableaux de NFABOR pour le calcul
!         de COEFA et COEFB de PHI

ifinia = idebia
icefap = idebra
icefbp = icefap + nfabor
ifinra = icefbp + nfabor
call rasize ('lageqp',ifinra)
!==========

do ifac = 1, nfabor
  iel = ifabor(ifac)

  if ( ia(iitypf+ifac-1) .eq. ientre ) then

!      Flux Nul

    ra(icefap+ifac-1) = zero
    ra(icefbp+ifac-1) = 1.d0

  else if ( ia(iitypf+ifac-1) .eq. iparoi) then

!      FLux nul

    ra(icefap+ifac-1) = zero
    ra(icefbp+ifac-1) = 1.d0

  else if ( ia(iitypf+ifac-1) .eq. iparug) then

!      FLux nul

    ra(icefap+ifac-1) = zero
    ra(icefbp+ifac-1) = 1.d0

  else if ( ia(iitypf+ifac-1) .eq. isymet) then

!      FLux nul

    ra(icefap+ifac-1) = zero
    ra(icefbp+ifac-1) = 1.d0

  else if ( ia(iitypf+ifac-1) .eq. isolib ) then

!      Valeur Imposee

    ra(icefap+ifac-1) = phia(iel)
    ra(icefbp+ifac-1) = zero

  else
    write(nfecra,1100) ia(iitypf+ifac-1)
    call csexit (1)
!              ======
  endif

enddo

!===============================================================================
! 3. RESOLUTION
!===============================================================================

! Pas de stationnaire
idtva0 = 0
! Pas de terme de convection
iconvp = 0
! Diffusion
idiffp = 1
! Methode de resolution : Gradient conjugue (pas de convection)
ireslp = 0
! Valeur par defaut
ndircp = 1
nitmap = 1000
nswrsp = 2
nswrgp = 10000
imligp = 1
ircflp = 1
ischcp = 1
isstpp = 0
imgrp  = 1
ncymxp = 100
nitmfp = 100
iwarnp = 10
blencp = 0.d0
epsilp = 1.d-8
epsrsp = 1.d-8
epsrgp = 1.d-5
climgp = 1.5d0
extrap = 0.d0
relaxp = 1.d0
iescap = 0

ipp  = 1
NOMVAR(IPP) = 'PoissonL'

!  IVAR = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

ivar = 0

! On annule les flux de masse

do ifac = 1,nfac
  fmala(ifac) = zero
enddo

do ifac = 1,nfabor
  fmalb(ifac) = zero
enddo

! Dans le cas d'un theta-schema on met theta = 1 (ordre 1)

thetap = 1.0d0

call codits                                                       &
!==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  ,                                              &
   idtva0 , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ia     ,                                                       &
   phia   , phia   , ra(icefap) , ra(icefbp) ,                    &
             ra(icefap) , ra(icefbp) ,                            &
             fmala       , fmalb       ,                          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbrs  , phi    ,                                     &
   rvoid  ,                                                       &
   ra     )


! Free memory
deallocate(w1, w2, w3)

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION POUR LA VARIABLE ',A19                       ,/,&
'      ---------------------------                            ',/)

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========                                               ',/,&
'@    ERREUR A LA RESOLUTION DE L''EQUATION DE POISSON :      ',/,&
'@      CONDITIONS AUX LIMITES SUR PHI NON PREVUES (LAGEQP).  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
