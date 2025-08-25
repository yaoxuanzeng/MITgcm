#ifdef ALLOW_SLANTWISE

CBOP
C     !ROUTINE: SLANTWISE_INIT
C     !INTERFACE:
      SUBROUTINE SLANTWISE_INIT( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SLANTWISE_INIT
C     | o Rountine to initialze SLANTWISE runtime fields
C     *==========================================================*
C     \ev

C     !USES
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SLANTWISE.h"
#include "SLANTWISE_OPTIONS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj

C     latitude calculation
      _RL d2r, sinlat, coslat
      d2r = 3.141592653589793D0 / 180.D0

C     Radius calculation
      _RL Rbot ! Bottom radius for spherical grids
      _RL Rloc ! Local radius for spherical grids

C--   Initialize fields, loop over tiles
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C--   Set latitudes: local latitude for spherical, const lat for cartesian
        
        DO k = 1, Nr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           IF ( usingSphericalPolarGrid ) THEN
            slcon_lat(i,j,k,bi,bj) = yC(i,j,bi,bj)
           ELSE
            slcon_lat(i,j,k,bi,bj) = slcon_constlat
           ENDIF
          END DO
         END DO
        END DO

C--   Set H: column depth for spherical, const for cartisian
        DO k = 1, Nr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           IF ( usingSphericalPolarGrid ) THEN
            sinlat = SIN(slcon_lat(i,j,k,bi,bj)*d2r)
            coslat = COS(slcon_lat(i,j,k,bi,bj)*d2r)
            Rbot = rSphere - slcon_constH
            Rloc = rSphere + rC(k)
            IF ( Rloc*coslat .GT. Rbot ) THEN ! Outside the TC, Hn = Rtop * sinlat
             slcon_Hn(i,j,k,bi,bj) = rSphere*sinlat
            ELSE ! Inside the TC
             slcon_Hn(i,j,k,bi,bj) = rSphere*sinlat - Rbot*sinlat
            END IF
           ELSE
            slcon_Hn(i,j,k,bi,bj) = slcon_constH/SIN(slcon_constlat*d2r)
           END IF
          END DO
         END DO
        END DO

C--   Set tensor eliments
        DO k = 1, Nr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           slcon_Kyy(i,j,k,bi,bj) = 0.0d0
           slcon_Kyz(i,j,k,bi,bj) = 0.0d0
           slcon_Kzz(i,j,k,bi,bj) = 0.0d0
          END DO
         END DO
        END DO

       END DO
      END DO

      RETURN
      END

#endif /* ALLOW_SLANTWISE */

