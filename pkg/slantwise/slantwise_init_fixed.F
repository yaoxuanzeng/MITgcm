#include "SLANTWISE_OPTIONS.h"

CBOP
C     !ROUTINE: SLANTWISE_INIT_FIXED
C     !INTERFACE:
      SUBROUTINE SLANTWISE_INIT_FIXED( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SLANTWISE_INIT_FIXED
C     | o Rountine to initialze SLANTWISE fixed fields
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

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj

C     latitude calculation
      _RL d2r, sinlat, coslat
      _RL Rbot, Rloc ! Bottom and local radius for spherical grids

      d2r = 3.141592653589793D0 / 180.D0

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
            slcon_Hn(i,j,k,bi,bj)=slcon_constH/SIN(slcon_constlat*d2r)
           END IF
          END DO
         END DO
        END DO
      
#ifdef SLANTWISE_DEBUG
        DO k = 1, Nr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           slcon_maskS(i,j,k,bi,bj) = 1
           slcon_maskW(i,j,k,bi,bj) = 1
           slcon_maskC(i,j,k,bi,bj) = 1
           slcon_maskFk(i,j,k,bi,bj) = 1
           slcon_sigmay(i,j,k,bi,bj) = 0
           slcon_sigmaz(i,j,k,bi,bj) = 0
          END DO
         END DO
        END DO
#endif

       END DO
      END DO

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL SLANTWISE_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

        RETURN
        END

