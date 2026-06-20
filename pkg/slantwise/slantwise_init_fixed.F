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
      _RL d2r, Rbot
      _RL sinlatw, coslatw, sinlatv, coslatv
      _RL ylocv, ylocw, Rlocv, Rlocw ! Bottom and local radius for spherical grids

      d2r = 3.141592653589793D0 / 180.D0

C--   Initialize H fields, loop over tiles
      IF (( slconHnvFile .NE. ' ') .AND. ( slconHnwFile .NE. ' ')) THEN
       _BARRIER
       CALL READ_REC_XYZ_RS( slconHnvFile, slcon_Hnv,
     &                      1, 0, myThid )
       CALL READ_REC_XYZ_RS( slconHnwFile, slcon_Hnw,
     &                      1, 0, myThid )
      ELSE
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)

C--    Set H: column depth for spherical, const for cartisian
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            IF ( usingSphericalPolarGrid ) THEN
             sinlatv = SIN(yG(i,j,bi,bj)*d2r)
             coslatv = COS(yG(i,j,bi,bj)*d2r)
             sinlatw = SIN(yC(i,j,bi,bj)*d2r)
             coslatw = COS(yC(i,j,bi,bj)*d2r)

             Rbot = rSphere - slcon_constH
             Rlocv = rSphere + rC(k)
             Rlocw = rSphere + rF(k)

             ylocv = Rlocv * coslatv
             ylocw = Rlocw * coslatw

             IF ( Rlocv*coslatv .GT. Rbot ) THEN ! Outside the TC, Hn = Rtop * sinlat
C              slcon_Hnv(i,j,k,bi,bj)=SQRT(rSphere*rSphere-ylocv*ylocv)
              slcon_Hnv(i,j,k,bi,bj)=SQRT(rSphere*rSphere-Rbot*Rbot)
             ELSE ! Inside the TC
              slcon_Hnv(i,j,k,bi,bj)=SQRT(rSphere*rSphere-ylocv*ylocv)
     &                              -SQRT(Rbot*Rbot-ylocv*ylocv)
             END IF
             IF ( Rlocw*coslatw .GT. Rbot ) THEN ! Outside the TC, Hn = Rtop * sinlat
C              slcon_Hnw(i,j,k,bi,bj)=SQRT(rSphere*rSphere-ylocw*ylocw)
              slcon_Hnw(i,j,k,bi,bj)=SQRT(rSphere*rSphere-Rbot*Rbot)
             ELSE ! Inside the TC
              slcon_Hnw(i,j,k,bi,bj) = SQRT(rSphere*rSphere-ylocw*ylocw)
     &                               - SQRT(Rbot*Rbot-ylocw*ylocw)
             END IF

            ELSE
             slcon_Hnv(i,j,k,bi,bj)=slcon_constH/SIN(slcon_constlat*d2r)
             slcon_Hnw(i,j,k,bi,bj)=slcon_constH/SIN(slcon_constlat*d2r)
            END IF
           END DO
          END DO
         END DO

        END DO
       END DO
      ENDIF

C--   Initialize H fields, loop over tiles
      IF ((slconvweightFile.NE.' ').AND.(slconwweightFile.NE.' ')) THEN
       _BARRIER
       CALL READ_REC_XYZ_RS( slconvweightFile, slcon_vweight,
     &                      1, 0, myThid )
       CALL READ_REC_XYZ_RS( slconwweightFile, slcon_wweight,
     &                      1, 0, myThid )
      ELSE
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)
         DO k = 1, Nr
          DO j = 1-OLy, sNy+OLy
           DO i = 1-OLx, sNx+OLx
            slcon_vweight(i,j,k,bi,bj) = 1.D0
            slcon_wweight(i,j,k,bi,bj) = 1.D0
           END DO
          END DO
         END DO

        END DO
       END DO
      ENDIF

#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
        CALL SLANTWISE_DIAGNOSTICS_INIT( myThid )
      ENDIF
#endif

        RETURN
        END

