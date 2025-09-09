#include "SLANTWISE_OPTIONS.h"

CBOP
C     !ROUTINE: SLANTWISE_INIT_VARIA
C     !INTERFACE:
      SUBROUTINE SLANTWISE_INIT_VARIA( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SLANTWISE_INIT_VARIA
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

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_SLANTWISE

C     !LOCAL VARIABLES:
      INTEGER i,j,k,bi,bj

C--   Initialize fields, loop over tiles
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)

C--   Set tensor eliments
        DO k = 1, Nr
         DO j = 1-OLy, sNy+OLy
          DO i = 1-OLx, sNx+OLx
           slcon_Kyy(i,j,k,bi,bj) = 0.0d0
           slcon_Kyz(i,j,k,bi,bj) = 0.0d0
           slcon_Kzz(i,j,k,bi,bj) = 0.0d0
           slcon_bn(i,j,k,bi,bj) = 0.0d0
          END DO
         END DO
        END DO

       END DO
      END DO

#endif /* ALLOW_SLANTWISE */

      RETURN
      END


