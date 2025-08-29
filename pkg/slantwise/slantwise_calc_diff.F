#include "SLANTWISE_OPTIONS.h"

      SUBROUTINE SLANTWISE_CALC_DIFF(
     I        bi,bj,iMin,iMax,jMin,jMax,kArg,kSize,
     U        KappaRx,
     I        tracerIdentity,myThid)

C     /==========================================================\
C     | SUBROUTINE SLANTWISE_CALC_DIFF                           |
C     | o Add contribution to net diffusivity from slantwise     |
C     \==========================================================/
      IMPLICIT NONE

C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "GAD.h"
#include "SLANTWISE.h"

C     == Routine arguments ==
C     bi, bj,   :: tile indices
C     iMin,iMax :: Range of points for which calculation is done
C     jMin,jMax :: Range of points for which calculation is done
C     kArg      :: = 0 -> do the k-loop here and treat all levels
C                  > 0 -> k-loop is done outside and treat only level
C                  k=kArg
C     kSize     :: 3rd Dimension of the vertical diffusivity array
C     KappaRx
C     KappaRx   :: vertical diffusivity array 
C     myThid    :: Instance number for this innvocation of
C     SLANTWISE_CALC_DIFF
C
      INTEGER bi,bj,iMin,iMax,jMin,jMax,kArg,kSize
      _RL KappaRx(1-Olx:sNx+Olx,1-Oly:sNy+Oly,kSize)
      INTEGER tracerIdentity
      INTEGER myThid

#ifdef ALLOW_SLANTWISE

C     == Local variables ==
C     i,j,k     :: Loop counters
      INTEGER i,j,k

C--   Add Vertical diffusivity contribution from slantwise
      IF ( kArg .EQ. 0 ) THEN
C-    do all levels :
        DO k=1,MIN(Nr,kSize)
         DO j=jMin,jMax
          DO i=iMin,iMax
           KappaRx(i,j,k) = KappaRx(i,j,k)+slcon_Kzz(i,j,k,bi,bj)
     &           *maskInC(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
      ELSE
C-    do level k=kArg only :
         k = MIN(kArg,kSize)
         DO j=jMin,jMax
          DO i=iMin,iMax
           KappaRx(i,j,k) = KappaRx(i,j,k)+slcon_Kzz(i,j,kArg,bi,bj)
     &           *maskInC(i,j,bi,bj)
          ENDDO
         ENDDO
      ENDIF

#endif /* ALLOW_SLANTWISE */

      RETURN
      END




