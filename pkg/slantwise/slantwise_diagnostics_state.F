#include "SLANTWISE_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SLANTWISE_DIAGNOSTICS_STATE

C !INTERFACE:
      SUBROUTINE SLANTWISE_DIAGNOSTICS_STATE(
     I                     myTime, myIter, myThid )

C     !DESCRIPTION:
C     Fill-in diagnostics array for SLANTWISE state variables

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SLANTWISE.h"

C     !INPUT PARAMETERS:
C     myTime   :: Current time of simulation ( s )
C     myIter   :: Current iteration number in simulation
C     myThid   :: my Thread Id number
C     bi, bj :: tile indices
      INTEGER myTime, myIter, myThid
      INTEGER bi, bj
CEOP

#ifdef ALLOW_SLANTWISE
#ifdef ALLOW_DIAGNOSTICS
      IF ( useDiagnostics ) THEN
       DO bj = myByLo(myThid), myByHi(myThid)
        DO bi = myBxLo(myThid), myBxHi(myThid)

         CALL DIAGNOSTICS_FILL( slcon_Kvy,'SL_KVY  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Kvz,'SL_KVZ  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Kwy,'SL_KWY  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Kwz,'SL_KWZ  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Hnv, 'SL_HNV  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Hnw, 'SL_HNW  ',
     &                          0,Nr, 1, bi, bj, myThid )

#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( slcon_Ty, 'SL_Ty   ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Tz, 'SL_Tz   ',
     &                          0,Nr, 1, bi, bj, myThid )
#endif
#ifdef SLANTWISE_SFLUX_DIAGNOSTICS
         CALL DIAGNOSTICS_FILL( slcon_Sy, 'SL_Sy   ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Sz, 'SL_Sz   ',
     &                          0,Nr, 1, bi, bj, myThid )
#endif

#ifdef SLANTWISE_DEBUG
         CALL DIAGNOSTICS_FILL( slcon_debug1, 'SL_bug1 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug2, 'SL_bug2 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug3, 'SL_bug3 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug4, 'SL_bug4 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug5, 'SL_bug5 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug6, 'SL_bug6 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug7, 'SL_bug7 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug8, 'SL_bug8 ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_debug9, 'SL_bug9 ',
     &                          0,Nr, 1, bi, bj, myThid )
#endif

        END DO
       END DO

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_SLANTWISE */

      RETURN
      END
