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

         CALL DIAGNOSTICS_FILL( slcon_Kyy,'SL_KYY  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Kyz,'SL_KYZ  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Kzz,'SL_KZZ  ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_Hn, 'SL_HN   ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_bn, 'SL_BN   ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_qykyy, 'SL_QyKyy',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_qykyz, 'SL_QyKyz',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_qzkyz, 'SL_QzKyz',
     &                          0,Nr, 1, bi, bj, myThid )
#ifdef SLANTWISE_DEBUG
         CALL DIAGNOSTICS_FILL( slcon_maskS, 'SL_mskS ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_maskW, 'SL_mskW ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_maskC, 'SL_mskC ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_maskFk, 'SL_mskFk',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_sigmay, 'SL_sigY ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_sigmaz, 'SL_sigZ ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_sy, 'SL_sY   ',
     &                          0,Nr, 1, bi, bj, myThid )
         CALL DIAGNOSTICS_FILL( slcon_sz, 'SL_sZ   ',
     &                          0,Nr, 1, bi, bj, myThid )
#endif
        END DO
       END DO

      ENDIF
#endif /* ALLOW_DIAGNOSTICS */

#endif /* ALLOW_SLANTWISE */

      RETURN
      END
