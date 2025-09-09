#include "SLANTWISE_OPTIONS.h"

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C !ROUTINE: SLANTWISE_DIAGNOSTICS_INIT

C !INTERFACE:
      SUBROUTINE SLANTWISE_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Initialize list of all available diagnostics

C     !USES:
      IMPLICIT NONE
C     == Global variables ===
#include "EEPARAMS.h"
#include "SIZE.h"
#include "PARAMS.h"
#include "SLANTWISE.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid ::  my Thread Id number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
C     !LOCAL VARIABLES:
C     diagNum   :: diagnostics number in the (long) list of available diag.
C     diagMate  :: diag. mate number in the (long) list of available diag.
C     diagName  :: local short name (8c) of a diagnostics
C     diagCode  :: local parser field with characteristics of the diagnostics
C              cf head of S/R DIAGNOSTICS_INIT_EARLY or DIAGNOSTICS_MAIN_INIT
C     diagUnits :: local string (16c): physical units of a diagnostic field
C     diagTitle :: local string (80c): description of field in diagnostic
      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      IF ( useDiagnostics ) THEN

        diagName  = 'SL_KYY  '
        diagTitle = 'Kyy from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_KYZ  '
        diagTitle = 'Kyz from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        
        diagName  = 'SL_KZZ  '
        diagTitle = 'Kzz from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_HN   '
        diagTitle = 'Axial height of the unstable column'
        diagUnits = 'm               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        
        diagName  = 'SL_BN   '
        diagTitle = 'Axial buoyancy gradient'
        diagUnits = 'm               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_QyKyy'
        diagTitle = 'Qy from Kyy'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_QyKyz'
        diagTitle = 'Qy from Kyz'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_QzKyz'
        diagTitle = 'Qz from Kyz'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef SLANTWISE_DEBUG
        diagName  = 'SL_mskW '
        diagTitle = 'maskW test'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_mskS '
        diagTitle = 'maskS test'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_mskC '
        diagTitle = 'maskC test'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_mskFk'
        diagTitle = 'maskFk test'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_sigY '
        diagTitle = 'drho/dy at C'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_sigZ '
        diagTitle = 'drho/dz at C'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_sY   '
        diagTitle = 'drho/dy raw'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_sZ   '
        diagTitle = 'drho/dz raw'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif
      ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
