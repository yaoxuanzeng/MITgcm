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

        diagName  = 'SL_KVY  '
        diagTitle = 'Kvy from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'VV P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_KVZ  '
        diagTitle = 'Kvz from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'VV      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        
        diagName  = 'SL_KWY  '
        diagTitle = 'Kwy from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'VM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_KWZ  '
        diagTitle = 'Kwz from slantwise convection'
        diagUnits = 'm^2/s           '
        diagCode  = 'SM P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_HNV  '
        diagTitle = 'Axial height of the unstable column at V points'
        diagUnits = 'm               '
        diagCode  = 'VV      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
        
        diagName  = 'SL_HNW  '
        diagTitle = 'Axial height of the unstable column at W points'
        diagUnits = 'm               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
        diagName  = 'SL_Ty   '
        diagTitle = 'Meridional T flux'
        diagUnits = '1               '
        diagCode  = 'VV P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_Tz   '
        diagTitle = 'Vertical T flux'
        diagUnits = '1               '
        diagCode  = 'VM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

#ifdef SLANTWISE_SFLUX_DIAGNOSTICS
        diagName  = 'SL_Sy   '
        diagTitle = 'Meridional S flux'
        diagUnits = '1               '
        diagCode  = 'VV P    MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_Sz   '
        diagTitle = 'Vertical S flux'
        diagUnits = '1               '
        diagCode  = 'VM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

#ifdef SLANTWISE_DEBUG
        diagName  = 'SL_bug1 '
        diagTitle = 'Debug variable 1'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug2 '
        diagTitle = 'Debug variable 2'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug3 '
        diagTitle = 'Debug variable 3'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug4 '
        diagTitle = 'Debug variable 4'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug5 '
        diagTitle = 'Debug variable 5'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug6 '
        diagTitle = 'Debug variable 6'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug7 '
        diagTitle = 'Debug variable 7'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug8 '
        diagTitle = 'Debug variable 8'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

        diagName  = 'SL_bug9 '
        diagTitle = 'Debug variable 9'
        diagUnits = '1               '
        diagCode  = 'SM      MR      '
        CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#endif

      ENDIF

#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
