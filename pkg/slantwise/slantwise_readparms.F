#ifdef ALLOW_SLANTWISE

CBOP
C     !ROUTINE: SLANTWISE_READPARMS
C     !INTERFACE:
      SUBROUTINE SLANTWISE_READPARMS( myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SLANTWISE_READPARMS
C     | o Routine to initialize slantwise variables and constants
C     *==========================================================*
C     | Initialize slantwise parameters, read in data.slantwise
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SLANTWISE.h"
#include "SLANTWISE_OPTIONS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
      INTEGER myThid

C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf     :: Informational/error message buffer
C     iUnit      :: Work variable for IO unit number
      CHARACTER*(MAX_LEN_MBUF) msgBuf
      INTEGER iUnit

CEOP


C--   Namelist
      NAMELIST /SLANTWISE_PARM01/
     &                          slcon_on, slcon_useconstH,
     &                          slcon_constlat, slcon_constH,
     &                          slcon_cn, slcon_cr

C     default values
      slcon_on                = .FALSE.
      slcon_useconstH         = .TRUE.
      slcon_constlat          = 45.d0
      slcon_constH            = 6.0d4
      slcon_cn                = 6.60d-1
      slcon_cr                = 3.26d-1
        
      WRITE(msgBuf,'(A)') ' SLANTWISE_READPARMS: opening data.slantwise'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit, SQUEEZE_RIGHT,1)
      CALL OPEN_COPY_DATA_FILE(
     I                        'data.slantwise','SLANTWISE_READPARMS',
     O                        iUnit,
     I                        myThid )

      READ(UNIT=iUnit,NML=SLANTWISE_PARM01)
      WRITE(msgBuf,'(A)') ' SLANTWISE_READPARMS: finish data.slantwise'
      CLOSE(iUnit)
      
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit,
     &                    SQUEEZE_RIGHT , 1)

      RETURN
      END

#endif /* ALLOW_SLANTWISE */

