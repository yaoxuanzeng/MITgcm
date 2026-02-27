#include "SLANTWISE_OPTIONS.h"

CBOP
      SUBROUTINE SLANTWISE_CHECK( myThid )
      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "SLANTWISE_OPTIONS.h"
      INTEGER myThid
      CHARACTER*(MAX_LEN_MBUF) msgBuf
CEOP

#ifdef ALLOW_SLANTWISE

      WRITE(msgBuf,'(A)') ' SLANTWISE_CHECK: package enabled.'
      CALL PRINT_MESSAGE( msgBuf, standardMessageUnit, SQUEEZE_RIGHT,1)
 
#endif /* ALLOW_SLANTWISE */

      RETURN
      END


