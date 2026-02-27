#include "SLANTWISE_OPTIONS.h"

CBOP
C     !ROUTINE: SLANTWISE_RTRANSPORT
C     !INTERFACE:
      SUBROUTINE SLANTWISE_RTRANSPORT(
     I     trIdentity, bi, bj, k,
     I     iMin, iMax, jMin, jMax,
     I     maskFk, Tracer,
     U     df,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | o SUBROUTINE SLANTWISE_RTRANSPORT
C     |   Add horizontal r transport terms from slantwise
C     |   convection parameterization.
C     *==========================================================*
C     \ev

      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SLANTWISE.h"
#include "SLANTWISE_OPTIONS.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     trIdentity   :: tracer Id number
C     bi, bj       :: current tile indices
C     k            :: current level index
C     iMin,iMax    :: Range of 1rst index where results will be set
C     jMin,jMax    :: Range of 2nd  index where results will be set
C     maskUp       :: 2-D array for mask at W points
C     Tracer       :: 3D Tracer field
C     df           :: Diffusive flux component work array.
      
      INTEGER trIdentity, bi, bj, k
      INTEGER iMin, iMax, jMin, jMax
      _RS maskFk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Tracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL df    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid

CEOP

#ifdef ALLOW_SLANTWISE

C     !LOCAL VARIABLES
      INTEGER i,j,km1 ! Local indices
C       _RL dCdz_W
       _RL dCdy_W ! Tracer gradients at W grids
#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
      _RL dTdy_W, Tflux
      _RL dTdz_W
#endif
#ifdef SLANTWISE_SFLUX_DIAGNOSTICS
      _RL dSdy_W, Sflux
      _RL dSdz_W
#endif

      km1  = MAX(k-1,1)

C--   Loop for each tile
      DO j = jMin, jMax
       DO i = iMin, iMax

C-    Calculate the tracer gradients at the W grids
C     Calculate dCdz_W at the W grids 
C        dCdz_W = ( Tracer(i,j,km1) - Tracer(i,j,k) ) * maskC(i,j,k,bi,bj) * recip_drF(k)
C     Calculate dCdy_W at the W grids, copied from gmredi
        IF (k.GT.1) THEN
         dCdy_W = 0.5 _d 0*(
     &     0.5 _d 0*recip_deepFacC(k)
     &        *( _maskS(i,j+1,k,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( Tracer(i,j+1,k) - Tracer(i,j,k) )
     &         + _maskS(i, j ,k,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( Tracer(i,j,k) - Tracer(i,j-1,k) )
     &         )
     &    +0.5 _d 0*recip_deepFacC(k-1)
     &        *( _maskS(i,j+1,k-1,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( Tracer(i,j+1,k-1) - Tracer(i,j,k-1) )
     &         + _maskS(i, j ,k-1,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( Tracer(i,j,k-1) - Tracer(i,j-1,k-1) )
     &         )          )
        END IF

C-    Calculate the flux
         df(i,j) = df(i,j) -_rA(i,j,bi,bj)*deepFac2F(k)
     &            *slcon_Kwy(i,j,k,bi,bj)*dCdy_W
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
C     &           -_rA(i,j,bi,bj)*deepFac2F(k)
C     &            *slcon_Kwz(i,j,k,bi,bj)*dCdz_W
C     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)

#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
        dTdz_W =  ( theta(i,j,km1,bi,bj)
     &        - theta(i,j,k,bi,bj) )
     &        * maskC(i,j,k,bi,bj) * recip_drF(k)
        IF (k.GT.1) THEN
         dTdy_W = 0.5 _d 0*(
     &     0.5 _d 0*recip_deepFacC(k)
     &        *( _maskS(i,j+1,k,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( theta(i,j+1,k,bi,bj) - theta(i,j,k,bi,bj) )
     &         + _maskS(i, j ,k,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( theta(i,j,k,bi,bj) - theta(i,j-1,k,bi,bj) )
     &         )
     &    +0.5 _d 0*recip_deepFacC(k-1)
     &        *( _maskS(i,j+1,k-1,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( theta(i,j+1,k-1,bi,bj) - theta(i,j,k-1,bi,bj) )
     &         + _maskS(i, j ,k-1,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( theta(i,j,k-1,bi,bj) - theta(i,j-1,k-1,bi,bj) )
     &         )          )
        END IF
         slcon_Tz(i,j,k,bi,bj) = -_rA(i,j,bi,bj)*deepFac2F(k)
     &            *slcon_Kwy(i,j,k,bi,bj)*dTdy_W
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
     &           -_rA(i,j,bi,bj)*deepFac2F(k)
     &            *slcon_Kwz(i,j,k,bi,bj)*dTdz_W
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
#endif

#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
        dSdz_W =  ( salt(i,j,km1,bi,bj)
     &        - salt(i,j,k,bi,bj) )
     &        * maskC(i,j,k,bi,bj) * recip_drF(k)
        IF (k.GT.1) THEN
         dSdy_W = 0.5 _d 0*(
     &     0.5 _d 0*recip_deepFacC(k)
     &        *( _maskS(i,j+1,k,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( salt(i,j+1,k,bi,bj) - salt(i,j,k,bi,bj) )
     &         + _maskS(i, j ,k,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( salt(i,j,k,bi,bj) - salt(i,j-1,k,bi,bj) )
     &         )
     &    +0.5 _d 0*recip_deepFacC(k-1)
     &        *( _maskS(i,j+1,k-1,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( salt(i,j+1,k-1,bi,bj) - salt(i,j,k-1,bi,bj) )
     &         + _maskS(i, j ,k-1,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( salt(i,j,k-1,bi,bj) - salt(i,j-1,k-1,bi,bj) )
     &         )          )
        END IF
         slcon_Sz(i,j,k,bi,bj) = -_rA(i,j,bi,bj)*deepFac2F(k)
     &            *slcon_Kwy(i,j,k,bi,bj)*dSdy_W
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
     &           -_rA(i,j,bi,bj)*deepFac2F(k)
     &            *slcon_Kwz(i,j,k,bi,bj)*dSdz_W
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
#endif

       END DO
      END DO

#endif /* ALLOW_SLANTWISE */

      RETURN
      END


