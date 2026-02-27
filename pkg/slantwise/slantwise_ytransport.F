#include "SLANTWISE_OPTIONS.h"

CBOP
C     !ROUTINE: SLANTWISE_YTRANSPORT
C     !INTERFACE:
      SUBROUTINE SLANTWISE_YTRANSPORT(
     I     trIdentity, bi, bj, k,
     I     iMin, iMax, jMin, jMax,
     I     yA, maskFk, Tracer,
     U     df,
     I     myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | o SUBROUTINE SLANTWISE_YTRANSPORT
C     |   Add horizontal y transport terms from slantwise
C     |   convection parameterization.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == GLobal variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SLANTWISE.h"
#include "SLANTWISE_OPTIONS.h"
#include "GAD.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     trIdentity :: tracer Id number
C     bi, bj     :: current tile indices
C     k          :: current level index
C     iMin,iMax  :: Range of 1rst index where results will be set
C     jMin,jMax  :: Range of 2nd  index where results will be set
C     yA         :: Area of Y face
C     maskFk     :: 2-D mask for vertical interface k (between level k-1 & k)
C     Tracer     :: 3D Tracer field
C     df         :: Diffusive flux component work array.
C     myThid     :: my Thread Id number
      INTEGER trIdentity
      INTEGER bi, bj, k
      INTEGER iMin, iMax, jMin, jMax
      _RS yA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskFk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Tracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL df    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid

CEOP

#ifdef ALLOW_SLANTWISE

C     !LOCAL VARIABLES:
      INTEGER i,j,kp1,km1 ! local indices
      _RL maskp1 ! mask
      _RL dCdz_S, dCdy_S ! Tracer gradients at the southern boundary / W grids
#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
      _RL dTdz_S, dTdy_S, Tflux
#endif
#ifdef SLANTWISE_SFLUX_DIAGNOSTICS
      _RL dSdz_S, dSdy_S, Sflux
#endif

      km1 = MAX(k-1,1)
      kp1 = MIN(k+1,Nr)

C--   Loop for each tile
      DO j = jMin, jMax
       DO i = iMin, iMax

C-    Calculate the tracer gradients at the southern boundary
C     Calculate dC/dy at the southern boundary
        dCdy_S = (Tracer(i,j,k)-Tracer(i,j-1,k))
     &          *_recip_dyC(i,j,bi,bj)*recip_deepFacC(k)

C     Calculate dC/dz at W grids, copied from gmredi
        maskp1 = 1. _d 0
        IF (k.GE.Nr) maskp1 = 0. _d 0
        dCdz_S = 0.5 _d 0*(
     &    +0.5 _d 0*recip_drC(k)*
     &        ( maskFk(i,j-1)*
     &            ( Tracer(i,j-1,km1)-Tracer(i,j-1,k) )
     &         +maskFk(i, j )*
     &            ( Tracer(i, j ,km1)-Tracer(i, j ,k) )
     &        )
     &    +0.5 _d 0*recip_drC(kp1)*
     &        ( maskC(i,j-1,k,bi,bj)*maskC(i,j-1,kp1,bi,bj)*maskp1*
     &            ( Tracer(i,j-1,k)-Tracer(i,j-1,kp1) )
     &         +maskC(i, j ,k,bi,bj)*maskC(i, j ,kp1,bi,bj)*maskp1*
     &            ( Tracer(i, j ,k)-Tracer(i, j ,kp1) )
     &        )          )

C-    Calculate the flux
        df(i,j) = df(i,j)
     &           -yA(i,j)*dCdy_S
     &           *slcon_Kvy(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)
     &           -yA(i,j)*dCdz_S
     &           *slcon_Kvz(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)

#ifdef SLANTWISE_TFLUX_DIAGNOSTICS
        dTdy_S = (theta(i,j,k,bi,bj)-theta(i,j-1,k,bi,bj))
     &          *_recip_dyC(i,j,bi,bj)*recip_deepFacC(k)
        dTdz_S = 0.5 _d 0*(
     &    +0.5 _d 0*recip_drC(k)*
     &        ( maskFk(i,j-1)*
     &            ( theta(i,j-1,km1,bi,bj)-theta(i,j-1,k,bi,bj) )
     &         +maskFk(i, j )*
     &            ( theta(i, j ,km1,bi,bj)-theta(i, j ,k,bi,bj) )
     &        )
     &    +0.5 _d 0*recip_drC(kp1)*
     &        ( maskC(i,j-1,k,bi,bj)*maskC(i,j-1,kp1,bi,bj)*maskp1*
     &            ( theta(i,j-1,k,bi,bj)-theta(i,j-1,kp1,bi,bj) )
     &         +maskC(i, j ,k,bi,bj)*maskC(i, j ,kp1,bi,bj)*maskp1*
     &            ( theta(i, j ,k,bi,bj)-theta(i, j ,kp1,bi,bj) )
     &        )          )
        slcon_Ty(i,j,k,bi,bj) =  -yA(i,j)*dTdy_S
     &           *slcon_Kvy(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)
     &           -yA(i,j)*dTdz_S
     &           *slcon_Kvz(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)
#endif

#ifdef SLANTWISE_SFLUX_DIAGNOSTICS
        dSdy_S = (salt(i,j,k,bi,bj)-salt(i,j-1,k,bi,bj))
     &          *_recip_dyC(i,j,bi,bj)*recip_deepFacC(k)
        dSdz_S = 0.5 _d 0*(
     &    +0.5 _d 0*recip_drC(k)*
     &        ( maskFk(i,j-1)*
     &            ( salt(i,j-1,km1,bi,bj)-salt(i,j-1,k,bi,bj) )
     &         +maskFk(i, j )*
     &            ( salt(i, j ,km1,bi,bj)-salt(i, j ,k,bi,bj) )
     &        )
     &    +0.5 _d 0*recip_drC(kp1)*
     &        ( maskC(i,j-1,k,bi,bj)*maskC(i,j-1,kp1,bi,bj)*maskp1*
     &            ( salt(i,j-1,k,bi,bj)-salt(i,j-1,kp1,bi,bj) )
     &         +maskC(i, j ,k,bi,bj)*maskC(i, j ,kp1,bi,bj)*maskp1*
     &            ( salt(i, j ,k,bi,bj)-salt(i, j ,kp1,bi,bj) )
     &        )          )
        slcon_Sy(i,j,k,bi,bj) =  -yA(i,j)*dSdy_S
     &           *slcon_Kvy(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)
     &           -yA(i,j)*dSdz_S
     &           *slcon_Kvz(i,j,k,bi,bj)
     &           *_maskS(i,j,k,bi,bj)
#endif

       ENDDO
      ENDDO

#endif /* ALLOW_SLANTWISE */

      RETURN
      END



