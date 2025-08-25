#include "CPP_OPTIONS.h"
#ifdef ALLOW_SLANTWISE

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

C     !LOCAL VARIABLES:
      INTEGER i,j,kp1,km1 ! local indices
      _RL maskp1 ! mask
      _RL Kyy_S, Kyz_S ! Tensor elements at the southern boundary
      _RL epsS ! Avoiding zero denominator
      _RL dR_W_jk, dR_W_jm1k, dR_W_jkm1, dR_W_jm1km1, dR_W_k, dR_W_km1 ! Gradients at W grids, used for calculating dCdz_S
      _RL dCdz_S, dCdy_S ! Tracer gradients at the southern boundary / W grids

      epsS = 1.D-20
      km1 = MAX(k-1,1)
      kp1 = MIN(k+1,Nr)

C--   Loop for each tile
      DO j = jMin, jMax
       DO i = iMin, iMax
C-    Calculate the tensor element at the southern boundary
        Kyy_S = ( _maskS(i,j  ,k,bi,bj)*slcon_Kyy(i,j  ,k,bi,bj)
     &          + _maskS(i,j-1,k,bi,bj)*slcon_Kyy(i,j-1,k,bi,bj) ) /
     &          (_maskS(i,j  ,k,bi,bj) + _maskS(i,j-1,k,bi,bj) + epsS)
        Kyz_S = ( _maskS(i,j  ,k,bi,bj)*slcon_Kyz(i,j  ,k,bi,bj)
     &          + _maskS(i,j-1,k,bi,bj)*slcon_Kyz(i,j-1,k,bi,bj) ) /
     &          (_maskS(i,j  ,k,bi,bj) + _maskS(i,j-1,k,bi,bj) + epsS)
C-    Calculate the tracer gradients at the southern boundary
C     Calculate dC/dy at the southern boundary
        dCdy_S = ( Tracer(i,j,k) - Tracer(i,j-1,k) ) *_maskS(i,j,k,bi,bj)
     &           * _recip_dyC(i,j,bi,bj)*recip_deepFacC(k)
C     Calculate dC/dz at W grids
        maskp1 = 1. _d 0
        IF (k.GE.Nr) maskp1 = 0. _d 0
        dR_W_jkm1   = _recip_drC(k) * ( Tracer(i,j  ,km1) - Tracer(i,j  ,k) )
        dR_W_jm1km1 = _recip_drC(k) * ( Tracer(i,j-1,km1) - Tracer(i,j-1,k) )
        dR_W_km1    = 0.5d0 * ( maskFk(i,j-1) * dR_W_jm1km1 + maskFk(i,j) * dR_W_jkm1 )
        dR_W_jk     = _recip_drC(kp1) * ( Tracer(i, j ,k) - Tracer(i, j ,kp1) )
        dR_W_jm1k   = _recip_drC(kp1) * ( Tracer(i,j-1,k) - Tracer(i,j-1,kp1) )
        dR_W_k      = 0.5d0*(_maskC(i,j-1,k ,bi,bj)*_maskC(i,j-1,kp1,bi,bj)*maskp1*dR_W_jm1k
     &                +      _maskC(i,j  ,k ,bi,bj)*_maskC(i,j  ,kp1,bi,bj)*maskp1*dR_W_jk)

        dCdz_S = 0.5d0 * ( dR_W_km1 + dR_W_k )
C-    Calculate the flux
        df(i,j)=df(i,j)-yA(i,j)*(Kyy_S*dCdy_S+Kyz_S*dCdz_S)*_maskS(i,j,k,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END

#endif /* ALLOW_SLANTWISE */








