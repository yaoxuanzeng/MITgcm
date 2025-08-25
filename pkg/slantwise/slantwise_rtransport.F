#include "CPP_OPTIONS.h"
#ifdef ALLOW_SLANTWISE

CBOP
C     !ROUTINE: SLANTWISE_RTRANSPORT
C     !INTERFACE:
      SUBROUTINE SLANTWISE_RTRANSPORT(
     I     trIdentity, bi, bj, k,
     I     iMin, iMax, jMin, jMax,
     I     rA, maskFk, Tracer,
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

C     !INPUT/OUTPUT PARAMETERS:
C     trIdentity   :: tracer Id number
C     bi, bj       :: current tile indices
C     k            :: current level index
C     iMin,iMax    :: Range of 1rst index where results will be set
C     jMin,jMax    :: Range of 2nd  index where results will be set
C     rA           :: Area of W face
C     maskUp       :: 2-D array for mask at W points
C     Tracer       :: 3D Tracer field
C     df           :: Diffusive flux component work array.
      
      INTEGER trIdentity, bi, bj, k
      INTEGER iMin, iMax, jMin, jMax
      _RS rA    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RS maskFk(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL Tracer(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL df    (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER myThid
CEOP

C     !LOCAL VARIABLES
      INTEGER i,j,km1 ! Local indices
      _RL epsW ! Avoiding zero denominator
      _RL dCdz_W, dCdy_W ! Tracer gradients at W grids
      _RL Kyz_W, Kzz_W ! Tensor elements at W grids
      _RL dY_S_jk, dY_S_jp1k, dY_S_jkm1, dY_S_jp1km1, dY_k, dY_km1 ! Gradients at the southern boundary, used to calculate dCdy_W 

      epsW  = 1.D-20
      km1  = MAX(k-1,1)

C--   Loop for each tile
      DO j = jMin, jMax
       DO i = iMin, iMax

C-    Calculate the tensor element at the W grids
        Kyz_W = (_maskW(i,j,k  ,bi,bj)*slcon_Kyz(i,j,k  ,bi,bj)
     &          +_maskW(i,j,km1,bi,bj)*slcon_Kyz(i,j,km1,bi,bj) ) 
     &          /(_maskW(i,j,k  ,bi,bj)+_maskW(i,j,km1,bi,bj)+epsW)
        Kzz_W = (_maskW(i,j,k  ,bi,bj)*slcon_Kzz(i,j,k  ,bi,bj)
     &          +_maskW(i,j,km1,bi,bj)*slcon_Kzz(i,j,km1,bi,bj) ) 
     &          /(_maskW(i,j,k  ,bi,bj)+_maskW(i,j,km1,bi,bj)+epsW)

C-    Calculate the tracer gradients at the W grids
C     Calculate dCdz_W at the W grids
        dCdz_W =  _recip_drF(k) * ( Tracer(i,j,k) - Tracer(i,j,km1) ) * _maskW(i,j,k,bi,bj)
C     Calculate dCdy_W at the W grids
        dY_S_jp1k   = (Tracer(i,j+1,k  )-Tracer(i,j  ,k  )) 
     &                *_recip_dyC(i,j+1,bi,bj)*recip_deepFacC(k  )*_maskS(i,j+1,k  ,bi,bj)
        dY_S_jk     = (Tracer(i,j  ,k  )-Tracer(i,j-1,k  ))     
     &                *_recip_dyC(i,j  ,bi,bj)*recip_deepFacC(k  )*_maskS(i,j  ,k  ,bi,bj)
        dY_S_jp1km1 = (Tracer(i,j+1,k-1)-Tracer(i,j  ,k-1))     
     &                *_recip_dyC(i,j+1,bi,bj)*recip_deepFacC(k-1)*_maskS(i,j+1,k-1,bi,bj)
        dY_S_jkm1   = (Tracer(i,j  ,k-1)-Tracer(i,j-1,k-1))
     &                *_recip_dyC(i,j  ,bi,bj)*recip_deepFacC(k-1)*_maskS(i,j  ,k-1,bi,bj)
        dY_k   = 0.5D0 * ( dY_S_jp1k   + dY_S_jk   )
        dY_km1 = 0.5D0 * ( dY_S_jp1km1 + dY_S_jkm1 )
        dCdy_W = 0.5D0 * ( _maskW(i,j,k  ,bi,bj) * dY_k
     &                 +   _maskW(i,j,k-1,bi,bj) * dY_km1 )

C-    Calculate the flux
        df(i,j) = df(i,j)-rA(i,j)*(Kyz_W*dCdy_W+Kzz_W*dCdz_W)*_maskW(i,j,k,bi,bj)*maskFk(i,j)

       END DO
      END DO

      RETURN
      END
#endif /* ALLOW_SLANTWISE */


