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
      _RL epsW ! Avoiding zero denominator
      _RL op5 ! 0.5
      _RL dCdz_W, dCdy_W ! Tracer gradients at W grids
      _RL Kyz_W, Kzz_W ! Tensor elements at W grids

      epsW  = 1.D-20
      op5 = 0.5 _d 0
      km1  = MAX(k-1,1)

C--   Loop for each tile
      DO j = jMin, jMax
       DO i = iMin, iMax

C-    Calculate the tensor element at the W grids
        Kyz_W = (maskC(i,j,k  ,bi,bj)*slcon_Kyz(i,j,k  ,bi,bj)
     &          +maskC(i,j,km1,bi,bj)*slcon_Kyz(i,j,km1,bi,bj) ) 
     &          /(maskC(i,j,k  ,bi,bj)+maskC(i,j,km1,bi,bj)+epsW)
C     Adding Kzz in the slantwise_calc_diff.F
C        Kzz_W = (maskC(i,j,k  ,bi,bj)*slcon_Kzz(i,j,k  ,bi,bj)
C     &          +maskC(i,j,km1,bi,bj)*slcon_Kzz(i,j,km1,bi,bj) ) 
C     &          /(maskC(i,j,k  ,bi,bj)+maskC(i,j,km1,bi,bj)+epsW)

C-    Calculate the tracer gradients at the W grids
C     Calculate dCdz_W at the W grids 
C        dCdz_W =  _recip_drF(k) * ( Tracer(i,j,k) - Tracer(i,j,km1) ) * maskC(i,j,k,bi,bj)
C     Calculate dCdy_W at the W grids, copied from gmredi
        IF (k.GT.1) THEN
         dCdy_W = op5*(
     &     op5*recip_deepFacC(k)
     &        *( _maskS(i,j+1,k,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( Tracer(i,j+1,k) - Tracer(i,j,k) )
     &         + _maskS(i, j ,k,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( Tracer(i,j,k) - Tracer(i,j-1,k) )
     &         )
     &    +op5*recip_deepFacC(k-1)
     &        *( _maskS(i,j+1,k-1,bi,bj)*_recip_dyC(i,j+1,bi,bj)
     &              *( Tracer(i,j+1,k-1) - Tracer(i,j,k-1) )
     &         + _maskS(i, j ,k-1,bi,bj)*_recip_dyC(i, j ,bi,bj)
     &              *( Tracer(i,j,k-1) - Tracer(i,j-1,k-1) )
     &         )          )
        END IF

C-    Calculate the flux
         slcon_qzkyz(i,j,k,bi,bj) = -_rA(i,j,bi,bj)*(Kyz_W*dCdy_W)
     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
C         slcon_qzkzz(i,j,k,bi,bj) = -_rA(i,j,bi,bj)*(Kyz_W*dCdy_W)
C     &            *maskC(i,j,k,bi,bj)*maskFk(i,j)
         df(i,j) = df(i,j)
     &            + slcon_qzkyz(i,j,k,bi,bj)

       END DO
      END DO

#endif /* ALLOW_SLANTWISE */

      RETURN
      END


