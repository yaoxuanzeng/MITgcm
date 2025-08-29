#include "SLANTWISE_OPTIONS.h"

CBOP
C     !ROUTINE: SLANTWISE_CALC_TENSOR
C     !INTERFACE:
      SUBROUTINE SLANTWISE_CALC_TENSOR(
     I      bi, bj, iMin, iMax, jMin, jMax,
     I      sigmaR, sigmaY,
     I      myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *===================================================================*
C     | Build (y,z) diffusivity tensor directly for slantwise convection.
C     | Scaling:
C     |   Kn = 0.660 * |bn|^(3/2) * (2*Omega)^(-2) * Hn^2
C     |   Kr = 0.326 * Kn
C     |
C     | Tensor entries:
C     |   Kyy = Kn cos^2θ + Kr sin^2θ
C     |   Kyz = (Kn - Kr) sinθ cosθ
C     |   Kzz = Kn sin^2θ + Kr cos^2θ
C     *===================================================================*
C     \ev

C     !USES:
      IMPLICIT NONE

C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SLANTWISE.h"

C     !INPUT/OUTPUT PARAMETERS
C     bi, bj    :: Tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number in simulation
C     myThid    :: My Thread Id. number
C     sigmaR_C   :: Vertical density gradient at tracer grids
C     sigmaY_C   :: Meridional density gradient at tracer grids

      INTEGER bi, bj, iMin, iMax, jMin, jMax, K
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY_C(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaR_C(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL myTime
      INTEGER myIter, myThid
CEOP

#ifdef ALLOW_SLANTWISE

C     !LOCAL VARIABLES
      INTEGER i, j
      _RL d2r, clat, slat
      _RL bn, Kn, Kr, inv_f2
      _RL maskp1, maskloc, epsmask

      d2r = 3.141592653589793D0 / 180.D0
      epsmask = 1.D-20
      inv_f2 = 1.D0 / (2.D0*omega) / (2.D0*omega) 

C     Calculating density gradients at the tracer grids

      DO K = 1, Nr
       DO i = 1-OLx, sNx+OLx
        DO j = 1-OLy, sNy+OLy-1
         maskloc = _maskS(i,j,K,bi,bj)
         maskp1 = _maskS(i,j+1,K,bi,bj)
         sigmaY_C(i,j,K) = (maskloc*sigmaY(i,j,K)
     &                    + maskp1*sigmaY(i,j+1,K))
     &                    / (maskloc+maskp1+epsmask)
        END DO
        sigmaY_C(i,sNy+OLy,K) = sigmaY_C(i,sNy+OLy-1,K)
       END DO
      END DO
      
      DO j = 1-OLy, sNy+OLy
       DO i = 1-OLx, sNx+OLx
        sigmaR_C(i,j,1) = sigmaR(i,j,1)*_maskW(i,j,1,bi,bj)
       END DO
      END DO

      DO K = 2, Nr
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
         maskloc = _maskW(i,j,K,bi,bj)
         maskp1  = _maskW(i,j,K-1,bi,bj)
         sigmaR_C(i,j,K) = (maskloc*sigmaR(i,j,K)
     &                    + maskp1*sigmaR(i,j,K-1))
     &                    / (maskloc+maskp1+epsmask)
        END DO
       END DO
      END DO

C     Calculating tensor elements

      DO K = 1, Nr
       DO j = 1-OLy, sNy+OLy
        DO i = 1-OLx, sNx+OLx
          IF ( maskC(i,j,K,bi,bj) .LE. 0.D0 ) THEN
           slcon_Kyy(i,j,K,bi,bj) = 0.D0
           slcon_Kyz(i,j,K,bi,bj) = 0.D0
           slcon_Kzz(i,j,K,bi,bj) = 0.D0
          ELSE
           clat = COS( slcon_lat(i,j,K,bi,bj) * d2r )
           slat = SIN( slcon_lat(i,j,K,bi,bj) * d2r )
           bn = (clat*sigmaY_C(i,j,K)+slat*sigmaR_C(i,j,K))
     &          * gravitySign*gravity/rhoConst
C     If bn>0, it is stable
           IF (bn .GT. 0.D0) THEN
            slcon_Kyy(i,j,K,bi,bj) = 0.D0
            slcon_Kyz(i,j,K,bi,bj) = 0.D0
            slcon_Kzz(i,j,K,bi,bj) = 0.D0
           ELSE
            Kn = slcon_cn * ABS(bn) * SQRT(ABS(bn)) * inv_f2 
     &          * (slcon_Hn(i,j,K,bi,bj)*slcon_Hn(i,j,K,bi,bj))
            Kr = slcon_cr * Kn
C           Rotate diag(Kn,Kr) to (y,z) and store
            slcon_Kyy(i,j,K,bi,bj) = Kn*clat*clat + Kr*slat*slat
            slcon_Kyz(i,j,K,bi,bj) = (Kn-Kr)*slat*clat
            slcon_Kzz(i,j,K,bi,bj) = Kn*slat*slat + Kr*clat*clat
           ENDIF
          slcon_Kyy(i,j,K,bi,bj) = slcon_Kyy(i,j,K,bi,bj) 
     &                            * maskC(i,j,K,bi,bj)
          slcon_Kyz(i,j,K,bi,bj) = slcon_Kyz(i,j,K,bi,bj) 
     &                            * maskC(i,j,K,bi,bj)
          slcon_Kzz(i,j,K,bi,bj) = slcon_Kzz(i,j,K,bi,bj) 
     &                            * maskC(i,j,K,bi,bj)
         ENDIF
        END DO
       END DO
      END DO

#endif /* ALLOW_SLANTWISE */

      RETURN
      END

