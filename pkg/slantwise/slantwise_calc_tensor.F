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
      INTEGER i, j, kp1
      _RL d2r, clat, slat
      _RL Kn, Kr, inv_f2
      _RL maskp1, maskloc, epsmask

      d2r = 3.141592653589793D0 / 180.D0
      epsmask = 1.D-20

C     Calculating density gradients at the tracer grids

      DO K = 1, Nr
       DO i = 1-OLx, sNx+OLx
        DO j = 1-OLy, sNy+OLy-1

C     Filling masks for debug
#ifdef SLANTWISE_DEBUG
         slcon_maskS(i,j,K,bi,bj) = _maskS(i,j,K,bi,bj)
         slcon_maskW(i,j,K,bi,bj) = _maskW(i,j,K,bi,bj)
         slcon_maskC(i,j,K,bi,bj) = maskC(i,j,K,bi,bj)
#endif
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
        DO k = 2, Nr
         kp1 = MIN(k+1,Nr)
         maskloc = maskC(i,j,k,bi,bj)
         maskp1  = maskC(i,j,kp1,bi,bj)
         sigmaR_C(i,j,K) = (maskloc*sigmaR(i,j,k)
     &                    + maskp1*sigmaR(i,j,kp1))
     &                    / (maskloc+maskp1+epsmask)
        END DO
C     SigmaR is zero at the surface, but maskC is not
        sigmaR_C(i,j,1) = sigmaR(i,j,2)*maskC(i,j,2,bi,bj)
       END DO
      END DO

#ifdef SLANTWISE_DEBUG
      DO K = 1, Nr
       DO i = 1-OLx, sNx+OLx
        DO j = 1-OLy, sNy+OLy-1
         slcon_sigmay(i,j,K,bi,bj) = sigmaY_C(i,j,k)
         slcon_sigmaz(i,j,K,bi,bj) = sigmaR_C(i,j,k)
         slcon_sy(i,j,K,bi,bj) = sigmaY(i,j,k)
         slcon_sz(i,j,K,bi,bj) = sigmaR(i,j,k)
        END DO
       END DO
      END DO
#endif



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
           inv_f2 = 1/( fCori(i,j,bi,bj)*fCori(i,j,bi,bj)
     &             + fCoriCos(i,j,bi,bj)*fCoriCos(i,j,bi,bj) )
           slcon_bn(i,j,K,bi,bj) = 
     &           (clat*sigmaY_C(i,j,K)+slat*sigmaR_C(i,j,K))
     &          * gravitySign*gravity/rhoConst
C     If bn>0, it is stable
           IF (slcon_bn(i,j,K,bi,bj) .GT. 0.D0) THEN
            slcon_Kyy(i,j,K,bi,bj) = 0.D0
            slcon_Kyz(i,j,K,bi,bj) = 0.D0
            slcon_Kzz(i,j,K,bi,bj) = 0.D0
           ELSE
            Kn = slcon_cn * ABS(slcon_bn(i,j,K,bi,bj)) 
     &          * SQRT(ABS(slcon_bn(i,j,K,bi,bj))) * inv_f2 
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

