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
C     sigmaRs   :: Vertical density gradient at southern boundary
C     sigmaYr   :: Meridional density gradient at upper boundary


      INTEGER bi, bj, iMin, iMax, jMin, jMax, k
      _RL sigmaR(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaY(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaRs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL sigmaYr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr)
      _RL myTime
      INTEGER myIter, myThid

CEOP

#ifdef ALLOW_SLANTWISE

C     !LOCAL VARIABLES
      INTEGER i, j, kp1, maskp1
      _RL d2r, clat, slat, lats
      _RL Kn, Kr, inv_f2
      _RL maskupper, bn, localf, localfp

      d2r = 3.141592653589793D0 / 180.D0

      IF ( usingSphericalPolarGrid ) THEN
       inv_f2 = rotationPeriod * rotationPeriod / 16. _d 0 
     &            / 3.141592653589793D0 / 3.141592653589793D0
      ELSE
       inv_f2 = 1 / (f0 * f0 + fPrime * fPrime)
      END IF

C     Calculating Kvy and Kvz at southern boundary

      DO k=Nr,1,-1
       kp1 = MIN(Nr,k+1)
       maskp1 = 1. _d 0
       IF (k.GE.Nr) maskp1 = 0. _d 0

       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
         IF ( _maskS(i,j,k,bi,bj) .LE. 0. _d 0 ) THEN
          slcon_Kvy(i,j,k,bi,bj) = 0. _d 0
          slcon_Kvz(i,j,k,bi,bj) = 0. _d 0
         ELSE
          sigmaRs(i,j,k) = 0.25 _d 0 *(sigmaR(i,j-1,k)+sigmaR(i,j,k)
     &                     +(sigmaR(i,j-1,kp1)+sigmaR(i,j,kp1))*maskp1
     &                      )*_maskS(i,j,k,bi,bj)
          IF ( usingSphericalPolarGrid ) THEN
           lats = yG(i,j,bi,bj)
          ELSE
           lats = slcon_constlat
          END IF
          clat = COS(lats * d2r)
          slat = SIN(lats * d2r)
          bn = (clat*sigmaY(i,j,k)+slat*sigmaRs(i,j,k))
     &          * gravitySign*gravity/rhoConst
          IF (bn*slat .GE. 0.D0) THEN
           slcon_Kvy(i,j,k,bi,bj) = 0. _d 0
           slcon_Kvz(i,j,k,bi,bj) = 0. _d 0
          ELSE
           Kn = MIN(slcon_cn * ABS(bn)
     &         * SQRT(ABS(bn)) * inv_f2
     &         * (slcon_Hnv(i,j,k,bi,bj)*slcon_Hnv(i,j,k,bi,bj)),
     &          slcon_maxkappa)
           Kr = slcon_cr * Kn / slcon_cn
           slcon_Kvy(i,j,k,bi,bj) = (Kn*clat*clat + Kr*slat*slat)
           slcon_Kvz(i,j,k,bi,bj) = (Kn-Kr)*slat*clat
          ENDIF
         ENDIF
        END DO
       END DO
      END DO

C     Calculating kwy and kwz at upper boundary

      DO k=Nr,2,-1
       DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
         maskUpper = maskC(i,j,k-1,bi,bj)*maskC(i,j,k,bi,bj)
         IF ( maskUpper .LE. 0. _d 0 ) THEN
          slcon_Kwy(i,j,k,bi,bj) = 0. _d 0
          slcon_Kwz(i,j,k,bi,bj) = 0. _d 0
         ELSE
          sigmaYr(i,j,k) = 0.25 _d 0 * (
     &                        sigmaY(i,j+1,k-1)+sigmaY(i,j,k-1)
     &                       +sigmaY(i,j+1, k )+sigmaY(i,j, k )
     &                      )*maskUpper
          IF ( usingSphericalPolarGrid ) THEN
           lats = yC(i,j,bi,bj)
          ELSE
           lats = slcon_constlat
          END IF
          clat = COS(lats * d2r)
          slat = SIN(lats * d2r)
          bn = (clat*sigmaYr(i,j,k)+slat*sigmaR(i,j,k))
     &          * gravitySign*gravity/rhoConst
          IF (bn*slat .GE. 0.D0) THEN
           slcon_Kwy(i,j,k,bi,bj) = 0. _d 0
           slcon_Kwz(i,j,k,bi,bj) = 0. _d 0
          ELSE
           Kn = MIN(slcon_cn * ABS(bn)
     &         * SQRT(ABS(bn)) * inv_f2
     &         * (slcon_Hnw(i,j,k,bi,bj)*slcon_Hnw(i,j,k,bi,bj)),
     &          slcon_maxkappa)
           Kr = slcon_cr * Kn / slcon_cn
           slcon_Kwy(i,j,k,bi,bj) = (Kn-Kr)*slat*clat
           slcon_Kwz(i,j,k,bi,bj) = Kn*slat*slat + Kr*clat*clat
          ENDIF
         ENDIF
        END DO
       END DO
      END DO



#endif /* ALLOW_SLANTWISE */

      RETURN
      END

