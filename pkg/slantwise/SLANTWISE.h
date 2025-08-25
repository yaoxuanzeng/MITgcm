#ifdef ALLOW_SLANTWISE
CBOP
C     !ROUTINE: SLANTWISE.h
C     !INTERFACE:
C     !#include "SLANTWISE.h"

C     !DESCRIPTION:
C     *================================================================*
C     | SLANTWISE.h
C     | o Header file defining slantwise convection scheme parameters and variables
C     | written by Yaoxuan Zeng at Aug, 2025
C     *================================================================*
CEOP

#include "SIZE.h"

C--   Scalar parameters
      COMMON /SLANTWISE_PARAMS/
     &   	  slcon_on, slcon_useconstH,
     &		  slcon_constlat,slcon_constH,
     &		  slcon_cn, slcon_cr

C--   Logical-type parameters
C     slcon_on		:: apply slantwise convection package
C     slcon_useconstH		:: use constant H (True) or stratification-defined H (False)

      LOGICAL slcon_on
      LOGICAL slcon_useconstH

C--   Other parameters
C     slcon_constlat		:: constant latitude
C     slcon_constH		:: constant height (H_n)
C     slcon_cn			:: pre-factor for axial diffusion
C     slcon_cr                  :: pre-factor for radial diffusion

      _RL slcon_constlat
      _RL slcon_constH
      _RL slcon_cn
      _RL slcon_cr

C--   Runtime fields
      COMMON /SLANTWISE_FIELDS/
     &                          slcon_lat, slcon_Hn, 
     &				slcon_Kyy, slcon_Kyz, slcon_Kzz

C     slcon_lat		:: local latitude
C     slcon_Hn		:: Hn field
C     slcon_Kyy/Kyz/Kzz :: Diffusivity tensor entries
     
      _RL slcon_lat(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Hn(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kyy(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kyz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kzz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)

#endif


