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

C--   Scalar parameters
      COMMON /SLANTWISE_PARAMS/
     &   	  slcon_on, slcon_useconstH,
     &		  slcon_constlat,slcon_constH,
     &		  slcon_cn, slcon_cr, slcon_maxkappa

C--   Logical-type parameters
C     slcon_on			:: apply slantwise convection package
C     slcon_useconstH		:: use constant H (True) or stratification-defined H (False)

      LOGICAL slcon_on
      LOGICAL slcon_useconstH

C--   Other parameters
C     slcon_constlat		:: constant latitude
C     slcon_constH		:: constant height (H_n)
C     slcon_cn			:: pre-factor for axial diffusion
C     slcon_cr                  :: pre-factor for radial diffusion
C     slcon_maxkappa		:: maximum kappa for stability

      _RL slcon_constlat
      _RL slcon_constH
      _RL slcon_cn
      _RL slcon_cr
      _RL slcon_maxkappa

      CHARACTER*(MAX_LEN_FNAM) slconHnvFile
      CHARACTER*(MAX_LEN_FNAM) slconHnwFile

C--   Input fields
      COMMON /SLANTWISE_INPUT_FIELDS/
     &				slconHnvFile, slconHnwFile

C--   Runtime fields
      COMMON /SLANTWISE_FIELDS/
     &                          slcon_Hnv, slcon_Hnw,
     &				slcon_Kvy, slcon_Kvz,
     &				slcon_Kwy, slcon_Kwz,
     &				slcon_Ty, slcon_Tz,
     &				slcon_Sy, slcon_Sz

C     slcon_Hnv/w		:: Hn field at v/w points
C     slcon_Kvy/Kvz/Kwy/Kwz	:: Diffusivity tensor entries
C     slcon_Ty/Tz/Sy/Sz		:: Temperature/Salinity flux at y/z directions

      _RL slcon_Hnv(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Hnw(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kvy(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kvz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kwy(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Kwz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Ty(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Tz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Sy(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_Sz(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)

C--   Debug list
      COMMON /SLANTWISE_DEBUG/
     &                        slcon_debug1, slcon_debug2, slcon_debug3,
     &                        slcon_debug4, slcon_debug5, slcon_debug6,
     &                        slcon_debug7, slcon_debug8, slcon_debug9

      _RL slcon_debug1(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug2(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug3(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug4(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug5(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug6(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug7(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug8(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)
      _RL slcon_debug9(1-OLx:sNx+OLx, 1-OLy:sNy+OLy, 1:Nr, nSx, nSy)


#endif


