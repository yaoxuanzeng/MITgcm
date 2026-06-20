#ifndef SLANTWISE_OPTIONS_H
#define SLANTWISE_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: SLANTWISE_OPTIONS.h
C !INTERFACE:
C #include "SLANTWISE_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for slantwise package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_SLANTWISE
C     Package-specific compile-time options for SLANTWISE.

C--   Diagnostics (enable registration/fill of fields via DIAGNOSTICS package)
#define SLANTWISE_ALLOW_DIAGNOSTICS

C--   Allow time-varying H computed from stratification, instead of constant H
C     (runtime flag still controls use; this switch compiles the code path).
#undef  SLANTWISE_ALLOW_TIMEVAR_H

C--   Debug: extra outputs under #ifdef SLANTWISE_DEBUG
#define  SLANTWISE_DEBUG

C--   Allow temperature flux output
#define SLANTWISE_TFLUX_DIAGNOSTICS

C--   Allow salinity flux output
#define SLANTWISE_SFLUX_DIAGNOSTICS

#endif /* ALLOW_SLANTWISE */
#endif /* SLANTWISE_OPTIONS_H */

