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

C--   Keep/intermediate tensor components (Kyy, Kyz, Kzz) for diagnostics
#define SLANTWISE_KEEP_TENSOR

C--   Cheap debug: extra checks/prints under #ifdef SLANTWISE_DEBUG
#undef  SLANTWISE_DEBUG

C--   Placeholder: allow coupling to cost function (not used by default)
#undef  ALLOW_SLANTWISE_COST

#endif /* ALLOW_SLANTWISE */
#endif /* SLANTWISE_OPTIONS_H */

