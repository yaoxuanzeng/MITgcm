# Slantwise Convection Package for MITgcm

Implementation of the slantwise convection parameterization described in Zeng & Jansen 2026: "Parameterizing slantwise convection in icy moon oceans".

To enable the package:

1. Add `slantwise` to `packages.conf` before compiling.
2. Set `useSlantwise=.TRUE.` in `data.pkg`.

Example `data.slantwise`:

&SLANTWISE_PARM01
 slcon_constlat = 60.,       # Reference latitude (in degree) used to calculate Hz in Cartesian coordinates
 slcon_constH   = 40000.,    # Ocean depth used to calculate Hz
 slcon_maxkappa = 10000.,    # Maximum kappa_SL cap for numerical stability
&

