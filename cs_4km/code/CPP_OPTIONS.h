#ifndef CPP_OPTIONS_H
#define CPP_OPTIONS_H

C CPP flags controlling particular source code features

C********* RELEVANT CHANGES *********

C o Allow the use of Non-Linear Free-Surface formulation
C   this implies that surface thickness (hFactors) vary with time
#define NONLIN_FRSURF

C o NEW OPTION to disable rStar (z*) code
#define DISABLE_RSTAR_CODE

C********* RELEVANT CHANGES *********

C o Shortwave heating as extra term in external_forcing.F
C Note: this should be a run-time option
#define SHORTWAVE_HEATING

C o Include/exclude phi_hyd calculation code
#define INCLUDE_PHIHYD_CALCULATION_CODE

C o Include/exclude call to S/R CONVECT
#define INCLUDE_CONVECT_CALL

C o Include/exclude call to S/R CALC_DIFFUSIVITY
#define INCLUDE_CALC_DIFFUSIVITY_CALL

C o Include/exclude Implicit vertical advection code
#undef INCLUDE_IMPLVERTADV_CODE

C o Include/exclude nonHydrostatic code
#undef ALLOW_NONHYDROSTATIC

C o Include pressure loading code
#define ATMOSPHERIC_LOADING

C o Allow latitudinally varying BryanLewis79 vertical diffusivity
#undef ALLOW_BL79_LAT_VARY

C o Allow full 3D specification of vertical viscosity
#define ALLOW_3D_VISCA4

C o Allow full 3D specification of vertical diffusivity
#define ALLOW_3D_DIFFKR

C o Use "Exact Convervation" of fluid in Free-Surface formulation
C   so that d/dt(eta) is exactly equal to - Div.Transport
#define EXACT_CONSERV

C o Use "OLD" UV discretisation near boundaries (*not* recommended)
C   Note - only works with  #undef NO_SLIP_LATERAL  in calc_mom_rhs.F
C          because the old code did not have no-slip BCs
#undef  OLD_ADV_BCS

C o Use LONG.bin, LATG.bin, etc., initialization for ini_curviliear_grid.F
C   Default is to use "new" grid files (OLD_GRID_IO undef) but OLD_GRID_IO
C   is still useful with, e.g., single-domain curvilinear configurations.
#define OLD_GRID_IO

C o Use thsice+seaice (old) call sequence: ice-Dyn,ice-Advect,ice-Thermo(thsice)
C              as opposed to new sequence: ice-Thermo(thsice),ice-Dyn,ice-Advect
#undef OLD_THSICE_CALL_SEQUENCE

C o Execution environment support options
#include "CPP_EEOPTIONS.h"

C o Include/exclude code specific to the ECCO/SEALION version.
C   AUTODIFF or EXF package.
C   Currently controled by a single header file
C   For this to work, PACKAGES_CONFIG.h needs to be included!
cph#if (defined (ALLOW_AUTODIFF) || \
cph     defined (ALLOW_ECCO) || \
cph     defined (ALLOW_EXF))
cph# include "ECCO_CPPOPTIONS.h"
cph#endif

#endif /* CPP_OPTIONS_H */

