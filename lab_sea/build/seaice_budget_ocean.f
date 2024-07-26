

















CBOP
C !ROUTINE: CPP_OPTIONS.h
C !INTERFACE:
C #include "CPP_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | main CPP options file for the model:
C | Control which optional features to compile in model/src code.
C *==================================================================*
CEOP

C CPP flags controlling particular source code features

C-- Forcing code options:

C o Shortwave heating as extra term in external_forcing.F
C Note: this should be a run-time option

C o Include/exclude Geothermal Heat Flux at the bottom of the ocean

C o Allow to account for heating due to friction (and momentum dissipation)

C o Allow mass source or sink of Fluid in the interior
C   (3-D generalisation of oceanic real-fresh water flux)

C o Include pressure loading code

C o Include/exclude balancing surface forcing fluxes code

C o Include/exclude balancing surface forcing relaxation code

C o Include/exclude checking for negative salinity

C-- Options to discard parts of the main code:

C o Exclude/allow external forcing-fields load
C   this allows to read & do simple linear time interpolation of oceanic
C   forcing fields, if no specific pkg (e.g., EXF) is used to compute them.
C   If defined, use same method (with pkg/autodiff compiled or not) for checking
C   when to load new reccord ; by default, use simpler method with pkg/autodiff.

C o Include/exclude phi_hyd calculation code

C o Include/exclude sound speed calculation code
C o (Note that this is a diagnostic from Del Grasso algorithm, not derived
C    from EOS)

C-- Vertical mixing code options:

C o Include/exclude calling S/R CONVECTIVE_ADJUSTMENT

C o Include/exclude calling S/R CONVECTIVE_ADJUSTMENT_INI, turned off by
C   default because it is an unpopular historical left-over

C o Include/exclude call to S/R CALC_DIFFUSIVITY

C o Allow full 3D specification of vertical diffusivity

C o Allow latitudinally varying BryanLewis79 vertical diffusivity

C o Exclude/allow partial-cell effect (physical or enhanced) in vertical mixing
C   this allows to account for partial-cell in vertical viscosity and diffusion,
C   either from grid-spacing reduction effect or as artificially enhanced mixing
C   near surface & bottom for too thin grid-cell

C o Exclude/allow to use isotropic 3-D Smagorinsky viscosity as diffusivity
C   for tracers (after scaling by constant Prandtl number)

C-- Time-stepping code options:

C o Include/exclude combined Surf.Pressure and Drag Implicit solver code

C o Include/exclude Implicit vertical advection code

C o Include/exclude AdamsBashforth-3rd-Order code

C o Include/exclude Quasi-Hydrostatic Stagger Time-step AdamsBashforth code

C-- Model formulation options:

C o Allow/exclude "Exact Convervation" of fluid in Free-Surface formulation
C   that ensures that d/dt(eta) is exactly equal to - Div.Transport

C o Allow the use of Non-Linear Free-Surface formulation
C   this implies that grid-cell thickness (hFactors) varies with time
C o Disable code for rStar coordinate and/or code for Sigma coordinate
c#define DISABLE_RSTAR_CODE
c#define DISABLE_SIGMA_CODE

C o Include/exclude nonHydrostatic code

C o Include/exclude GM-like eddy stress in momentum code

C-- Algorithm options:

C o Include/exclude code for Non Self-Adjoint (NSA) conjugate-gradient solver

C o Include/exclude code for single reduction Conjugate-Gradient solver

C o Choices for implicit solver routines solve_*diagonal.F
C   The following has low memory footprint, but not suitable for AD
C   The following one suitable for AD but does not vectorize

C   Implementation alternative (might be faster on some platforms ?)

C-- Retired code options:

C o ALLOW isotropic scaling of harmonic and bi-harmonic terms when
C   using an locally isotropic spherical grid with (dlambda) x (dphi*cos(phi))
C *only for use on a lat-lon grid*
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
C   The preferred method is specifying a value for viscAhGrid or viscA4Grid
C   in data which is then automatically scaled by the grid size;
C   the old method of specifying viscAh/viscA4 and this flag is provided
C   for completeness only (and for use with the adjoint).
c#define ISOTROPIC_COS_SCALING

C o This flag selects the form of COSINE(lat) scaling of bi-harmonic term.
C *only for use on a lat-lon grid*
C   Has no effect if ISOTROPIC_COS_SCALING is undefined.
C   Has no effect on vector invariant momentum equations.
C   Setting this flag here affects both momentum and tracer equation unless
C   it is set/unset again in other header fields (e.g., GAD_OPTIONS.h).
C   The definition of the flag is commented to avoid interference with
C   such other header files.
c#define COSINEMETH_III

C o Use LONG.bin, LATG.bin, etc., initialization for ini_curviliear_grid.F
C   Default is to use "new" grid files (OLD_GRID_IO undef) but OLD_GRID_IO
C   is still useful with, e.g., single-domain curvilinear configurations.

C-- Other option files:

C o Execution environment support options

CBOP
C     !ROUTINE: CPP_EEOPTIONS.h
C     !INTERFACE:
C     include "CPP_EEOPTIONS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP\_EEOPTIONS.h                                         |
C     *==========================================================*
C     | C preprocessor "execution environment" supporting        |
C     | flags. Use this file to set flags controlling the        |
C     | execution environment in which a model runs - as opposed |
C     | to the dynamical problem the model solves.               |
C     | Note: Many options are implemented with both compile time|
C     |       and run-time switches. This allows options to be   |
C     |       removed altogether, made optional at run-time or   |
C     |       to be permanently enabled. This convention helps   |
C     |       with the data-dependence analysis performed by the |
C     |       adjoint model compiler. This data dependency       |
C     |       analysis can be upset by runtime switches that it  |
C     |       is unable to recoginise as being fixed for the     |
C     |       duration of an integration.                        |
C     |       A reasonable way to use these flags is to          |
C     |       set all options as selectable at runtime but then  |
C     |       once an experimental configuration has been        |
C     |       identified, rebuild the code with the appropriate  |
C     |       options set at compile time.                       |
C     *==========================================================*
CEOP

C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C=== Macro related options ===
C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working set size.
C     However, on vector CRAY systems this degrades performance.
C     Enable to switch REAL4_IS_SLOW from genmake2 (with LET_RS_BE_REAL4):

C--   Control use of "double" precision constants.
C     Use D0 where it means REAL*8 but not where it means REAL*16

C=== IO related options ===
C--   Flag used to indicate whether Fortran formatted write
C     and read are threadsafe. On SGI the routines can be thread
C     safe, on Sun it is not possible - if you are unsure then
C     undef this option.

C--   Flag used to indicate whether Binary write to Local file (i.e.,
C     a different file for each tile) and read are thread-safe.

C--   Flag to turn off the writing of error message to ioUnit zero

C--   Alternative formulation of BYTESWAP, faster than
C     compiler flag -byteswapio on the Altix.

C--   Flag to turn on old default of opening scratch files with the
C     STATUS='SCRATCH' option. This method, while perfectly FORTRAN-standard,
C     caused filename conflicts on some multi-node/multi-processor platforms
C     in the past and has been replace by something (hopefully) more robust.

C--   Flag defined for eeboot_minimal.F, eeset_parms.F and open_copy_data_file.F
C     to write STDOUT, STDERR and scratch files from process 0 only.
C WARNING: to use only when absolutely confident that the setup is working
C     since any message (error/warning/print) from any proc <> 0 will be lost.

C=== MPI, EXCH and GLOBAL_SUM related options ===
C--   Flag turns off MPI_SEND ready_to_receive polling in the
C     gather_* subroutines to speed up integrations.

C--   Control MPI based parallel processing
CXXX We no longer select the use of MPI via this file (CPP_EEOPTIONS.h)
CXXX To use MPI, use an appropriate genmake2 options file or use
CXXX genmake2 -mpi .
CXXX #undef  ALLOW_USE_MPI

C--   Control use of communication that might overlap computation.
C     Under MPI selects/deselects "non-blocking" sends and receives.
C--   Control use of communication that is atomic to computation.
C     Under MPI selects/deselects "blocking" sends and receives.

C--   Control XY periodicity in processor to grid mappings
C     Note: Model code does not need to know whether a domain is
C           periodic because it has overlap regions for every box.
C           Model assume that these values have been
C           filled in some way.

C--   disconnect tiles (no exchange between tiles, just fill-in edges
C     assuming locally periodic subdomain)

C--   Always cumulate tile local-sum in the same order by applying MPI allreduce
C     to array of tiles ; can get slower with large number of tiles (big set-up)

C--   Alternative way of doing global sum without MPI allreduce call
C     but instead, explicit MPI send & recv calls. Expected to be slower.

C--   Alternative way of doing global sum on a single CPU
C     to eliminate tiling-dependent roundoff errors. Note: This is slow.

C=== Other options (to add/remove pieces of code) ===
C--   Flag to turn on checking for errors from all threads and procs
C     (calling S/R STOP_IF_ERROR) before stopping.

C--   Control use of communication with other component:
C     allow to import and export from/to Coupler interface.

C--   Activate some pieces of code for coupling to GEOS AGCM

C=== And define Macros ===
CBOP
C     !ROUTINE: CPP_EEMACROS.h
C     !INTERFACE:
C     include "CPP_EEMACROS.h"
C     !DESCRIPTION:
C     *==========================================================*
C     | CPP_EEMACROS.h
C     *==========================================================*
C     | C preprocessor "execution environment" supporting
C     | macros. Use this file to define macros for  simplifying
C     | execution environment in which a model runs - as opposed
C     | to the dynamical problem the model solves.
C     *==========================================================*
CEOP


C     In general the following convention applies:
C     ALLOW  - indicates an feature will be included but it may
C     CAN      have a run-time flag to allow it to be switched
C              on and off.
C              If ALLOW or CAN directives are "undef'd" this generally
C              means that the feature will not be available i.e. it
C              will not be included in the compiled code and so no
C              run-time option to use the feature will be available.
C
C     ALWAYS - indicates the choice will be fixed at compile time
C              so no run-time option will be present

C     Flag used to indicate which flavour of multi-threading
C     compiler directives to use. Only set one of these.
C     USE_SOLARIS_THREADING  - Takes directives for SUN Workshop
C                              compiler.
C     USE_KAP_THREADING      - Takes directives for Kuck and
C                              Associates multi-threading compiler
C                              ( used on Digital platforms ).
C     USE_IRIX_THREADING     - Takes directives for SGI MIPS
C                              Pro Fortran compiler.
C     USE_EXEMPLAR_THREADING - Takes directives for HP SPP series
C                              compiler.
C     USE_C90_THREADING      - Takes directives for CRAY/SGI C90
C                              system F90 compiler.






C--   Define the mapping for the _BARRIER macro
C     On some systems low-level hardware support can be accessed through
C     compiler directives here.

C--   Define the mapping for the BEGIN_CRIT() and  END_CRIT() macros.
C     On some systems we simply execute this section only using the
C     master thread i.e. its not really a critical section. We can
C     do this because we do not use critical sections in any critical
C     sections of our code!

C--   Define the mapping for the BEGIN_MASTER_SECTION() and
C     END_MASTER_SECTION() macros. These are generally implemented by
C     simply choosing a particular thread to be "the master" and have
C     it alone execute the BEGIN_MASTER..., END_MASTER.. sections.

CcnhDebugStarts
C      Alternate form to the above macros that increments (decrements) a counter each
C      time a MASTER section is entered (exited). This counter can then be checked in barrier
C      to try and detect calls to BARRIER within single threaded sections.
C      Using these macros requires two changes to Makefile - these changes are written
C      below.
C      1 - add a filter to the CPP command to kill off commented _MASTER lines
C      2 - add a filter to the CPP output the converts the string N EWLINE to an actual newline.
C      The N EWLINE needs to be changes to have no space when this macro and Makefile changes
C      are used. Its in here with a space to stop it getting parsed by the CPP stage in these
C      comments.
C      #define IF ( a .EQ. 1 ) THEN  IF ( a .EQ. 1 ) THEN  N EWLINE      CALL BARRIER_MS(a)
C      #define ENDIF    CALL BARRIER_MU(a) N EWLINE        ENDIF
C      'CPP = cat $< | $(TOOLSDIR)/set64bitConst.sh |  grep -v '^[cC].*_MASTER' | cpp  -traditional -P'
C      .F.f:
C      $(CPP) $(DEFINES) $(INCLUDES) |  sed 's/N EWLINE/\n/' > $@
CcnhDebugEnds

C--   Control storage of floating point operands
C     On many systems it improves performance only to use
C     8-byte precision for time stepped variables.
C     Constant in time terms ( geometric factors etc.. )
C     can use 4-byte precision, reducing memory utilisation and
C     boosting performance because of a smaller working
C     set size. However, on vector CRAY systems this degrades
C     performance.
C- Note: global_sum/max macros were used to switch to  JAM routines (obsolete);
C  in addition, since only the R4 & R8 S/R are coded, GLOBAL RS & RL macros
C  enable to call the corresponding R4 or R8 S/R.



C- Note: a) exch macros were used to switch to  JAM routines (obsolete)
C        b) exch R4 & R8 macros are not practically used ; if needed,
C           will directly call the corrresponding S/R.

C--   Control use of JAM routines for Artic network (no longer supported)
C     These invoke optimized versions of "exchange" and "sum" that
C     utilize the programmable aspect of Artic cards.
CXXX No longer supported ; started to remove JAM routines.
CXXX #ifdef LETS_MAKE_JAM
CXXX #define CALL GLOBAL_SUM_R8 ( a, b) CALL GLOBAL_SUM_R8_JAM ( a, b)
CXXX #define CALL GLOBAL_SUM_R8 ( a, b ) CALL GLOBAL_SUM_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RS ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XY_RL ( a, b ) CALL EXCH_XY_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RS ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #define CALL EXCH_XYZ_RL ( a, b ) CALL EXCH_XYZ_R8_JAM ( a, b )
CXXX #endif

C--   Control use of "double" precision constants.
C     Use d0 where it means REAL*8 but not where it means REAL*16

C--   Substitue for 1.D variables
C     Sun compilers do not use 8-byte precision for literals
C     unless .Dnn is specified. CRAY vector machines use 16-byte
C     precision when they see .Dnn which runs very slowly!

C--   Set the format for writing processor IDs, e.g. in S/R eeset_parms
C     and S/R open_copy_data_file. The default of I9.9 should work for
C     a long time (until we will use 10e10 processors and more)

C--   Set the format for writing ensemble task IDs in S/R eeset_parms
C     and S/R open_copy_data_file.

C--   Set ACTION= in OPEN instruction for input file (before doing IO)
C     leave it empty (if EXCLUDE_OPEN_ACTION) or set it to proper value



C o Include/exclude single header file containing multiple packages options
C   (AUTODIFF, COST, CTRL, ECCO, EXF ...) instead of the standard way where
C   each of the above pkg get its own options from its specific option file.
C   Although this method, inherited from ECCO setup, has been traditionally
C   used for all adjoint built, work is in progress to allow to use the
C   standard method also for adjoint built.
c#ifdef 
c# include "ECCO_CPPOPTIONS.h"
c#endif


C     *==========================================================*
C     | SEAICE_OPTIONS.h
C     | o CPP options file for sea ice package.
C     *==========================================================*
C     | Use this file for selecting options within the sea ice
C     | package.
C     *==========================================================*

C---  Package-specific Options & Macros go here

C--   Write "text-plots" of certain fields in STDOUT for debugging.

C--   By default, the sea-ice package uses its own integrated bulk
C     formulae to compute fluxes (fu, fv, EmPmR, Qnet, and Qsw) over
C     open-ocean.  When this flag is set, these variables are computed
C     in a separate external package, for example, pkg/exf, and then
C     modified for sea-ice effects by pkg/seaice.

C--   This CPP flag has been retired.  The number of ice categories
C     used to solve for seaice flux is now specified by run-time
C     parameter SEAICE_multDim.
C     Note: be aware of pickup_seaice.* compatibility issues when
C     restarting a simulation with a different number of categories.
c#define SEAICE_MULTICATEGORY

C--   run with sea Ice Thickness Distribution (ITD);
C     set number of categories (nITD) in SEAICE_SIZE.h

C--   Since the missing sublimation term is now included
C     this flag is needed for backward compatibility

C--   Suspected missing term in coupled ocn-ice heat budget (to be confirmed)

C--   Default is constant seaice salinity (SEAICE_salt0); Define the following
C     flag to consider (space & time) variable salinity: advected and forming
C     seaice with a fraction (=SEAICE_saltFrac) of freezing seawater salinity.
C- Note: SItracer also offers an alternative way to handle variable salinity.

C--   Enable grease ice parameterization (requires to define ALLOW_SITRACER):
C     The grease ice parameterization delays formation of solid sea ice from
C     frazil ice by a time constant and provides a dynamic calculation of the
C     initial solid sea ice thickness HO as a function of winds, currents and
C     available grease ice volume. Grease ice does not significantly reduce heat
C     loss from the ocean in winter and area covered by grease is thus handled
C     like open water (For details see Smedsrud and Martin, 2014, Ann.Glac.).
C     Set SItrName(1) = 'grease' in namelist SEAICE_PARM03 in data.seaice
C     then output SItr01 is SItrNameLong(1) = 'grease ice volume fraction',
C     with SItrUnit(1) = '[0-1]', which needs to be multiplied by SIheff to
C     yield grease ice volume. Additionally, the actual grease ice layer
C     thickness (diagnostic SIgrsLT) can be saved.

C--   Tracers of ice and/or ice cover.
C-    To try avoid 'spontaneous generation' of tracer maxima by advdiff.

C-    Include code to diagnose sea ice tracer budgets in
C     seaice_advdiff.F and seaice_tracer_phys.F. Diagnostics are
C     computed the "call diagnostics_fill" statement is commented out.

C--   Allow sea-ice dynamic code. These options are provided so that,
C     if turned off (#undef), to compile (and process with TAF) only the
C     the thermodynamics component of the code. Note that, if needed,
C     sea-ice dynamics can be turned off at runtime (SEAICEuseDYNAMICS=F).

C--   Historically, the seaice model was discretized on a B-Grid. This
C     discretization should still work but it is not longer actively
C     tested and supported. Define this flag to compile it. It cannot be
C     defined together with SEAICE_CGRID

C--   The following flag should always be set in order to use C the
C--   operational C-grid discretization.

C--   Options for the C-grid version only:

C     enable advection of sea ice momentum

C     enable JFNK code by defining the following flag

C     enable Krylov code by defining the following flag

C--   Use a different order when mapping 2D velocity arrays to 1D vector
C     before passing it to FGMRES.

C     to reproduce old verification results for JFNK

C     enable LSR to use global (multi-tile) tri-diagonal solver

C     enable EVP code by defining the following flag
C-    When set use SEAICE_zetaMin and SEAICE_evpDampC to limit viscosities
C     from below and above in seaice_evp: not necessary, and not recommended

C     Include code to avoid underflows in EVP-code (copied from CICE).
C     Many compilers can handle this more efficiently with the help of a flag.

C     Include code to print residual of EVP iteration for debugging/diagnostics

C     smooth regularization (without max-function) of delta for
C     better differentiability

C     regularize zeta to zmax with a smooth tanh-function instead
C     of a min(zeta,zmax). This improves convergence of iterative
C     solvers (Lemieux and Tremblay 2009, JGR). No effect on EVP

C--   Different yield curves within the VP rheology framework
C     allow the truncated ellipse rheology (runtime flag SEAICEuseTEM)

C     allow the use of the Mohr Coulomb rheology (runtime flag SEAICEuseMCS)
C     as defined in (Ip 1991) /!\ This is known to give unstable results,
C     use with caution

C     allow the use of Mohr Coulomb with elliptical plastic potential
C     (runtime flag SEAICEuseMCE)

C     allow the teardrop and parabolic lens  rheology
C     (runtime flag SEAICEuseTD and SEAICEusePL)

C--   LSR solver settings
C     Use LSR vector code; not useful on non-vector machines, because it
C     slows down convergence considerably, but the extra iterations are
C     more than made up by the much faster code on vector machines. For
C     the only regularly test vector machine these flags a specified
C     in the build options file SUPER-UX_SX-8_sxf90_awi, so that we comment
C     them out here.

C     Use zebra-method (alternate lines) for line-successive-relaxation
C     This modification improves the convergence of the vector code
C     dramatically, so that is may actually be useful in general, but
C     that needs to be tested. Can be used without vectorization options.

C     Include code to print residual of nonlinear outer loop of LSR

C     This flag is also required for an actual adjoint of seaice_lsr;
C     increases memory requirements a lot.

C     Use parameterisation of grounding ice for a better representation
C     of fastice in shallow seas



C--   Some regularisations
C-    When set limit the Ice-Loading to mass of 1/5 of Surface ocean grid-box

C-    When set use SEAICE_clipVelocties = .true., to clip U/VICE at 40cm/s,
C     not recommended

C-    When set cap the sublimation latent heat flux in solve4temp according
C     to the available amount of ice+snow. Otherwise this term is treated
C     like all of the others -- residuals heat and fw stocks are passed to
C     the ocean at the end of seaice_growth in a conservative manner.
C     SEAICE_CAP_SUBLIM is not needed as of now, but kept just in case.

C--   AD flags
C-    TAF related flag, currently only used in seaice_ad_check_lev[1-4]_dir.h;
C     it is unclear if this is ever needed.

C-    Reset fields to zero to stabilise AD code of dynamics solver
C     (resulting in wrong gradients)

C-    Another flag to simplify dependencies for TAF-generated AD-code
C     the thermodynamic part, mostly by resetting variables to zero

C-    Special seaice flag for AD testing

C--   Use the adjointable sea-ice thermodynamic model
C     in seaice_growth_adx.F instead of seaice_growth.F
C     This options excludes more complex physics such
C     as sublimation, ITD, and frazil.

C--   These flags are not strictly AD-related but may help obtaining
C     simpler AD-code:
C-    Do not compile code that resets AREA (or AREAITD) to a mininum value
C     of SEAICE_area_floor (=SIeps with default of 1e-5) if there is
C     some finite sea ice thickness

C-    Do not compile growth/thermodynamics code (avoiding this code can
C     also be done by setting runtime parameter usePWthermodynamics=F)

C-    Do not compile/use seaice-related obcs code when using obcs.

C--   Enable free drift code

C--   pkg/seaice cost functions compile flags
C-    Sea-ice volume (requires pkg/cost)


CBOP
C !ROUTINE: EXF_OPTIONS.h
C !INTERFACE:
C #include "EXF_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for EXternal Forcing (EXF) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP


C-- Package-specific Options & Macros go here

C   --------------------
C   pkg/exf CPP options:
C   (see also table below on how to combine options)

C   > ( EXF_VERBOSE ) < replaced with run-time integer parameter "exf_debugLev"
C
C   >>> ALLOW_ATM_WIND <<<
C       If defined, set default value of run-time param. "useAtmWind" to True.
C       If useAtmWind=True, read-in and use wind vector (uwind/vwind)
C       to compute surface wind stress.
C
C   >>> ALLOW_ATM_TEMP <<<
C       This is the main EXF option controlling air-sea buoyancy fluxes:
C      If undefined, net heat flux (Qnet) and net fresh water flux
C       (EmP or EmPmR) are set according to hfluxfile & sfluxfile setting.
C      If defined, net heat flux and net fresh water flux are computed
C       from sum of various components (radiative SW,LW + turbulent heat
C       fluxes SH,LH ; Evap, Precip and optionally RunOff) thus ignoring
C       hfluxfile & sfluxfile.
C      In addition, it allows to read-in from files atmospheric temperature
C       and specific humidity, net radiative fluxes, and precip.
C       Also enable to read-in Evap (if EXF_READ_EVAP is defined) or
C       turbulent heat fluxes (if ALLOW_READ_TURBFLUXES is defined).
C
C   >>> ALLOW_DOWNWARD_RADIATION <<<
C       If defined, downward long-wave and short-wave radiation
C       can be read-in form files to compute net lwflux and swflux.
C
C   >>> ALLOW_ZENITHANGLE <<<
C       If defined, ocean albedo varies with the zenith angle, and
C       incoming fluxes at the top of the atmosphere are computed
C
C   >>> ALLOW_BULKFORMULAE <<<
C       Allows the use of bulk formulae in order to estimate
C       turbulent fluxes (Sensible,Latent,Evap) at the ocean surface.
C
C   >>> EXF_CALC_ATMRHO
C       Calculate the local air density as function of temp, humidity
C       and pressure
C
C   >>> EXF_READ_EVAP <<<
C       If defined, evaporation field is read-in from file;
C     Note: if ALLOW_BULKFORMULAE is defined, evap that is computed from
C       atmospheric state will be replaced by read-in evap but computed
C       latent heat flux will be kept.
C
C   >>> ALLOW_READ_TURBFLUXES <<<
C       If defined, turbulent heat fluxes (sensible and latent) can be read-in
C       from files (but overwritten if ALLOW_BULKFORMULAE is defined).
C
C   >>> ALLOW_RUNOFF <<<
C       If defined, river and glacier runoff can be read-in from files.
C
C   >>> ALLOW_SALTFLX <<<
C       If defined, upward salt flux can be read-in from files.
C
C   >>> ALLOW_RUNOFTEMP <<<
C       If defined, river and glacier runoff temperature
C       can be read-in from files.
C
C   >>>  <<<
C       If defined, atmospheric pressure can be read-in from files.
C   WARNING: this flag is set (define/undef) in CPP_OPTIONS.h
C            and cannot be changed here (in EXF_OPTIONS.h)
C
C   >>> EXF_ALLOW_TIDES <<<
C       If defined, 2-D tidal geopotential can be read-in from files
C
C   >>> EXF_SEAICE_FRACTION <<<
C       If defined, seaice fraction can be read-in from files (areaMaskFile)
C
C   >>> ALLOW_CLIMSST_RELAXATION <<<
C       Allow the relaxation of surface level temperature to SST (climatology),
C       e.g. the Reynolds climatology.
C
C   >>> ALLOW_CLIMSSS_RELAXATION <<<
C       Allow the relaxation of surface level salinity to SSS (climatology),
C       e.g. the Levitus climatology.
C
C   >>> USE_EXF_INTERPOLATION <<<
C       Allows to provide input field on arbitrary Lat-Lon input grid
C       (as specified in EXF_NML_04) and to interpolate to model grid.
C     Note: default is to interpolate unless {FLD}_interpMethod is set to 0
C
C   ====================================================================
C
C    The following CPP options:
C       ALLOW_ATM_WIND / useAtmWind (useWind)
C       ALLOW_ATM_TEMP               (TEMP)
C       ALLOW_DOWNWARD_RADIATION     (DOWN)
C       ALLOW_BULKFORMULAE           (BULK)
C       EXF_READ_EVAP                (EVAP)
C       ALLOW_READ_TURBFLUXES        (TURB)
C
C    permit all ocean-model forcing configurations listed in the 2 tables below.
C    The first configuration (A1,B1) is the flux-forced, ocean model.
C    Configurations A2,B3 and A2,B4 use pkg/exf open-water bulk formulae
C    to compute, from atmospheric variables, the missing surface fluxes.
C    The forcing fields in the rightmost column are defined in EXF_FIELDS.h
C    (ocean-model surface forcing field are defined in model/inc/FFIELDS.h)
C
C    (A) Surface momentum flux: [model: fu,fv ; exf: ustress,vstress]
C
C    # |useWind|        actions
C   ---|-------|-------------------------------------------------------------
C   (1)| False | Read-in ustress,vstress (if needed in B, compute wind-speed)
C      |       |
C   (2)| True  | Read-in uwind,vwind ; compute wind stress ustress,vstress.
C   ---|-------|-------------------------------------------------------------
C
C    (B) Surface buoyancy flux:
C        [ net heat flux: Qnet (exf: hflux), net short-wave: Qsw (exf: swflux)
C          fresh-water flux: EmPmR (exf: sflux) and saltFlux (exf: saltflx) ]
C
C    # |TEMP |DOWN |BULK |EVAP |TURB |            actions
C   ---|-----|-----|-----|-----|-----|-------------------------------------
C   (1)|  -  |  -  |  -  |  -  |  -  | Read-in hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (2)|  -  | def |  -  |  -  |  -  | Read-in hflux, swdown and sflux.
C      |     |     |     |     |     | Compute swflux.
C      |     |     |     |     |     |
C   (3)| def | def | def |  -  |  -  | Read-in atemp, aqh, swdown, lwdown,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (4)| def |  -  | def |  -  |  -  | Read-in atemp, aqh, swflux, lwflux,
C      |     |     |     |     |     |  precip, and runoff.
C      |     |     |     |     |     | Compute hflux and sflux.
C      |     |     |     |     |     |
C   (5)| def | def |  -  | def | def | Read-in hs, hl, swdown, lwdown,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute hflux, swflux and sflux.
C      |     |     |     |     |     |
C   (6)| def |  -  |  -  | def | def | Read-in hs, hl, swflux, lwflux,
C      |     |     |     |     |     |  evap, precip and runoff.
C      |     |     |     |     |     | Compute  hflux and sflux.
C
C   =======================================================================

C-  Bulk formulae related flags.
C Note: To use ALLOW_BULKFORMULAE or EXF_READ_EVAP, needs #define 
C use Large and Yeager (2004) modification to Large and Pond bulk formulae
C use drag formulation of Large and Yeager (2009), Climate Dyn., 33, pp 341-364

C-  Other forcing fields

C Note: To use EXF_CALC_ATMRHO, both 
C       and  need to be defined

C-  Zenith Angle/Albedo related flags.

C-  Use ocean_emissivity*lwdown in lwFlux. This flag should be defined
C   unless to reproduce old results (obtained with inconsistent old code)

C-  Surface level relaxation to prescribed fields (e.g., climatologies)

C-  Allows to read-in (2-d) tidal geopotential forcing

C-  Allows to read-in seaice fraction from files (areaMaskFile)

C-  Use spatial interpolation to interpolate
C   forcing files from input grid to model grid.
C   for interpolated vector fields, rotate towards model-grid axis
C   using old rotation formulae (instead of grid-angles)
C   for interpolation around N & S pole, use the old formulation
C   (no pole symmetry, single vector-comp interp, reset to 0 zonal-comp @ N.pole)


C-  Not recommended (not tested nor maintained) and un-documented Options:


CStartOfInterface
      SUBROUTINE SEAICE_BUDGET_OCEAN(
     I     UG,
     I     TSURF,
     O     netHeatFlux, SWHeatFlux,
     I     bi, bj, myTime, myIter, myThid )
C     *================================================================*
C     | SUBROUTINE seaice_budget_ocean
C     | o Calculate surface heat fluxes over open ocean
C     |   see Hibler, MWR, 108, 1943-1973, 1980
C     |   If  is defined this routine simply
C     |   copies the global fields to the seaice-local fields.
C     *================================================================*
      IMPLICIT NONE

C     === Global variables ===
CBOP
C    !ROUTINE: SIZE.h
C    !INTERFACE:
C    include SIZE.h
C    !DESCRIPTION: \bv
C     *==========================================================*
C     | SIZE.h Declare size of underlying computational grid.
C     *==========================================================*
C     | The design here supports a three-dimensional model grid
C     | with indices I,J and K. The three-dimensional domain
C     | is comprised of nPx*nSx blocks (or tiles) of size sNx
C     | along the first (left-most index) axis, nPy*nSy blocks
C     | of size sNy along the second axis and one block of size
C     | Nr along the vertical (third) axis.
C     | Blocks/tiles have overlap regions of size OLx and OLy
C     | along the dimensions that are subdivided.
C     *==========================================================*
C     \ev
C
C     Voodoo numbers controlling data layout:
C     sNx :: Number of X points in tile.
C     sNy :: Number of Y points in tile.
C     OLx :: Tile overlap extent in X.
C     OLy :: Tile overlap extent in Y.
C     nSx :: Number of tiles per process in X.
C     nSy :: Number of tiles per process in Y.
C     nPx :: Number of processes to use in X.
C     nPy :: Number of processes to use in Y.
C     Nx  :: Number of points in X for the full domain.
C     Ny  :: Number of points in Y for the full domain.
C     Nr  :: Number of points in vertical direction.
CEOP
      INTEGER sNx
      INTEGER sNy
      INTEGER OLx
      INTEGER OLy
      INTEGER nSx
      INTEGER nSy
      INTEGER nPx
      INTEGER nPy
      INTEGER Nx
      INTEGER Ny
      INTEGER Nr
      PARAMETER (
     &           sNx =  10,
     &           sNy =   8,
     &           OLx =   4,
     &           OLy =   4,
     &           nSx =   2,
     &           nSy =   2,
     &           nPx =   1,
     &           nPy =   1,
     &           Nx  = sNx*nSx*nPx,
     &           Ny  = sNy*nSy*nPy,
     &           Nr  =  23)

C     MAX_OLX :: Set to the maximum overlap region size of any array
C     MAX_OLY    that will be exchanged. Controls the sizing of exch
C                routine buffers.
      INTEGER MAX_OLX
      INTEGER MAX_OLY
      PARAMETER ( MAX_OLX = OLx,
     &            MAX_OLY = OLy )

CBOP
C     !ROUTINE: EEPARAMS.h
C     !INTERFACE:
C     include "EEPARAMS.h"
C
C     !DESCRIPTION:
C     *==========================================================*
C     | EEPARAMS.h                                               |
C     *==========================================================*
C     | Parameters for "execution environemnt". These are used   |
C     | by both the particular numerical model and the execution |
C     | environment support routines.                            |
C     *==========================================================*
CEOP

C     ========  EESIZE.h  ========================================

C     MAX_LEN_MBUF  :: Default message buffer max. size
C     MAX_LEN_FNAM  :: Default file name max. size
C     MAX_LEN_PREC  :: Default rec len for reading "parameter" files

      INTEGER MAX_LEN_MBUF
      PARAMETER ( MAX_LEN_MBUF = 512 )
      INTEGER MAX_LEN_FNAM
      PARAMETER ( MAX_LEN_FNAM = 512 )
      INTEGER MAX_LEN_PREC
      PARAMETER ( MAX_LEN_PREC = 200 )

C     MAX_NO_THREADS  :: Maximum number of threads allowed.
CC    MAX_NO_PROCS    :: Maximum number of processes allowed.
CC    MAX_NO_BARRIERS :: Maximum number of distinct thread "barriers"
      INTEGER MAX_NO_THREADS
      PARAMETER ( MAX_NO_THREADS =  4 )
c     INTEGER MAX_NO_PROCS
c     PARAMETER ( MAX_NO_PROCS   =  70000 )
c     INTEGER MAX_NO_BARRIERS
c     PARAMETER ( MAX_NO_BARRIERS = 1 )

C     Particularly weird and obscure voodoo numbers
C     lShare :: This wants to be the length in
C               [148]-byte words of the size of
C               the address "window" that is snooped
C               on an SMP bus. By separating elements in
C               the global sum buffer we can avoid generating
C               extraneous invalidate traffic between
C               processors. The length of this window is usually
C               a cache line i.e. small O(64 bytes).
C               The buffer arrays are usually short arrays
C               and are declared REAL ARRA(lShare[148],LBUFF).
C               Setting lShare[148] to 1 is like making these arrays
C               one dimensional.
      INTEGER cacheLineSize
      INTEGER lShare1
      INTEGER lShare4
      INTEGER lShare8
      PARAMETER ( cacheLineSize = 256 )
      PARAMETER ( lShare1 =  cacheLineSize )
      PARAMETER ( lShare4 =  cacheLineSize/4 )
      PARAMETER ( lShare8 =  cacheLineSize/8 )

CC    MAX_VGS  :: Maximum buffer size for Global Vector Sum
c     INTEGER MAX_VGS
c     PARAMETER ( MAX_VGS = 8192 )

C     ========  EESIZE.h  ========================================

C     Symbolic values
C     precXXXX :: precision used for I/O
      INTEGER precFloat32
      PARAMETER ( precFloat32 = 32 )
      INTEGER precFloat64
      PARAMETER ( precFloat64 = 64 )

C     Real-type constant for some frequently used simple number (0,1,2,1/2):
      Real*8     zeroRS, oneRS, twoRS, halfRS
      PARAMETER ( zeroRS = 0.0D0 , oneRS  = 1.0D0 )
      PARAMETER ( twoRS  = 2.0D0 , halfRS = 0.5D0 )
      Real*8     zeroRL, oneRL, twoRL, halfRL
      PARAMETER ( zeroRL = 0.0D0 , oneRL  = 1.0D0 )
      PARAMETER ( twoRL  = 2.0D0 , halfRL = 0.5D0 )

C     UNSET_xxx :: Used to indicate variables that have not been given a value
      Real*8  UNSET_FLOAT8
      PARAMETER ( UNSET_FLOAT8 = 1.234567D5 )
      Real*4  UNSET_FLOAT4
      PARAMETER ( UNSET_FLOAT4 = 1.234567E5 )
      Real*8     UNSET_RL
      PARAMETER ( UNSET_RL     = 1.234567D5 )
      Real*8     UNSET_RS
      PARAMETER ( UNSET_RS     = 1.234567D5 )
      INTEGER UNSET_I
      PARAMETER ( UNSET_I      = 123456789  )

C     debLevX  :: used to decide when to print debug messages
      INTEGER debLevZero
      INTEGER debLevA, debLevB,  debLevC, debLevD, debLevE
      PARAMETER ( debLevZero=0 )
      PARAMETER ( debLevA=1 )
      PARAMETER ( debLevB=2 )
      PARAMETER ( debLevC=3 )
      PARAMETER ( debLevD=4 )
      PARAMETER ( debLevE=5 )

C     SQUEEZE_RIGHT      :: Flag indicating right blank space removal
C                           from text field.
C     SQUEEZE_LEFT       :: Flag indicating left blank space removal
C                           from text field.
C     SQUEEZE_BOTH       :: Flag indicating left and right blank
C                           space removal from text field.
C     PRINT_MAP_XY       :: Flag indicating to plot map as XY slices
C     PRINT_MAP_XZ       :: Flag indicating to plot map as XZ slices
C     PRINT_MAP_YZ       :: Flag indicating to plot map as YZ slices
C     commentCharacter   :: Variable used in column 1 of parameter
C                           files to indicate comments.
C     INDEX_I            :: Variable used to select an index label
C     INDEX_J               for formatted input parameters.
C     INDEX_K
C     INDEX_NONE
      CHARACTER*(*) SQUEEZE_RIGHT
      PARAMETER ( SQUEEZE_RIGHT = 'R' )
      CHARACTER*(*) SQUEEZE_LEFT
      PARAMETER ( SQUEEZE_LEFT = 'L' )
      CHARACTER*(*) SQUEEZE_BOTH
      PARAMETER ( SQUEEZE_BOTH = 'B' )
      CHARACTER*(*) PRINT_MAP_XY
      PARAMETER ( PRINT_MAP_XY = 'XY' )
      CHARACTER*(*) PRINT_MAP_XZ
      PARAMETER ( PRINT_MAP_XZ = 'XZ' )
      CHARACTER*(*) PRINT_MAP_YZ
      PARAMETER ( PRINT_MAP_YZ = 'YZ' )
      CHARACTER*(*) commentCharacter
      PARAMETER ( commentCharacter = '#' )
      INTEGER INDEX_I
      INTEGER INDEX_J
      INTEGER INDEX_K
      INTEGER INDEX_NONE
      PARAMETER ( INDEX_I    = 1,
     &            INDEX_J    = 2,
     &            INDEX_K    = 3,
     &            INDEX_NONE = 4 )

C     EXCH_IGNORE_CORNERS :: Flag to select ignoring or
C     EXCH_UPDATE_CORNERS    updating of corners during an edge exchange.
      INTEGER EXCH_IGNORE_CORNERS
      INTEGER EXCH_UPDATE_CORNERS
      PARAMETER ( EXCH_IGNORE_CORNERS = 0,
     &            EXCH_UPDATE_CORNERS = 1 )

C     FORWARD_SIMULATION
C     REVERSE_SIMULATION
C     TANGENT_SIMULATION
      INTEGER FORWARD_SIMULATION
      INTEGER REVERSE_SIMULATION
      INTEGER TANGENT_SIMULATION
      PARAMETER ( FORWARD_SIMULATION = 0,
     &            REVERSE_SIMULATION = 1,
     &            TANGENT_SIMULATION = 2 )

C--   COMMON /EEPARAMS_L/ Execution environment public logical variables.
C     eeBootError    :: Flags indicating error during multi-processing
C     eeEndError     :: initialisation and termination.
C     fatalError     :: Flag used to indicate that the model is ended with an error
C     debugMode      :: controls printing of debug msg (sequence of S/R calls).
C     useSingleCpuIO :: When useSingleCpuIO is set, MDS_WRITE_FIELD outputs from
C                       master MPI process only. -- NOTE: read from main parameter
C                       file "data" and not set until call to INI_PARMS.
C     useSingleCpuInput :: When useSingleCpuInput is set, EXF_INTERP_READ
C                       reads forcing files from master MPI process only.
C                       -- NOTE: read from main parameter file "data"
C                          and defaults to useSingleCpuInput = useSingleCpuIO
C     printMapIncludesZeros  :: Flag that controls whether character constant
C                               map code ignores exact zero values.
C     useCubedSphereExchange :: use Cubed-Sphere topology domain.
C     useCoupler     :: use Coupler for a multi-components set-up.
C     useNEST_PARENT :: use Parent Nesting interface (pkg/nest_parent)
C     useNEST_CHILD  :: use Child  Nesting interface (pkg/nest_child)
C     useNest2W_parent :: use Parent 2-W Nesting interface (pkg/nest2w_parent)
C     useNest2W_child  :: use Child  2-W Nesting interface (pkg/nest2w_child)
C     useOASIS       :: use OASIS-coupler for a multi-components set-up.
      COMMON /EEPARAMS_L/
c    &  eeBootError, fatalError, eeEndError,
     &  eeBootError, eeEndError, fatalError, debugMode,
     &  useSingleCpuIO, useSingleCpuInput, printMapIncludesZeros,
     &  useCubedSphereExchange, useCoupler,
     &  useNEST_PARENT, useNEST_CHILD,
     &  useNest2W_parent, useNest2W_child, useOASIS,
     &  useSETRLSTK, useSIGREG
      LOGICAL eeBootError
      LOGICAL eeEndError
      LOGICAL fatalError
      LOGICAL debugMode
      LOGICAL useSingleCpuIO
      LOGICAL useSingleCpuInput
      LOGICAL printMapIncludesZeros
      LOGICAL useCubedSphereExchange
      LOGICAL useCoupler
      LOGICAL useNEST_PARENT
      LOGICAL useNEST_CHILD
      LOGICAL useNest2W_parent
      LOGICAL useNest2W_child
      LOGICAL useOASIS
      LOGICAL useSETRLSTK
      LOGICAL useSIGREG

C--   COMMON /EPARAMS_I/ Execution environment public integer variables.
C     errorMessageUnit    :: Fortran IO unit for error messages
C     standardMessageUnit :: Fortran IO unit for informational messages
C     maxLengthPrt1D :: maximum length for printing (to Std-Msg-Unit) 1-D array
C     scrUnit1      :: Scratch file 1 unit number
C     scrUnit2      :: Scratch file 2 unit number
C     eeDataUnit    :: Unit # for reading "execution environment" parameter file
C     modelDataUnit :: Unit number for reading "model" parameter file.
C     numberOfProcs :: Number of processes computing in parallel
C     pidIO         :: Id of process to use for I/O.
C     myBxLo, myBxHi :: Extents of domain in blocks in X and Y
C     myByLo, myByHi :: that each threads is responsble for.
C     myProcId      :: My own "process" id.
C     myPx          :: My X coord on the proc. grid.
C     myPy          :: My Y coord on the proc. grid.
C     myXGlobalLo   :: My bottom-left (south-west) x-index global domain.
C                      The x-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     myYGlobalLo   :: My bottom-left (south-west) y-index in global domain.
C                      The y-coordinate of this point in for example m or
C                      degrees is *not* specified here. A model needs to
C                      provide a mechanism for deducing that information
C                      if it is needed.
C     nThreads      :: No. of threads
C     nTx, nTy      :: No. of threads in X and in Y
C                      This assumes a simple cartesian gridding of the threads
C                      which is not required elsewhere but that makes it easier
C     ioErrorCount  :: IO Error Counter. Set to zero initially and increased
C                      by one every time an IO error occurs.
      COMMON /EEPARAMS_I/
     &  errorMessageUnit, standardMessageUnit, maxLengthPrt1D,
     &  scrUnit1, scrUnit2, eeDataUnit, modelDataUnit,
     &  numberOfProcs, pidIO, myProcId,
     &  myPx, myPy, myXGlobalLo, myYGlobalLo, nThreads,
     &  myBxLo, myBxHi, myByLo, myByHi,
     &  nTx, nTy, ioErrorCount
      INTEGER errorMessageUnit
      INTEGER standardMessageUnit
      INTEGER maxLengthPrt1D
      INTEGER scrUnit1
      INTEGER scrUnit2
      INTEGER eeDataUnit
      INTEGER modelDataUnit
      INTEGER ioErrorCount(MAX_NO_THREADS)
      INTEGER myBxLo(MAX_NO_THREADS)
      INTEGER myBxHi(MAX_NO_THREADS)
      INTEGER myByLo(MAX_NO_THREADS)
      INTEGER myByHi(MAX_NO_THREADS)
      INTEGER myProcId
      INTEGER myPx
      INTEGER myPy
      INTEGER myXGlobalLo
      INTEGER myYGlobalLo
      INTEGER nThreads
      INTEGER nTx
      INTEGER nTy
      INTEGER numberOfProcs
      INTEGER pidIO

CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***
CBOP
C     !ROUTINE: FFIELDS.h
C     !INTERFACE:
C     include "FFIELDS.h"
C     !DESCRIPTION:
C     \bv
C     *==========================================================*
C     | FFIELDS.h
C     | o Model forcing fields
C     *==========================================================*
C     | More flexible surface forcing configurations are
C     | available via, e.g., pkg/exf
C     *==========================================================*
C     \ev
CEOP
C
C     fu    :: Zonal surface wind stress in N/m^2
C              > 0 for increase in uVel, which is west to
C                  east for cartesian and spherical polar grids
C              Typical range: -0.5 < fu < 0.5
C              Southwest C-grid U point
C
C     fv    :: Meridional surface wind stress in N/m^2
C              > 0 for increase in vVel, which is south to
C                  north for cartesian and spherical polar grids
C              Typical range: -0.5 < fv < 0.5
C              Southwest C-grid V point
C
C     EmPmR :: Net upward freshwater flux in kg/m2/s
C              EmPmR = Evaporation - precipitation - runoff
C              > 0 for increase in salt (ocean salinity)
C              Typical range: -1e-4 < EmPmR < 1e-4
C              Southwest C-grid tracer point
C           NOTE: for backward compatibility EmPmRfile is specified in
C                 m/s when using external_fields_load.F.  It is converted
C                 to kg/m2/s by multiplying by rhoConstFresh.
C
C  saltFlux :: Net upward salt flux in g/kg.kg/m^2/s = g/m^2/s
C              flux of Salt taken out of the ocean per time unit (second).
C              Note: only used when salty sea-ice forms or melts.
C              > 0 for decrease in SSS.
C              Southwest C-grid tracer point
C
C     Qnet  :: Net upward surface heat flux (including shortwave) in W/m^2
C              Qnet = latent + sensible + net longwave + net shortwave
C              > 0 for decrease in theta (ocean cooling)
C              Typical range: -250 < Qnet < 600
C              Southwest C-grid tracer point
C
C     Qsw   :: Net upward shortwave radiation in W/m^2
C              Qsw = - ( downward - ice and snow absorption - reflected )
C              > 0 for decrease in theta (ocean cooling)
C              Typical range: -350 < Qsw < 0
C              Southwest C-grid tracer point
C
C     SST   :: Sea surface temperature in degrees C for relaxation
C              Southwest C-grid tracer point
C
C     SSS   :: Sea surface salinity in g/kg for relaxation
C              Southwest C-grid tracer point
C
C     lambdaThetaClimRelax :: Inverse time scale for SST relaxation ( 1/s ).
C
C     lambdaSaltClimRelax  :: Inverse time scale for SSS relaxation ( 1/s ).

C     phiTide2d :: vertically uniform (2d-map), time-dependent geopotential
C                  anomaly (e.g., tidal forcing); Units are m^2/s^2
C     pLoad :: for the ocean:      atmospheric pressure anomaly (relative to
C                                   "surf_pRef") at z=eta
C                Units are           Pa=N/m^2
C              for the atmosphere (hack): geopotential anomaly of the orography
C                Units are           m^2/s^2
C     sIceLoad :: sea-ice loading, expressed in Mass of ice+snow / area unit
C                Units are           kg/m^2
C              Note: only used with Sea-Ice & RealFreshWater formulation
C     addMass  :: source (<0: sink) of fluid in the domain interior
C                 (generalisation of oceanic real fresh-water flux)
C                Units are           kg/s  (mass per unit of time)
C     frictionHeating :: heating caused by friction and momentum dissipation
C                Units are           in W/m^2 (thickness integrated)
C     eddyPsiX -Zonal Eddy Streamfunction in m^2/s used in taueddy_external_forcing.F
C     eddyPsiY -Meridional Streamfunction in m^2/s used in taueddy_external_forcing.F
C     EfluxY - y-component of Eliassen-Palm flux vector
C     EfluxP - p-component of Eliassen-Palm flux vector

      COMMON /FFIELDS_fu/ fu
      COMMON /FFIELDS_fv/ fv
      COMMON /FFIELDS_Qnet/ Qnet
      COMMON /FFIELDS_Qsw/ Qsw
      COMMON /FFIELDS_EmPmR/ EmPmR
      COMMON /FFIELDS_saltFlux/ saltFlux
      COMMON /FFIELDS_SST/ SST
      COMMON /FFIELDS_SSS/ SSS
      COMMON /FFIELDS_lambdaThetaClimRelax/ lambdaThetaClimRelax
      COMMON /FFIELDS_lambdaSaltClimRelax/ lambdaSaltClimRelax
      COMMON /FFIELDS_phiTide/ phiTide2d
      COMMON /FFIELDS_pLoad/ pLoad
      COMMON /FFIELDS_sIceLoad/ sIceLoad

      Real*8  fu       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  fv       (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qnet     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qsw      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  EmPmR    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  saltFlux (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SST      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SSS      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  lambdaThetaClimRelax(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  lambdaSaltClimRelax(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  phiTide2d(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  pLoad    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  sIceLoad (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C- jmc: commented out until corresponding (ghost-like) code apparition
C     dQdT  :: Thermal relaxation coefficient in W/m^2/degrees
C              Southwest C-grid tracer point
c     COMMON /FFIELDS_dQdT/ dQdT
c     Real*8  dQdT   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
c#ifdef ALLOW_EP_FLUX
c     COMMON /FFIELDS_eflux/ EfluxY,EfluxP
c     Real*8  EfluxY (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
c     Real*8  EfluxP (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
c#endif


C     loadedRec     :: time-record currently loaded (in temp arrays *[1])
C     taux[0,1]     :: Temp. for zonal wind stress
C     tauy[0,1]     :: Temp. for merid. wind stress
C     Qnet[0,1]     :: Temp. for heat flux
C     EmPmR[0,1]    :: Temp. for fresh water flux
C     saltFlux[0,1] :: Temp. for isurface salt flux
C     SST[0,1]      :: Temp. for theta climatalogy
C     SSS[0,1]      :: Temp. for theta climatalogy
C     Qsw[0,1]      :: Temp. for short wave component of heat flux
C     pLoad[0,1]    :: Temp. for atmospheric pressure at z=eta
C     [0,1]         :: End points for interpolation

      COMMON /FFIELDS_I/ loadedRec
      INTEGER loadedRec(nSx,nSy)

      COMMON /TDFIELDS/
     &                 taux0, tauy0, Qnet0, EmPmR0, SST0, SSS0,
     &                 taux1, tauy1, Qnet1, EmPmR1, SST1, SSS1,
     &                 saltFlux0, saltFlux1
     &               , Qsw0, Qsw1
     &               , pLoad0, pLoad1

      Real*8  taux0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  tauy0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qnet0    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  EmPmR0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  saltFlux0(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SST0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SSS0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  taux1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  tauy1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qnet1    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  EmPmR1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  saltFlux1(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SST1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  SSS1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qsw0     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  Qsw1     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  pLoad0   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  pLoad1   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     surfaceForcingU     units are  r_unit.m/s^2 (=m^2/s^2 if r=z)
C                -> usage in gU:     gU = gU + surfaceForcingU/drF [m/s^2]
C     surfaceForcingV     units are  r_unit.m/s^2 (=m^2/s^-2 if r=z)
C                -> usage in gU:     gV = gV + surfaceForcingV/drF [m/s^2]
C
C     surfaceForcingS     units are  r_unit.g/kg/s (=g/kg.m/s if r=z)
C            - EmPmR * S_surf plus salinity relaxation*drF(1)
C                -> usage in gS:     gS = gS + surfaceForcingS/drF [g/kg/s]
C
C     surfaceForcingT     units are  r_unit.Kelvin/s (=Kelvin.m/s if r=z)
C            - Qnet (+Qsw) plus temp. relaxation*drF(1)
C                -> calculate        -lambda*(T(model)-T(clim))
C            Qnet assumed to be net heat flux including ShortWave rad.
C                -> usage in gT:     gT = gT + surfaceforcingT/drF [K/s]
C     adjustColdSST_diag :: diagnostic field for how much too cold (below
C              Tfreezing) SST has been adjusted (with allowFreezing=T).
C              > 0 for increase of SST (up to Tfreezing).
C              Units are r_unit.K/s (=Kelvin.m/s if r=z).
C        Note: 1) allowFreezing option is a crude hack to fix too cold SST that
C              results from missing seaice component. It should never be used
C              with any seaice component, neither current seaice pkg (pkg/seaice
C              or pkg/thsice) nor a seaice component from atmos model when
C              coupled to it.
C              2) this diagnostic is currently used by KPP package (kpp_calc.F
C              and kpp_transport_t.F) although it is not very clear it should.

      COMMON /SURFACE_FORCING/
     &                         surfaceForcingU,
     &                         surfaceForcingV,
     &                         surfaceForcingT,
     &                         surfaceForcingS,
     &                         adjustColdSST_diag
      Real*8  surfaceForcingU   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  surfaceForcingV   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  surfaceForcingT   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  surfaceForcingS   (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  adjustColdSST_diag(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C     botDragU :: bottom stress (for diagnostics), Zonal component
C                Units are N/m^2 ;   > 0 increase uVel @ bottom
C     botDragV :: bottom stress (for diagnostics), Merid. component
C                Units are N/m^2 ;   > 0 increase vVel @ bottom
      COMMON /FFIELDS_bottomStress/ botDragU, botDragV
      Real*8  botDragU (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  botDragV (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     === Routine arguments ===
C     INPUT:
C     UG      :: thermal wind of atmosphere
C     TSURF   :: ocean surface temperature in Kelvin
C     bi,bj   :: loop indices
C     myTime  :: Simulation time
C     myIter  :: Simulation timestep number
C     myThid  :: Thread no. that called this routine.
C     OUTPUT:
C     netHeatFlux :: net surface heat flux over open water or under ice
C     SWHeatFlux  :: short wave heat flux over open water or under ice
      Real*8 UG         (1:sNx,1:sNy)
      Real*8 TSURF      (1:sNx,1:sNy)
      Real*8 netHeatFlux(1:sNx,1:sNy)
      Real*8 SWHeatFlux (1:sNx,1:sNy)
      Real*8 myTime
      INTEGER bi, bj, myIter, myThid
CEndOfInterface

C     === Local variables ===
C     i,j - Loop counters
      INTEGER i, j

C NOW DETERMINE OPEN WATER HEAT BUD. ASSUMING TSURF=WATER TEMP.
C WATER ALBEDO IS ASSUMED TO BE THE CONSTANT SEAICE_waterAlbedo
      DO j=1,sNy
       DO i=1,sNx
        netHeatFlux(i,j) = Qnet(i,j,bi,bj)
        SWHeatFlux (i,j) =  Qsw(i,j,bi,bj)
       ENDDO
      ENDDO

      RETURN
      END
