










C CPP options file for salt_plume package
C Use this file for selecting options within the salt_plume package








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


C Place CPP define/undef flag here

C SALT_PLUME_IN_LEADS
C   Motivation: As ice concentration AREA -> 1, leads occur -> ice
C     production is no longer uniform in grid box -> assumptions
C     which motivate KPP no longer holds -> treat overturn more
C     realistic with this flag.
C   if defined: Activate pkg/salt_plume only when seaice AREA exceeds
C               a certain value representative of lead opening AND only
C               if seaice growth dh is from atmospheric cooling.
C   if undefined: Activate pkg/salt_plume whenever seaice forms.
C                 This is the default of pkg/salt_plume.


CBOP
C     !ROUTINE: SALT_PLUME_FRAC
C     !INTERFACE:
      SUBROUTINE SALT_PLUME_FRAC(
     I                  imax, fact,SPDepth,
     U                  plumek,
     I                  myTime, myIter, myThid )

C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE SALT_PLUME_FRAC
C     | o Compute saltplume penetration.
C     *==========================================================*
C     | Compute fraction of saltplume (flux) penetrating to
C     | specified depth, plumek, due to rejected salt
C     | during freezing.
C     | For example, if surface value is Saltplume0,
C     | and each level gets equal fraction 1/5 down to SPDepth=5,
C     | SALT_PLUME_FRAC will report plumek = 1/5,2/5,3/5,4/5,5/5 on
C     | output if input plumek = 1,2,3,4,5. Else, output plumek = 0.
C     | Reference : Duffy et al, (GRL 1999)
C     |
C     | =====
C     | Written by   : ATN (based on SWFRAC)
C     | Date         : Sep 13, 2007
C     | Modified     : Mar 16, 2014 by atn to improve/clean up
C     | -> replace 1-[cummulative plumefrac] code which was based
C     |    on swfrac with cleaner [cummulative plumefrac] on output
C     |    in order to avoid 1-[1-[cummulative_plumefrac]] whenever
C     |    SALT_PLUME_FRAC is called and used from outside.
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
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


C--   SALT_PLUME parameters
C     Find surface where the potential density (ref.lev=surface) is
C     larger than surface density plus SaltPlumeCriterion.

C     SaltPlumeSouthernOcean: TRUE  = apply salt plume globally
C                             FALSE = apply salt plume in Arctic Ocean only
      LOGICAL SaltPlumeSouthernOcean
      COMMON /SALT_PLUME_PARAMS_L/ SaltPlumeSouthernOcean

C     CriterionType: 1=delta_rho, 2=drhodz, default is 1
C     PlumeMethod: method of distributing salt plume vertically
C       1=power, 2=exp, 3=overshoot, 5=dump_at_top, 6=reverse of 1
C       default is 1
C     Npower: choices of distributing salt uniformly (0), linear (1),
C       or higher power (Npower>1); default is 0 when PlumeMethod = 1
      INTEGER CriterionType, PlumeMethod
     & , Npower
      COMMON /SALT_PLUME_PARAMS_I/ CriterionType, PlumeMethod, Npower

C     SaltPlumeCriterion
C       for CriterionType=1, default is 0.4 kg/m^3 of Duffy et al 1999
C       for CriterionType=2, default is 0.005 kg/m^3/m
C     SPovershoot: overshooting depth of penetrating salt plume,
C       so that 1.0 = no-overshoot, 1.2 = 20% overshoot.
C       default is 1.0
C     SPsalFRAC: fraction of the salt by-product of seaice growth (not melt) that
C       will be re-distributed vertically according to the salt_plume_frac.F
C       Its default is 1. (for 100% effect), and its range is [0. 1.]
C     SPinflectionPoint: the inflection point of a nonlinear function
C       f(AREA) controlling saltPlumeFlux. f(AREA) is a logistic curve
C       (sigmoid) with range [0. 1.] and f(SPinflectionPoint) == 0.5.
C       Usage: pkg/salt_plume activates when AREA>=SPinflectionPoint.
C       To assure only narrow leads generate plumes:
C       set SPinflectionPoint >= 0.8.
C     SPalpha :: fraction of grid volume designated to be brine, [0. 1.]
C       If grid cell 18km x 18km x 10m, take SPalpha=0.001 gives
C       volume of 0.001*drF(1)*[dx*dy] of brine. Thus SPbrineSalt
C       can be calc as adding SaltPlumeFlux into this fractional vol.
C       Default: 0.008 -> SPbrineSalt ~37 if SSS is ~32.
C     SPbrineSconst :: salinity of brine pocket (g/kg)

      Real*8 SaltPlumeCriterion, SPovershoot
     &   , SPsalFRAC
      COMMON /SALT_PLUME_PARAMS_R/
     &   SPsalFRAC, SaltPlumeCriterion, SPovershoot
C--   SALT_PLUME 2-dim. fields
C     SaltPlumeDepth :: depth of penetration of salt plumes
C                       rejected during sea ice growth
C     saltPlumeFlux :: Net downward salt flux in g/m^2/s
C              Note: only used when salty sea-ice forms.
C              > 0 for increasing in SSS.
C              Southwest C-grid tracer point
C     dSPvolSurf2kLev :: downward volume frac from klev=1 associated w/ saltPlumeFlux
C     dSPvolBelow2kLev:: upward volume frac from grid below (RETIRED)
C     dSPvolkLev2Above:: upward volume frac to grid above
C     SPbrineVolFlux  :: brine Vol associated w/ SPbrineSconst & saltPlumeFlux
C     SPforcingS      :: 3D forcingS associated w/ saltPlumeFlux in g/m^2/s
C     SPforcingT      :: 3D forcingT associated w/ saltPlumeFlux in W/m^2
C     SPkBottom       :: bottom kLev associated with SaltPlumeDepth

      Real*8 SaltPlumeDepth (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      Real*8  saltPlumeFlux (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      COMMON /DYNVARS_SALT_PLUME/ SaltPlumeDepth
      COMMON /FFIELDS_saltPlumeFlux/ saltPlumeFlux

C     !INPUT/OUTPUT PARAMETERS:
C     input arguments
C     imax    :: number of vertical grid points
C     fact    :: scale  factor to apply to depth array
C     SPDpeth :: corresponding SaltPlumeDepth(i,j) at this grid point
C     myTime  :: Current time in simulation
C     myIter  :: Current iteration number in simulation
C     myThid  :: My Thread Id. number
      INTEGER imax
      Real*8     fact
      Real*8     myTime
      INTEGER myIter
      INTEGER myThid
C     input/output arguments
C     plumek :: on input: vertical depth for desired plume fraction
C               (fact*plumek) is negative distance (m) from surface
C     plumek :: on output: saltplume contribution fraction
      Real*8     plumek(imax), SPDepth(imax)
CEOP

C     !LOCAL VARIABLES:
      Real*8 facz, dd, dd20
      INTEGER i, Npowerloc
c     INTEGER kk
      Real*8     one, two, three, S, So, zero
      parameter( one = 1.D0, two = 2.D0, three = 3.D0 )
      parameter( zero = 0.D0 )
C     This is an abbreviation of 1./(exp(1.)-1.)
      Real*8     recip_expOneM1
      parameter( recip_expOneM1 = 0.581976706869326343 )

      DO i = 1,imax
       facz = abs(fact*plumek(i))
         Npowerloc = Npower

       IF (SPDepth(i).GE.facz .and. SPDepth(i) .GT. zero) THEN

C     Default: uniform distribution, PlumeMethod=1, Npower=0
        IF (PlumeMethod .EQ. 1) THEN
         dd20 = (abs(SPDepth(i)))
c#ifdef TARGET_NEC_SX
C     This rewriting (originally for TARGET_NEC_SX) is better for all platforms
         IF ( dd20 .GT. zero ) THEN
          S   = (facz/dd20)
C     crazy attempt to make the code faster and raise S to (Npower+1)
          IF (Npowerloc .GT. 0) S = S*S**Npowerloc
         ELSE
          S = zero
         ENDIF
         plumek(i) = max(zero,S)
c#else
c        S  = one                  !input depth temp
c        So = one
c        DO kk=1,Npowerloc+1
c         S  = facz*S              !raise to the Npowerloc+1
c         So = dd20*So
c        ENDDO
c        plumek(i) = max(zero,S/So)
c#endif 

        ELSEIF (PlumeMethod .EQ. 2) THEN !exponential distribution
         dd = abs(SPDepth(i))
         S  = exp(facz/dd)-one
         So = recip_expOneM1       !So = exp(one)-one
         plumek(i) = max(zero,S*So)

C     PlumeMethod = 3, distribute salt LINEARLY between SPDepth and
C     SPDepth/SPovershoot
C     (1-SPovershoot)percent has already been taken into account in
C     SPDepth calculation, i.e., SPDepth = SPovershoot*SPDepth.
        ELSEIF (PlumeMethod .EQ. 3) THEN !overshoot 20%
         dd20 = (abs(SPDepth(i)))
         dd   = dd20/SPovershoot
         So=dd20-dd
         S =facz-dd
         IF( (facz.GE.dd).AND.(facz.LT.dd20) ) THEN
          plumek(i) = max(zero,S/So)
         ELSE
          plumek(i) = zero
         ENDIF

C     PlumeMethod = 5, dumping all salt at the top layer
        ELSEIF (PlumeMethod .EQ. 5) THEN
         dd   = zero
         dd20 = one
         IF( (facz.GE.dd).AND.(facz.LT.dd20) ) THEN
          plumek(i) = zero
         ELSE
          plumek(i) = one
         ENDIF
        ELSEIF (PlumeMethod .EQ. 6) THEN
C     PLumeMethod = 6, currently only works for Npower = 1 and 2.
         dd20 = (abs(SPDepth(i)))
c#ifdef TARGET_NEC_SX
C     This rewriting (originally for TARGET_NEC_SX) is better for all platforms
         IF ( dd20 .GT. zero ) THEN
          S  = (facz/dd20)
C     crazy attempt to make the code faster and raise S to (Npower+1)
          IF (Npowerloc .GT. 0) S = S*S**Npowerloc
          So = 1.D0/dd20
         ELSE
          S  = zero
          So = zero
         ENDIF
         IF(Npowerloc.EQ.1) THEN   !Npower=1
          plumek(i) = max(zero,two*So*facz-S)
         ELSE                      !Npower=2
          plumek(i) = max(zero,
     &         three*So*facz - three*So*So*facz*facz + S)
         ENDIF
c#else
c        S  = one                  !input depth temp
c        So = one
c        DO kk=1,Npowerloc+1
c         S  = facz*S              !raise to the Npower+1
c         So = dd20*So
c        ENDDO
c        IF(Npowerloc.EQ.1) THEN   !Npower=1
c         plumek(i) = max(zero,two/dd20*facz-S/So)
c        ELSE                      !Npower=2
c         plumek(i) = max(zero,
c    &         three/dd20*facz - three/(dd20*dd20)*facz*facz + S/So)
c        ENDIF
c#endif 

catn: 15.Mar.2014
catn: this is a new method by atn. After fixing adjoint compiling error,
catn: will switch this on.  Currently commenting out for purpose of
catn: comparing old (1-abc) vs new (abc) way of coding
c        ELSEIF (PlumeMethod .EQ. 7) THEN
cC     PLumeMethod = 7, dump an offset parabolla with more salt at surface
cC        tapered to zero at depth SPDepth/2, then increased until SPDepth
cC        need input SPDepth, NPower = percentage of SPDepth
cC        Ex: if Npower = 10 -> (10/2)=5% of SPDepth
cC        NPower can be negative integer here.
cC        0 -> parabola centered at SPDepth/2;
cC        + -> parabola offset, salt @ surface < @ SPDepth
cC        - -> parabola offset, salt @ surface > @ SPDepth
cC        S and So are dummy variables
c         dd   = (abs(SPDepth(i)))
c         dd20 = dd*(one/two-Npower/200.D0)
c         So   = (dd*dd*dd/three)
c     &            -(dd*dd)      *dd20
c     &            +(dd)         *dd20*dd20
c         S    = (facz*facz *facz/three)
c     &             - (facz*facz)*dd20
c     &             + (facz)     *dd20*dd20
c         plumek(i) = max(zero,(S/So))
c
        ELSE
         plumek(i) = one
         WRITE(*,*) 'salt_plume_frac: PLumeMethod =', PLumeMethod,
     &        'not implemented'
         STOP 'ABNORMAL in S/R SALT_PLUME_FRAC'
        ENDIF
       ELSE
        plumek(i) = one
       ENDIF
      ENDDO


      RETURN
      END
