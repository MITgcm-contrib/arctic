

















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


C     Package-specific Options & Macros go here

C allow to define specific regions and the corresponding mask ;
C  used to perform regional statistics over a limited area

C allow to stop & restart at any time (i.e. not at a multiple of
C  the diagnostics frequency) reading diagnostics storage arrays
C  from pickup file.
C Note: Use with cautious since it does not work for all restart
C  cases (e.g., changing data.diagnostics).


CEH3 ;;; Local Variables: ***
CEH3 ;;; mode:fortran ***
CEH3 ;;; End: ***

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP 0
C     !ROUTINE: DIAGNOSTICS_INTERP_P2P

C     !INTERFACE:
      SUBROUTINE DIAGNOSTICS_INTERP_P2P(
     O                       qprs,
     I                       qinp,pkz,pksrf,pktop,pk,
     I                       undef, pInc,ijm,lm, myThid )

C     !DESCRIPTION:
C***********************************************************************
C
C PURPOSE
C   To interpolate an arbitrary quantity to Specified Pressure Levels
C
C INPUT
C   QINP .. QINP (ijm,lm) Arbitrary Input Quantity
C   PKZ ... PKZ  (ijm,lm) Pressure to the Kappa at Input Levels
C   PKSRF . PKSRF(ijm) Surface Pressure to the Kappa
C   PKTOP . Pressure to the Kappa at Input-Level-Edge (1) (top of model)
C   PK .... Output Pressure to the Kappa Level (mb)
C   pInc .. if True, assume pressure increases with level index
C   IJM ... Horizontal Dimension of Input
C   LM .... Vertical  Dimension of Input
C
C OUTPUT
C   QPRS .. QPRS (ijm) Arbitrary Quantity at Pressure p
C
C NOTE
C   Quantity is interpolated Linear in P**Kappa.
C   Between  PTOP**Kappa and PKZ(1),  quantity is extrapolated.
C   Between PKSRF**Kappa and PKZ(LM), quantity is extrapolated.
C   Undefined Input quantities are not used.
C   Finally: This routine assumes that pressure levels are counted
C            top down -- ie, level 1 is the top, level lm is the bottom
C
C***********************************************************************
C     !USES:
      IMPLICIT NONE

C     !INPUT PARAMETERS:
      INTEGER  ijm,lm,myThid
      Real*8  qinp (ijm,lm)
      Real*8  pkz  (ijm,lm)
      Real*8  pksrf(ijm)
      Real*8  pktop,pk
      Real*8  undef
      LOGICAL pInc

C     !OUTPUT PARAMETERS:
      Real*8  qprs (ijm)
CEOP

C     !LOCAL VARIABLES:
      INTEGER  i,l
      Real*8  pkmin,pkmax,temp

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

c Initialize to UNDEFINED
c -----------------------
      DO i=1,ijm
       qprs(i) = undef
      ENDDO

      IF ( pInc ) THEN
C---  Case where Levels are orderd by increasing pressure

c Interpolate to Pressure Between Input Levels
c --------------------------------------------
      DO L=1,lm-1
       pkmin = pkz(1,L)
       pkmax = pkz(1,L+1)

       DO i=2,ijm
         IF ( pkz(i,L)  .LT.pkmin ) pkmin = pkz(i,L)
         IF ( pkz(i,L+1).GT.pkmax ) pkmax = pkz(i,L+1)
       ENDDO

       IF ( pk.LE.pkmax .AND. pk.GE.pkmin ) THEN
         DO i=1,ijm
           IF ( pk.GE.pkz(i,L) .AND. pk.LE.pkz(i,L+1) ) THEN
             temp = ( pk-pkz(i,L) ) / ( pkz(i,L+1)-pkz(i,L) )

             IF     ( qinp(i,L)  .NE.undef  .AND.
     &                qinp(i,L+1).NE.undef ) THEN
               qprs(i) = qinp(i,L+1)*temp + qinp(i,L)*(1.-temp)
             ELSEIF ( qinp(i,L+1).NE.undef  .AND. temp.GE.0.5 ) THEN
               qprs(i) = qinp(i,L+1)
             ELSEIF ( qinp(i,L)  .NE.undef  .AND. temp.LE.0.5 ) THEN
               qprs(i) = qinp(i,L)
             ENDIF
           ENDIF
         ENDDO
       ENDIF

      ENDDO

      DO i=1,ijm
c Extrapolate to Pressure between Ptop and Highest Input Level
c ----------------------------------------------------------
       IF ( pk.LE.pkz(i,1) .AND. pk.GE.pktop ) THEN
         temp = ( pk-pkz(i,1) ) / ( pkz(i,2)-pkz(i,1) )

         IF     ( qinp(i,1).NE.undef  .AND.
     &           qinp(i,2).NE.undef ) THEN
           qprs(i) = qinp(i,2)*temp + qinp(i,1)*(1.-temp)
         ELSEIF ( qinp(i,1).NE.undef ) THEN
           qprs(i) = qinp(i,1)
         ENDIF

       ENDIF

c Extrapolate to Pressure between Psurf and Lowest Input Level
c ------------------------------------------------------------
       IF ( pk.GE.pkz(i,lm) .AND. pk.LE.pksrf(i) ) THEN
         temp = ( pk-pkz(i,lm) ) / ( pkz(i,lm-1)-pkz(i,lm) )

         IF     ( qinp(i,lm)  .NE.undef  .AND.
     &            qinp(i,lm-1).NE.undef ) THEN
            qprs(i) = qinp(i,lm-1)*temp + qinp(i,lm)*(1.-temp)
         ELSEIF ( qinp(i,lm)  .NE.undef ) THEN
            qprs(i) = qinp(i,lm)
         ENDIF

       ENDIF
      ENDDO

      ELSE
C---  Case where Levels are orderd by decreasing pressure

c Interpolate to Pressure Between Input Levels
c --------------------------------------------
      DO L=1,lm-1
       pkmin = pkz(1,L+1)
       pkmax = pkz(1,L)

       DO i=2,ijm
         IF ( pkz(i,L+1).LT.pkmin ) pkmin = pkz(i,L+1)
         IF ( pkz(i,L)  .GT.pkmax ) pkmax = pkz(i,L)
       ENDDO

       IF ( pk.LE.pkmax .AND. pk.GE.pkmin ) THEN
         DO i=1,ijm
           IF ( pk.LE.pkz(i,L) .AND. pk.GE.pkz(i,L+1) ) THEN
             temp = ( pk-pkz(i,L) ) / ( pkz(i,L+1)-pkz(i,L) )

             IF     ( qinp(i,L)  .NE.undef  .AND.
     &                qinp(i,L+1).NE.undef ) THEN
               qprs(i) = qinp(i,L+1)*temp + qinp(i,L)*(1.-temp)
             ELSEIF ( qinp(i,L+1).NE.undef  .AND. temp.GE.0.5 ) THEN
               qprs(i) = qinp(i,L+1)
             ELSEIF ( qinp(i,L)  .NE.undef  .AND. temp.LE.0.5 ) THEN
               qprs(i) = qinp(i,L)
             ENDIF
           ENDIF
         ENDDO
       ENDIF

      ENDDO

      DO i=1,ijm
c Extrapolate to Pressure between Ptop and Highest Input Level
c ----------------------------------------------------------
       IF ( pk.LE.pkz(i,lm) .AND. pk.GE.pktop ) THEN
         temp = ( pk-pkz(i,lm) ) / ( pkz(i,lm-1)-pkz(i,lm) )

         IF     ( qinp(i,lm)  .NE.undef  .AND.
     &            qinp(i,lm-1).NE.undef ) THEN
            qprs(i) = qinp(i,lm-1)*temp + qinp(i,lm)*(1.-temp)
         ELSEIF ( qinp(i,lm)  .NE.undef ) THEN
            qprs(i) = qinp(i,lm)
         ENDIF

       ENDIF

c Extrapolate to Pressure between Psurf and Lowest Input Level
c ------------------------------------------------------------
       IF ( pk.GE.pkz(i,1) .AND. pk.LE.pksrf(i) ) THEN
         temp = ( pk-pkz(i,1) ) / ( pkz(i,2)-pkz(i,1) )

         IF     ( qinp(i,1).NE.undef  .AND.
     &            qinp(i,2).NE.undef ) THEN
           qprs(i) = qinp(i,2)*temp + qinp(i,1)*(1.-temp)
         ELSEIF ( qinp(i,1).NE.undef ) THEN
           qprs(i) = qinp(i,1)
         ENDIF

       ENDIF
      ENDDO

C---  End case increasing/decreasing pressure
      ENDIF

      RETURN
      END
