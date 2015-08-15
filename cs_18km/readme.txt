Instructions for running optimized (A1) simulation from Nguyen et al (2011)
http://ecco2.org/manuscripts/2011/NguyenJGR2011.pdf

===============

DM 9-Jul-2015
Updated for checkpoint65m (2015/06/15) on pleiades

-> get a copy of ftp://ecco2.jpl.nasa.gov/data6/arctic/run_template_cube81
   also available on pfe:/nobackup/dmenemen/arctic18km/run_template_cube81
   and harbor.mit.edu:/net/nares/raid8/ecco-shared/arctic18km/run_template_cube81

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/atmos/jra25
   also available on pfe:/nobackup/hzhang1/forcing/jra25

cvs co MITgcm_contrib/arctic/cs_18km
cvs co -r checkpoint65m MITgcm_code
cd MITgcm
mkdir build run
cd build
module purge
module load comp-intel/2015.0.090 mpi-sgi/mpt.2.11r13 netcdf/4.0
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas -mods ../../MITgcm_contrib/arctic/cs_18km/code
make depend
make -j 16

cd ../run
cp ../../MITgcm_contrib/arctic/cs_18km/input/* .
rm data.obcs_*
ln -sf ../../run_template_cube81/*.bin .
ln -sf ../../run_template_cube81/BATHY_cube81_420x384_arctic .
ln -sf ../../run_template_cube81/WOA05_* .
ln -sf ../../run_template_cube81/*cube81 .
ln -sf ../../run_template_cube81/+ob1979_2014_merge/* .
ln -sf ../../jra25/jra25_*_199[1-9] .
ln -sf ../../jra25/jra25_*_2* .
rm jra25_rain_n*
ln -sf ../build/mitgcmuv .
qsub run_arctic_cs_18km_nas

-> ftp://ecco2.jpl.nasa.gov/data1/arctic/arctic_18km/GRACErunoff/matlab/
   ftp://ecco2.jpl.nasa.gov/data1/arctic/arctic_18km/GRACErunoff/matlab/lookat.m
   contains some routines for looking at this solution, in particular:
   ftp://ecco2.jpl.nasa.gov/data1/arctic/arctic_18km/GRACErunoff/figs/FIG9.jpg
   is a recreation of Figure 9 in Nguyen et al (2011):
   the results are similar but not identical, e.g., thinner sea ice near
   Canadian Archipelago for TBD reasons.

===============

ATN 23-Dec-2013

SEAICE_OPTIONS.h_atn_May2011: Paper: Nguyen et al 2011

Cpp_OPTIONS.h_oldcode_before2012: code similar to what is used in Nguyen et al 2011, but compatible with codes around Sep2011
SEAICE_OPTIONS.h_Ian_Sep2011: Equivalence to above but updated to be compatible with Ian's code on that date

After 2011, pkg/seaice has updated significantly with the merge of Ian's code into the main branch.
These followings _phAug2012 files include appropriate flag adjustments + param name change for new code.
They should reproduce almost identical results to Nguyen et al 2011 paper, except for whatever that has
improved in ocean/sea-ice physics.  An offline comparison done by ATN on Aug2012 showed the two runs Aug2012
and Nguyen et al 2011 papers are VERY close to each other, thus ATN considered the two solutions are equivalent.

CPP_OPTIONS.h_phAug2012
SEAICE_OPTIONS.h_phAug2012

To run this set up, use binary inputs from:
ftp://ecco2.jpl.nasa.gov/data6/arctic/run_template_cube81/
obcs for the years 1979--2014:  +ob1979_2014_merge

There are two bathymetries:
1) original from Nguyen et al 2011: BATHY_cube81_420x384_arctic
2) an updated version to blank out some CAA passages (set bathymetry to zero)
  because 18km can not resolve these passages resutling in excessive sea ice thickness
