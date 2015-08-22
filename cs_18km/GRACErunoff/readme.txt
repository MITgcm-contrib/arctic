Instructions for running optimized (A1) simulation from Nguyen et al (2011)
http://ecco2.org/manuscripts/2011/NguyenJGR2011.pdf
but using climatological monthly mean river runoff from Whitefield et al (2015)
http://ecco2.org/manuscripts/2015/Whitefield2015.pdf
for 1992-2002 followed GRACE-derived runoff for 2003-2014

===============

DM 22-Aug-20215
Instructions for checkpoint65m (2015/06/15) on pleiades
runARDAT: Whitefield et al (2015) river discharge and temperature
runGRACE: includes interranual discharge variability derived from GRACE
          for 2002-2014 for six largest rivers

-> get a copy of ftp://ecco2.jpl.nasa.gov/data6/arctic/run_template_cube81
   also available on pfe:/nobackup/dmenemen/arctic18km/run_template_cube81
   and harbor.mit.edu:/net/nares/raid8/ecco-shared/arctic18km/run_template_cube81

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/atmos/jra25
   also available on pfe:/nobackup/hzhang1/forcing/jra25

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/runoff/ARDAT
   also available on pfe:/nobackup/dmenemen/forcing/runoff/ARDAT

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/runoff/GRACE
   also available on pfe:/nobackup/dmenemen/forcing/runoff/GRACE

cvs co MITgcm_contrib/arctic/cs_18km
cvs co -r checkpoint65m MITgcm_code
cd MITgcm
mkdir buildARDAT runARDAT runGRACE
cd buildARDAT
module purge
module load comp-intel/2015.0.090 mpi-sgi/mpt.2.11r13 netcdf/4.0
cp ../../MITgcm_contrib/arctic/cs_18km/code/* .
cp ../../MITgcm_contrib/arctic/cs_18km/GRACE/EXF_OPTIONS.h .
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas
make depend
make -j 16

cd ../runARDAT
cp ../../MITgcm_contrib/arctic/cs_18km/input/* .
rm data.obcs_*
cp ../../MITgcm_contrib/arctic/cs_18km/GRACErunoff/data.exf .
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1991
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1992
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1993
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1994
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1995
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1996
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1997
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1998
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_1999
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2000
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2001
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2002
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2003
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2004
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2005
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2006
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2007
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2008
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2009
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2010
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2011
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2012
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2013
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2014
cp ../../ARDAT/ARDAT_temperature.bin ARDAT_temperature_2015
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1991
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1992
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1993
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1994
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1995
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1996
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1997
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1998
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_1999
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2000
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2001
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2002
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2003
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2004
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2005
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2006
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2007
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2008
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2009
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2010
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2011
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2012
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2013
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2014
cp ../../ARDAT/ARDAT_discharge_mps.bin GRACErunoff420x384x12_2015
ln -sf ../../run_template_cube81/*.bin .
rm runoff1p2472-360x180x12.bin
ln -sf ../../run_template_cube81/BATHY_cube81_420x384_arctic .
ln -sf ../../run_template_cube81/WOA05_* .
ln -sf ../../run_template_cube81/*cube81 .
ln -sf ../../run_template_cube81/+ob1979_2014_merge/* .
ln -sf ../../jra25/jra25_*_199[1-9] .
ln -sf ../../jra25/jra25_*_2* .
rm jra25_rain_n*
cp ../buildARDAT/mitgcmuv .
qsub run_arctic_cs_18km_nas

cd ../runGRACE
cp ../../MITgcm_contrib/arctic/cs_18km/input/* .
rm data.obcs_*
cp ../../MITgcm_contrib/arctic/cs_18km/GRACErunoff/data.exf .
cp ../runARDAT/ARDAT_temperature_* .
cp ../runARDAT/GRACErunoff420x384x12_* .
cp ../../GRACE/GRACErunoff420x384x12_20* .
ln -sf ../../run_template_cube81/*.bin .
rm runoff1p2472-360x180x12.bin
ln -sf ../../run_template_cube81/BATHY_cube81_420x384_arctic .
ln -sf ../../run_template_cube81/WOA05_* .
ln -sf ../../run_template_cube81/*cube81 .
ln -sf ../../run_template_cube81/+ob1979_2014_merge/* .
ln -sf ../../jra25/jra25_*_199[1-9] .
ln -sf ../../jra25/jra25_*_2* .
rm jra25_rain_n*
cp ../buildARDAT/mitgcmuv .
qsub run_arctic_cs_18km_nas
