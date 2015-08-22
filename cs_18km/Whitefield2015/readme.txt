Instructions for running optimized (A1) simulation from Nguyen et al (2011)
http://ecco2.org/manuscripts/2011/NguyenJGR2011.pdf
but using river runoff from Whitefield et al (2015)
http://ecco2.org/manuscripts/2015/Whitefield2015.pdf

===============

DM 10-Aug-20215
Instructions for checkpoint65m (2015/06/15) on pleiades
with Whitefield et al (2015) river runoff discharge (runARD)
and river runoff discharge and temperature (runARDAT)

-> get a copy of ftp://ecco2.jpl.nasa.gov/data6/arctic/run_template_cube81
   also available on pfe:/nobackup/dmenemen/arctic18km/run_template_cube81
   and harbor.mit.edu:/net/nares/raid8/ecco-shared/arctic18km/run_template_cube81

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/atmos/jra25
   also available on pfe:/nobackup/hzhang1/forcing/jra25

-> get a copy of ftp://ecco2.jpl.nasa.gov/data2/data/runoff/ARDAT
   also available on pfe:/nobackup/dmenemen/forcing/runoff/ARDAT

cvs co MITgcm_contrib/arctic/cs_18km
cvs co -r checkpoint65m MITgcm_code
cd MITgcm
mkdir buildARDAT runARD runARDAT
cd buildARDAT
module purge
module load comp-intel/2015.0.090 mpi-sgi/mpt.2.11r13 netcdf/4.0
cp ../../MITgcm_contrib/arctic/cs_18km/code/* .
cp ../../MITgcm_contrib/arctic/cs_18km/Whitefield2015/EXF_OPTIONS.h .
../tools/genmake2 -of ../tools/build_options/linux_amd64_ifort+mpi_ice_nas
make depend
make -j 16

cd ../runARD
cp ../../MITgcm_contrib/arctic/cs_18km/input/* .
rm data.obcs_*
cp ../../MITgcm_contrib/arctic/cs_18km/Whitefield2015/data.exf_ARD data.exf
ln -sf ../../ARDAT/ARDAT_discharge_mps.bin .
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

cd ../runARDAT
cp ../../MITgcm_contrib/arctic/cs_18km/input/* .
rm data.obcs_*
cp ../../MITgcm_contrib/arctic/cs_18km/Whitefield2015/data.exf_ARDAT data.exf
ln -sf ../../ARDAT/ARDAT_discharge_mps.bin .
ln -sf ../../ARDAT/ARDAT_temperature.bin .
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
