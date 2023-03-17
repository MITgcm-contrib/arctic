MITgcm, 36km resolution, Arctic Ocean regional simulation ‘How to’ guide

Dr Jonny Williams, School of Geographical Sciences, Bristol University, UK
jonny.williams@bristol.ac.uk

1.  cd to your MITgcm central directory and check out the code and (some of the) required
    input files
       a. cvs co -d cs_36km_tutorial MITgcm_contrib/arctic/cs_36km/

2.  cd to this new directory and make a ‘build’ directory

3.  Personally I recommend outputting NetCDF files. You will need to do the following if you
    want to do this
       a. add a line with “mnc” on it to cs_36km_tutorial/code/packages.conf
       b. add a line with “useMNC=.TRUE.,” to cs_36km_tutorial/input/data.pkg
       c. add a “data.mnc” file which you can copy from the input directory of any of the
          verification experiments that use NetCDF output

4.  Compile in the usual way (system dependent)

NOTE: the ftp sites listed below are no longer available but the files have
been moved to: https://ecco.jpl.nasa.gov/drive/files/ECCO2/Arctic/arctic_36km

5.  The following files are all required in the run directory
       a. the executable you get when compiling “mitgcmuv”
       b. ftp://ecco2.jpl.nasa.gov/data2/data/runoff/runoff-360x180x12.bin
       c. ftp://ecco2.jpl.nasa.gov/data2/data/atmos/jra25/*1979*
          i. note that the following files may need their “.txt” suffixes removed
             1. jra25_v10m_1979
             2. jra25_u10m_1979
             3. jra25_spfh2m_1979
             4. jra25_temp2m_1979
       d. Also, the run requires the following files which are not included in the
          ftp://ecco2.jpl.nasa.gov/data2/data/atmos/jra25/ directory (I have simply created
          cymbolic links to the 1979 equivalents).
          i. jra25_dlw_1978
         ii. jra25_dsw_1978
        iii. jra25_rain_1978
       e. ftp://ecco2.jpl.nasa.gov/data1/arctic/run_template2/*

6. There are some duplicate file names in
   ftp://ecco2.jpl.nasa.gov/data1/arctic/run_template2/ and cs_36km_tutorial/input and it is
   crucial that those in the latter are used. I have done this by copying them to my run
   directory (platform dependent) second to automatically overwrite the other ones.

7. The model also requires several files which need to be renamed and so I have created
   symbolic links to them.
       a. ln -s OBWv_arctic_210x192.bin OBWv_arctic_210x192_19792014m.bin
       b. ln -s OBWu_arctic_210x192.bin OBWu_arctic_210x192_19792014m.bin
       c. ln -s OBWt_arctic_210x192.stable OBWt_arctic_210x192_19792014m.stable
       d. ln -s OBWt_arctic_210x192.bin OBWt_arctic_210x192_19792014m.bin
       e. ln -s OBWs_arctic_210x192.stable OBWs_arctic_210x192_19792014m.stable
       f. ln -s OBWs_arctic_210x192.bin OBWs_arctic_210x192_19792014m.bin
       g. ln -s OBNv_arctic_210x192.bin OBNv_arctic_210x192_19792014m.bin
       h. ln -s OBNu_arctic_210x192.bin OBNu_arctic_210x192_19792014m.bin
       i. ln -s OBNt_arctic_210x192.stable OBNt_arctic_210x192_19792014m.stable
       j. ln -s OBNt_arctic_210x192.bin OBNt_arctic_210x192_19792014m.bin
       k. ln -s OBNs_arctic_210x192.stable OBNs_arctic_210x192_19792014m.stable
       l. ln -s OBNs_arctic_210x192.bin OBNs_arctic_210x192_19792014m.bin
       m. ln -s OBEv_arctic_210x192.bin OBEv_arctic_210x192_19792014m.bin
       n. ln -s OBEu_arctic_210x192.bin OBEu_arctic_210x192_19792014m.bin
       o. ln -s OBEt_arctic_210x192.stable OBEt_arctic_210x192_19792014m.stable
       p. ln -s OBEt_arctic_210x192.bin OBEt_arctic_210x192_19792014m.bin
       q. ln -s OBEs_arctic_210x192.stable OBEs_arctic_210x192_19792014m.stable
       r. ln -s OBEs_arctic_210x192.bin OBEs_arctic_210x192_19792014m.bin

8. I suggest changing the following in the “data” file just to get you going
       a. endtime=86400.,
          i. i.e. 1 day in seconds
       b. dumpFreq = 2400.,
          i. i.e. outputting data on every timestep, i.e. every 40 minutes in
       seconds

9. Submit the job (platform dependent)

10. Note that the NetCDF output from this simulation will consist of one file for each processor
    (80 by default here, nPx × nPy in the cs_36km_tutorial/code/SIZE.h file).

Acknowledgments
Getting this simulation to run would not have been possible without Fanny Monteiro in the author’s
home institution, members of the MITgcm-support online community and in particular, Dimitris
Menemenlis, Renske Gelderloos and Patrick Heimbach.

=======================

Directory contents for
ftp://ecco2.jpl.nasa.gov/data1/arctic/run_template2/

Matlab script "mk_run_template2.m" was used to generate many of the model
input files described below.

Grid information files needed when OLD_GRID_IO option is defined in
CPP_OPTIONS.h.  See GRID.h for a digram of various distances.
LONC.bin  longitude east of cell center
LATC.bin  latitude north of cell center
LONG.bin  longitude east of southwest corner of cell
LATG.bin  latitude north of southwest corner of cell
DYF.bin   meridional distance in m between V-points
DXF.bin   zonal distance in m between U-points
DYU.bin   meridional distance in m between U-points
DXV.bin   zonal distance in m between V-points
DYC.bin   meridional distance in m between tracer points
DXC.bin   zonal distance in m between tracer points
DYG.bin   meridional distance in m between cell corners
DXG.bin   zonal distance in m between cell corners
RAZ.bin   vertical face area in m^2 for vorticity points
RAW.bin   vertical face area in m^2 for u cells
RAS.bin   vertical face area in m^2 for v cells
RA.bin    vertical face area in m^2 for tracer cells

Model bathymetry, initial, and surface boundary condition files.
ETOPO2_210x192_arctic        model bathymetry (m)
WGHC_S_210x192x50_arctic     initial salinity (g/kg)
WGHC_T_210x192x50_arctic     initial potential temperature (deg C)
AREA_210x192_arctic.cube81   initial ice concentration (fractional >=0, <=1)
HEFF_210x192_arctic.cube81   initial effective sea ice thickness (m)
HSALT_210x192_arctic.cube81  initial effective sea ice salinity (g/m^2)
HSNOW_210x192_arctic.cube81  initial effective snow thickness (m)

Open boundary condition files.
OBNs_arctic_210x192.stable  North open boundary conditions, salinity (g/kg)
OBNt_arctic_210x192.stable  North open boundary conditions, temperature (deg C)
OBNu_arctic_210x192.bin     North open boundary conditions, U-velocity (m/s)
OBNv_arctic_210x192.bin     North open boundary conditions, V-velocity (m/s)
OBEs_arctic_210x192.stable  South open boundary conditions, salinity (g/kg)
OBEt_arctic_210x192.stable  South open boundary conditions, temperature (deg C)
OBEu_arctic_210x192.bin     South open boundary conditions, U-velocity (m/s)
OBEv_arctic_210x192.bin     South open boundary conditions, V-velocity (m/s)
OBWs_arctic_210x192.stable  West open boundary conditions, salinity (g/kg)
OBWt_arctic_210x192.stable  West open boundary conditions, temperature (deg C)
OBWu_arctic_210x192.bin     West open boundary conditions, U-velocity (m/s)
OBWv_arctic_210x192.bin     West open boundary conditions, V-velocity (m/s)

Runtime parameter files.
eedata
data
data.pkg
data.cal
data.exf
data.gmredi
data.seaice
data.obcs
data.kpp
data.salt_plume
data.diagnostics
