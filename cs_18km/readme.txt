ATN 23-Dec-2013

SEAICE_OPTIONS.h_atn_May2011: Paper: Nguyen et al 2011

CPP_OPTIONS.h_oldcode_before2012: code similar to what is used in Nguyen et al 2011, but compatible with codes around Sep2011
SEAICE_OPTIONS.h_Ian_Sep2011: Equivalence to above but updated to be compatible with Ian's code on that date

After 2011, pkg/seaice has updated significantly with the merge of Ian's code into the main branch.
These followings _phAug2012 files include appropriate flag adjustments + param name change for new code.
They should reproduce almost identical results to Nguyen et al 2011 paper, except for whatever that has
improved in ocean/sea-ice physics.  An offline comparison done by ATN on Aug2012 showed the two runs Aug2012
and Nguyen et al 2011 papers are VERY close to each other, thus ATN considered the two solutions are equivalent.

CPP_OPTIONS.h_phAug2012
SEAICE_OPTIONS.h_phAug2012

To run this set up, use binary inputs from:
http://ecco2.jpl.nasa.gov/data6/arctic/run_template_cube81/
obcs for the years 1979--2014:  +ob1979_2014_merge

There are two bathymetries:
1) original from Nguyen et al 2011: BATHY_cube81_420x384_arctic
2) an updated version to blank out some CAA passages (set bathymetry to zero)
  because 18km can not resolve these passages resutling in excessive sea ice thickness

