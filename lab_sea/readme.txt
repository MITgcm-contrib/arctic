# Instructions for building and running MITgcm-contrib/arctic
# See arctic/doc for additional instructions.

==============
# 1. Get code
  Ensure you have fortran compiler installed: https://github.com/fxcoudert/gfortran-for-macOS/releases
  git clone git@github.com:MITgcm-contrib/arctic
  git clone git@github.com:MITgcm/MITgcm.git
  cd MITgcm
  mkdir build run

==============
# 2. Build executable
  cd build
  ../tools/genmake2 -mo ../../arctic/lab_sea/code
  make depend
  make -j
  note: if using macbook air, use -j8 instead. Otherwise it will try to use more cores than available. 
==============
# 3. Instructions for running simulation
  cd ../run
  ln -sf ../build/mitgcmuv .
  cp ../../arctic/lab_sea/input/* .
  ./mitgcmuv


