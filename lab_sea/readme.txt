# Instructions for building and running MITgcm-contrib/arctic
# See arctic/doc for additional instructions.

Ensure you have fortran compiler installed: https://github.com/fxcoudert/gfortran-for-macOS/releases

My preferreed method:
Install Homebrew installation manager: /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
Follow output instructions to put in path. 
Install gfortran using brew: brew install gcc

==============
# 1. Get code
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
  note: The number here will determine how many processesors to use in next section to run model.

==============
# 3. Instructions for running simulation
  cd ../run
  ln -sf ../build/mitgcmuv .
  cp ../../arctic/lab_sea/input/* .
  ./mitgcmuv

 note: on Wind: mpirun -np 2 mitgcmuv 
