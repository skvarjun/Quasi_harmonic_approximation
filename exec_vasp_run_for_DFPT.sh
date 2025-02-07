#!/bin/bash -l
#export PATH="/home/arjun/Softwares/vasp.6.2/bin:$PATH"
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/include/hdf5-1.8.12/lib
ml intel-oneapi-compilers/2022.2.1  intel-oneapi-mpi/2021.1.1
ml  vasp/6.3.2
export OMP_NUM_THREADS=1
folders=$(find . -maxdepth 1 -type d ! -name '.')

for folder in $folders; do
   cd "$folder" || { echo "Failed to enter directory $folder"; continue; }

   if [[ $(basename "$folder") == "vasp_input" || $(basename "$folder") == "figures" ]]; then
      echo "vasp run not expected in $folder as dir. contains just input files for vasp DFPT/static run"
   else
      pwd
      mv POSCAR POSCAR-unitcell
      phonopy -d --dim 2 2 1 -c POSCAR-unitcell
      phonopy -d --dim 1 1 1 --pa auto -c SPOSCAR
      mv SPOSCAR POSCAR
      mpirun -np 64 vasp_std
   fi
   cd ..
done

echo "Job done"