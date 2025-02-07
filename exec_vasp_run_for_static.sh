#!/bin/bash -l
export PATH="/home/arjun/Softwares/vasp.6.2/bin:$PATH"
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/include/hdf5-1.8.12/lib
export OMP_NUM_THREADS=1
folders=$(find . -maxdepth 1 -type d ! -name '.')

for folder in $folders; do
   cd "$folder" || { echo "Failed to enter directory $folder"; continue; }

   if [[ $(basename "$folder") == "vasp_input" || $(basename "$folder") == "figures" ]]; then
      echo "vasp run not expected in $folder as dir. contains just input files for vasp DFPT/static run"
   else
      pwd
      mkdir "static_run"
      cd "static_run" || { echo "Failed to enter directory static_run"; continue; }
      cp "../../vasp_input/INCAR_static" "INCAR"
      cp "../../vasp_input/POTCAR" "POTCAR"
      cp "../../vasp_input/KPOINTS" "KPOINTS"
      cp "../POSCAR" "POSCAR"
      mpirun -np 64 vasp_std
      cd ..
   fi

   cd ..
done

echo "Job done"
