#!/bin/bash -l
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
      #mkdir "static_run"
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
