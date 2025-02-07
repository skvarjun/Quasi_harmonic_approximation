#!/bin/bash -l
conda activate

folders=$(find . -maxdepth 1 -type d ! -name '.')
for folder in $folders; do
    cd "$folder" || { echo "Failed to enter directory $folder"; continue; }
    if [[ $(basename "$folder") = "vasp_input" ]]; then
      echo "thermal_properties.yaml not expecting in $folder"
    else
      mv POSCAR POSCAR_
      phonopy --fc vasprun.xml
      phonopy-load --mesh 31 31 31 -t
      mv POSCAR_ POSCAR
    fi
    cd ..
done

echo "Job done"

