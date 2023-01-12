# PaCS-MD

The shell scripts to performe PaCS-MD
  
# DEMO  


# Features

 
# Requirement
- Gromacs (2018 or later)
- AmberTools
- Python 3.x
 - MDAnalysis
 - numpy
 
# Installation
 Copy the lib directory on the your each working directory.
 ```
 mkdir protein_X
 cd ./protein_X
 git clone https://github.com/riquri/PaCS-MD.git
 ```
 
 Place the initial and target structure
 ```
 cp {your_initial_structure.pdb} ./input/initial.pdb
 cp {your_target_structure.pdb} ./input/initial.pdb
 ```
 
 Modify the preference file (./lib/preference.sh)
 
 The `GMX` and `GMX_MPI` is set for the full path for the Gromacs binary.
 
 You can specify scoring script as `SCORING_SCRIPT` in ./lib directory.
 If you need your original scoring script, place it in the ./lib directory.
 
# Usage
 ```
 nohup ./lib/run_jobs.sh {TRIAL_NAME} &
 ```
 
 To restart PaCS-MD,
 ```
 cd {TRIAL_directory}
 nohup ./lib/run_pacs {START_CYCLE} {LAST_CYCLE} {NUMBER_OF_CANDIDATE}
 ```

 
# Note
Please cite this paper.


Ryuhei Harada and Akio Kitao. Parallel cascade selection molecular dynamics (PaCS-MD) to generate conformational transition pathway. J. Chem. Phys. (2021)DOI:[https://aip.scitation.org/doi/10.1063/1.4813023]

 
# Contributers
- Rikuri Morita*, Ryuhei Harada*.
- Center for Computational Sciences, University of Tsukuba
- morita@ccs.tsukuba.ac.jp
 
# License
 
This scripts are under MIT license.
