# PaCS-MD

The scripts to performe PaCS-MD
  
# DEMO  


# Features
The current implementation is for local workstation.
The parallel version for super computer is under construction.
 
# Requirement
- Gromacs (2018 or later)
- AmberTools
- Python 3.x
 - MDAnalysis
 - numpy
 - matplotlib (for analysis)
 
 ```
 conda install -c conda-forge mdanalysis
 ```
 
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
 cp {your_target_structure.pdb} ./input/target.pdb
 ```
 
 Modify the leap input file (./lib/run_leap.sh)
 ```
 set mol box {50 50 50 } # Size of the PBC box
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
 nohup ./lib/run_pacs.sh {START_CYCLE(not include the exist cycle)} {LAST_CYCLE} {NUMBER_OF_CANDIDATE=10}
 ```

 To visualize the progress with score,
 ```
 cd {TRIAL_directory}
 python3 ./lib/plot_score.py
 ```
 This script draw a plot and save a pdf file.
 
# Note
Please cite this paper.


Ryuhei Harada and Akio Kitao. Parallel cascade selection molecular dynamics (PaCS-MD) to generate conformational transition pathway. J. Chem. Phys. (2021)DOI:[https://doi.org/10.1063/1.4813023]

 
# Contributers
- Rikuri Morita*, Ryuhei Harada*.
- Center for Computational Sciences, University of Tsukuba
- morita@ccs.tsukuba.ac.jp, ryuhei@ccs.tsukuba.ac.jp
 
# License
 
This scripts are under MIT license.
