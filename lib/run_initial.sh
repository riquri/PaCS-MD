#!/bin/bash

if [ -z "$1" ]
then
   echo :Usage ./$( basename $0 ) name
   exit
fi

name=$1
RUNID=$3
WDIR=$(pwd)

source $WDIR/lib/preference.sh

########################################
# Make input files
########################################
cd $WDIR/input
$WDIR/lib/run_leap.sh

# Check
if [ -e ./initial.gro ]; then
    echo "Input structure is successfully built."
else
    echo "ERROR on Input file making."
    exit 1
fi



########################################
# Energy minimization
########################################
echo "########################################"
echo "Start energy minimization."
echo "########################################"

##########
# Define .mdp file
##########
echo "########################################"
echo "Make .mdp file for energy minimization for all atoms."
echo "########################################"
cat << eof0 > ./min.mdp
; minim.mdp - used as input into grompp to generate em.tpr
integrator	= steep		; Algorithm (steep = steepest descent minimization)
emtol		= 1000.0  	; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01      ; Energy step size
nsteps		= 10000	  	; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist		    = 1		    ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet
ns_type		    = grid		; Method to determine neighbor list (simple, grid)
coulombtype	    = PME		; Treatment of long range electrostatic interactions
rcoulomb	    = 1.0		; Short-range electrostatic cut-off
rvdw		    = 1.0		; Short-range Van der Waals cut-off
pbc		        = xyz 		; Periodic Boundary Conditions (yes/no)
eof0


##########
# Define .mdp file
##########
echo "########################################"
echo "Make .mdp file for NVT Equilibration."
echo "########################################"


cat << eof1 > nvt.mdp
title		= $name NVT equilibration 
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 0		; save coordinates every 1.0 ps
nstvout		= 0		; save velocities every 1.0 ps
nstenergy	= 0		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
nstxout-compressed  = 500
; Bond parameters
continuation	        = no		; first dynamics run
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = h-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10		; 20 fs, largely irrelevant with Verlet
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		; cubic interpolation
fourierspacing	= 0.16	; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1           ; time constant, in ps
ref_t		= 300 	  300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl		= no 		; no pressure coupling in NVT
; Periodic boundary conditions
pbc		= xyz		    ; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= yes		; assign velocities from Maxwell distribution
gen_temp	= 300		; temperature for Maxwell distribution
gen_seed	= $RANDOM		; generate a random seed

eof1

########################################
# NPT Equilibrations
########################################
echo "########################################"
echo "Start NPT Equilibrations."
echo "########################################"

##########
# Define .mdp file
##########
echo "########################################"
echo "Make .mdp file for NPT Equilibrations."
echo "########################################"

cat << eof2 > ./npt.mdp
title		= $name NPT equilibration
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000		; 2 * 50000 = 100 ps
dt		    = 0.002		; 2 fs
; Output control
nstxout		= 500		; save coordinates every 1.0 ps
nstvout		= 500		; save velocities every 1.0 ps
nstenergy	= 500		; save energies every 1.0 ps
nstlog		= 500		; update log file every 1.0 ps
; Bond parameters
continuation	        = yes		; Restarting after NVT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = h-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Berendsen	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
refcoord_scaling    = com
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel		= no		; Velocity generation is off 

eof2



########################################
# Production MD
########################################
echo "########################################"
echo "Start Production MD."
echo "########################################"


##########
# Define .mdp file
##########
echo "Make .mdp file for Production MD."
cat << eof1 > ./sample.mdp
title		= $name MD simulation (NPT)
; Run parameters
integrator	= md		; leap-frog integrator
nsteps		= 50000	; 2 * 50000 = 100 ps (100 ps)
dt		    = 0.002		; 2 fs
; Output control
nstxout		        = 0		; save coordinates every 10.0 ps
nstvout		        = 0		; save velocities every 10.0 ps
nstenergy	        = 0		; save energies every 1.0 ps
nstlog		        = 500		; update log file every 10.0 ps
nstxout-compressed  = 500       ; save compressed coordinates every 10.0 ps
                                ; nstxout-compressed replaces nstxtcout
compressed-x-grps   = System    ; replaces xtc-grps
; Bond parameters
;continuation	        = yes		; Restarting after NPT 
constraint_algorithm    = lincs	    ; holonomic constraints 
constraints	            = h-bonds	; all bonds (even heavy atom-H bonds) constrained
lincs_iter	            = 1		    ; accuracy of LINCS
lincs_order	            = 4		    ; also related to accuracy
; Neighborsearching
cutoff-scheme   = Verlet
ns_type		    = grid		; search neighboring grid cells
nstlist		    = 10	    ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb	    = 1.0		; short-range electrostatic cutoff (in nm)
rvdw		    = 1.0		; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype	    = PME		; Particle Mesh Ewald for long-range electrostatics
pme_order	    = 4		    ; cubic interpolation
fourierspacing	= 0.16		; grid spacing for FFT
; Temperature coupling is on
tcoupl		= V-rescale	            ; modified Berendsen thermostat
tc-grps		= Protein Non-Protein	; two coupling groups - more accurate
tau_t		= 0.1	  0.1	        ; time constant, in ps
ref_t		= 300 	  300	        ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl		        = Parrinello-Rahman	    ; Pressure coupling on in NPT
pcoupltype	        = isotropic	            ; uniform scaling of box vectors
tau_p		        = 2.0		            ; time constant, in ps
ref_p		        = 1.0		            ; reference pressure, in bar
compressibility     = 4.5e-5	            ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc		= xyz		; 3-D PBC
; Dispersion correction
DispCorr	= EnerPres	; account for cut-off vdW scheme
; Velocity generation
gen_vel     = yes       ; assign velocities from Maxwell distribution
gen_temp    = 300       ; temperature for Maxwell distribution
gen_seed    = $RANDOM        ; generate a random seed

eof1

##########
# 
##########

cd $WDIR/input

name=initial

$GMX make_ndx -f ${name}.gro << eof
q
eof

$GMX grompp -f ./min.mdp -c ./${name}.gro -p ./${name}.top -o ./em.tpr -n ./index.ndx
$GMX_MPI mdrun -v -deffnm ./em

$GMX grompp -maxwarn 1 -f nvt.mdp -c ./em.gro -p ./${name}.top -o nvt.tpr -r ./em.gro
$GMX_MPI mdrun -deffnm nvt -ntomp ${OMP} -v -cpo nvt.cpt

$GMX grompp -maxwarn 1 -f npt.mdp -c nvt.gro -t nvt.cpt -p ${name}.top -o npt.tpr -r nvt.gro
$GMX_MPI mdrun -deffnm npt -ntomp ${OMP} -v -cpo npt.cpt

$GMX grompp -maxwarn -1 -f sample.mdp -c ./npt.gro -t ./npt.cpt -p ./${name}.top -o sample.tpr -r ./npt.gro -n index.ndx
$GMX_MPI mdrun -deffnm sample -ntomp ${OMP}

$GMX trjconv -f em.gro -s em.gro -o em_protein.gro -n ./index.ndx<<EOF
1
EOF

$GMX trjconv -f sample.xtc -s em.gro -o sample_protein.xtc -n ./index.ndx<<EOF
1
EOF

echo "All candicates are ready."

##########
# Check
##########
if [ -e ./sample.xtc ]; then
    echo "Production MD has done."
else
    echo "ERROR on Production MD"
    exit 1
fi


cd $WDIR
mkdir ./cyc0
cd ./cyc0
mkdir ./candi1
cp $WDIR/input/sample.xtc ./candi1
cp $WDIR/input/sample_protein.xtc ./candi1
cd candi1


#Calculate score
python3 $WDIR/lib/$SCORING_SCRIPT

#Ranking
cd $WDIR/cyc0
python3 $WDIR/lib/merge_score.py 0 10

echo "Finished the initial cycle."

exit 0
