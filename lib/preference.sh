#!/bin/bash


#Run option
SCORING_SCRIPT=scoring_by_RMSD.py

#Environmet

OMP=$(nproc)
OMP_NUM_THREADS=$(nproc)

GMX=/home/software/gmxgpu_serial-2020.6/bin/gmx
GMX_MPI=/home/software/gmxgpu-2020.6/bin/gmx_mpi

