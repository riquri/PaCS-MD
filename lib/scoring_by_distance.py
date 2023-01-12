import os
import glob
import re
import csv
import sys
import copy

import numpy as np
import math

import MDAnalysis as mda
import MDAnalysis.transformations
import MDAnalysis.analysis.rms

import itertools


traj_name = f"./sample.xtc"

scores = []

ref = mda.Universe("../../input/target.pdb")
u = mda.Universe(traj_name.replace("sample.xtc", "../../input/em.gro"), traj_name)

u = mda.Universe(traj_name.replace("sample_protein.xtc", "em_protein.gro"), traj_name)
c_alpha = u.select_atoms("protein and name CA")
transform = mda.transformations.fit_rot_trans(c_alpha, ref_c_alpha)
u.trajectory.add_transformations(transform)

sele1 = u.select_atoms("resid 3 and name N")
sele2 = u.select_atoms("resid 8 and name O")

for ts in u.trajectory:
    score.append([ts.frame, np.linarg.norm(sele1.positions[0] - sele2.positions[1])])

np.savetxt("score.txt", score, fmt ="%.4f", delimiter=",")


