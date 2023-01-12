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

R = MDAnalysis.analysis.rms.RMSD(u, ref, select="protein and name CA", groupselections=["protein and name CA"])
R.run()
np.savetxt("score.txt", R.results.rmsd[:,1:3], fmt ="%.4f", delimiter=",")


