import os
import glob
import re
import csv
import sys
import copy

import numpy as np
import math

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.size'] = 8
matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['font.family'] = "Helvetica"
matplotlib.style.use("tableau-colorblind10")

matplotlib.rcParams.update({'figure.figsize': [83/24.5,83/24.5]})



file_list = sorted(glob.glob(f"./cyc*/rank.csv"), key=lambda x: int(x.split("/")[1].replace("cyc", "")))

rank_list = []
for cyc_name in file_list:
    cyc_index = int(cyc_name.split("/")[1].replace("cyc", ""))
    cyc_rank = np.loadtxt(cyc_name, delimiter=",")
    rank_list.append(cyc_rank[0][-1])
print(rank_list)

plt.figure(figsize=(83/24.5,83/24.5))
plt.plot(range(1, len(rank_list)+1), rank_list)
plt.xlim(0,len(rank_list)+1)
plt.xlabel("Cycle")
plt.ylabel("Score")
plt.tight_layout()
plt.savefig("cycle_score.pdf")




