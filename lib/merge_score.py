import os
import glob
import re
import csv
import sys
import copy

import numpy as np
import math

max_rank = int(sys.argv[2])

score_file_list = sorted(glob.glob(f"./candi*/score.txt"), key=lambda x: int(x.split("/")[1].replace("candi", "")))

score_list = []
for candi_name in score_file_list:
    candi_index = int(candi_name.split("/")[1].replace("candi", ""))
    candi_score = np.loadtxt(candi_name, delimiter=",")
    score_list.extend(np.hstack([np.full((np.shape(candi_score)[0],1), candi_index) , candi_score]))

#score_list = candi, frame, score 
#pre_rank pre_cycle pre_candi pre_frame pre_score
output = ""
for rank, c in enumerate(sorted(score_list, key=lambda x: x[2], reverse=False)):
    if rank >= max_rank:
        break
    output += f"{rank+1},{sys.argv[1]},{int(c[0])},{int(c[1])},{c[2]}\n"

with open("rank.csv", "w") as f:
    f.write(output)


