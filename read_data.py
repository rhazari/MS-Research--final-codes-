# -*- coding: utf-8 -*-
"""
Created on Fri Aug 16 00:23:30 2013

@author: ryan
"""

import numpy as np

file_name = []
for i in range(1,101):
    file_name.append("ECDS_Fwd_"+str(i)+".txt")

count = 25   
avg_MaxLF = np.zeros(count)
num_files = len(file_name)
for m in range(num_files):
    load_data = np.loadtxt(file_name[m])
    avg_MaxLF += load_data[0:count]

avg_MaxLF = avg_MaxLF/num_files