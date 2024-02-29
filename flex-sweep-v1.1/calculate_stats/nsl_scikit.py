#!/bin/python3.6

import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import allel
import numpy as np
import sys

# load in data
data = np.genfromtxt(sys.argv[1],skip_header=6,dtype='str')
#data = np.genfromtxt(sys.argv[1],skip_header=7,dtype='str') # if the header is longer, check that there aren't echo statements in run_msms and/or run_discoal
map = np.genfromtxt(sys.argv[2])
out = sys.argv[3]

# adjust
iter = sys.argv[1].split("_")[-1].split(".")[0]
split = [list(s) for s in data]
print(np.transpose(split).shape)
h = allel.HaplotypeArray(np.transpose(split),dtype='i1')
h_temp = h[int(map[:,1][0]):,:]
h_rest = h_temp[:map.shape[0],:] # to match the map file, which cuts off the first and last 100000 bases
daf = (h_rest.count_alleles()/h_rest.count_called(axis=1)[0])[:,1]
h_maf = h_rest[daf>=0.05,:] # because selscan only uses these
ids=map[:,2]
h_maf_ids = ids[daf>=0.05]

# calculate
nsl = allel.nsl(h_maf)
#print(nsl)
with open(out,"a") as f:
        for i,j,k in zip(h_maf_ids.astype(int),daf[daf>=0.05],nsl):
                f.write("{} {} {}\n".format(i,j,k))
