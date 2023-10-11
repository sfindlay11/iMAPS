#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import RBPamp
import os
import sys
import pandas as pd

user_seq_file = sys.argv[1] # path/to/seq_file.txt
user_motif_dir = sys.argv[2] # path/to/motif_files/

def load_model(fname):
    from RBPamp.params import ModelSetParams
    params = ModelSetParams.load(fname, 1)
    return params

def eval_model(params, seq):
    if len(seq) < params.k:
        return np.zeros(len(seq), dtype=np.float32)

    from RBPamp.cyska import seq_to_bits, PSAM_partition_function
    seqm = seq_to_bits(seq)
    seqm = seqm.reshape((1, len(seq)) )
    accm = np.ones(seqm.shape, dtype=np.float32)
    Z = np.array([PSAM_partition_function(seqm, accm, par.psam_matrix, single_thread=True) * par.A0/params.A0 for par in params])
    return Z.sum(axis=0)[0, :]  # one sequence, sum over all sub-motifs

# read in eCLIP data
data_in = pd.read_csv(user_seq_file, sep = "\t")
seqs = data_in['seq']
rbps = data_in['RBP_peak']

print("loaded seqs data")

# initialize the lists to store results
affs = [None] * len(rbps)
affs_str = [None] * len(rbps)
k_all = [None] * len(rbps)

# start with first RBP (assumes sorted input file)
start_RBP = rbps[0]
params = load_model(user_motif_dir + start_RBP + '.txt')
#print("loaded parameters of width", params.k)

for curr_entry in range(0, len(rbps)):
    # check current RBP for line and reload model if new
    if curr_entry > 0:
        if rbps[curr_entry] != rbps[curr_entry-1]:
            params = load_model(user_motif_dir + rbps[curr_entry] + '.txt')
            print("processing" + rbps[curr_entry])
    # fetch affinities
    affs[curr_entry] = eval_model(params, seqs[curr_entry])
    affs_str[curr_entry] = ','.join(map(str, affs[curr_entry]))
    k_all[curr_entry] = params.k

# add to DF
data_in['aff_string'] = affs_str
data_in['k_size'] = k_all

# write to file
data_in.to_csv(r'rbp_amp_eclip_fetch_general_OUT.txt', index = False, header = True, sep = "\t")