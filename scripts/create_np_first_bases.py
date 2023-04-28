#!/usr/bin/env python3

"""
Produces NP first base usage files from AIRR tsv.
"""

import pandas as pd
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(allow_abbrev=False)

parser.add_argument('-i', '--input_file',
                    action="store",
                    required=True)

args = parser.parse_args()

input_df = pd.read_csv(args.input_file, sep = "\t")

if set(['np1_length','np2_length','sequence',"v_sequence_end", "d_sequence_start","j_sequence_start"]).issubset(input_df.columns) == False:
    raise ValueError("One or more of the required columns are missing. np1_length,np2_length,sequence,v_sequence_end,d_sequence_start and j_sequence_start are required.")

np1 = input_df.copy()
np1 = np1[np1["np1_length"] <= 40] #only take regions up to 40 nuc in length

np1 = np1.dropna(axis=0, subset=['sequence', "v_sequence_end", "d_sequence_start"])
np1['d_sequence_start'] = np1['d_sequence_start'].astype(int)

#Get the np1 sequences 

np1['np1_seq']=[a[b : c - 1 ]for a, b, c, in zip(np1.sequence, np1.v_sequence_end, np1.d_sequence_start)]

np1_seqs=np1.np1_seq.to_list()

np2 = input_df.copy()
np2 = np2[np2["np2_length"] <= 40]

np2 = np2.dropna(subset=['sequence', "d_sequence_end", "j_sequence_start"])
np2['d_sequence_end'] = np2['d_sequence_end'].astype(int) 
np2['j_sequence_start'] = np2['j_sequence_start'].astype(int) 

#Get the np2 sequences

np2['np2_seq']=[a[b : c - 1 ]for a, b, c, in zip(np2.sequence, np2.d_sequence_end, np2.j_sequence_start)]

np2_seqs=np2.np2_seq.to_list()


# Get nucleotide usage at first base of NP1

first_base_probs_np1 = {"A": 0, "C": 0, "G":0, "T": 0}

for seq in np1_seqs:
    if len(seq) >= 1:
        first_base = seq[0]
        if first_base in first_base_probs_np1.keys():
            first_base_probs_np1[first_base] +=1


total = sum(first_base_probs_np1.values())
for k,v in first_base_probs_np1.items():
    first_base_probs_np1[k] = v/total

first_base_probs_np1_df = pd.DataFrame.from_dict(first_base_probs_np1, orient='index')
first_base_probs_np1_df.columns = ["proportion"]
first_base_probs_np1_df.to_csv("np1_first_base_probs.csv")

# Get nucleotide usage at first base of NP2

first_base_probs_np2 = {"A": 0, "C": 0, "G":0, "T": 0}

for seq in np2_seqs:
    if len(seq) >= 1:
        first_base = seq[0]
        if first_base in first_base_probs_np2.keys():
            first_base_probs_np2[first_base] +=1


total = sum(first_base_probs_np2.values())
for k,v in first_base_probs_np2.items():
    first_base_probs_np2[k] = v/total

first_base_probs_np2_df = pd.DataFrame.from_dict(first_base_probs_np2, orient='index')
first_base_probs_np2_df.columns = ["proportion"]
first_base_probs_np2_df.to_csv("np2_first_base_probs.csv")
