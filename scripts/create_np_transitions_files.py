#!/usr/bin/env python3

"""
Produces NP transition matrix files from AIRR tsv.
"""

import pandas as pd
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser(allow_abbrev=False)

parser.add_argument('-i', '--input_file',
                    action="store",
                    required=True)

args = parser.parse_args()

def unroll(data):
    if isinstance(data, dict):
        for key, value in data.items():
            # Recursively unroll the next level and prepend the key to each row.
            for row in unroll(value):
                yield [key] + row
    if isinstance(data, list):
        # This is the bottom of the structure (defines exactly one row).
        yield data

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
      
#Get transition matrices

transition_probs_np1_position = defaultdict(lambda: defaultdict(lambda: [0,0,0,0]))
bases = ["A", "C", "G", "T"]

for seq in np1_seqs:
    for i, base in enumerate(seq[:len(seq) - 1]):
        if base in bases:
            next_base = seq[i + 1]
            if next_base in bases:
                if next_base == "A":
                     transition_probs_np1_position[i][base][0] += 1
                elif next_base == "C":
                    transition_probs_np1_position[i][base][1] += 1
                elif next_base == "G":
                    transition_probs_np1_position[i][base][2] += 1
                elif next_base == "T":
                    transition_probs_np1_position[i][base][3] += 1

for length, dictionary in transition_probs_np1_position.items():
    for key in dictionary:
        count = sum(dictionary[key])
        for i, base in enumerate(dictionary[key]):
            dictionary[key][i] = base/count


transition_probs_np1_position_df = pd.DataFrame(list(unroll(transition_probs_np1_position)))
transition_probs_np1_position_df.columns = ["Length", "Base", "A", "C", "G", "T"]
transition_probs_np1_position_df = transition_probs_np1_position_df.sort_values(by=['Length'])
transition_probs_np1_position_df.to_csv("np1_transition_probs_per_position.csv", index = False)


transition_probs_np2_position = defaultdict(lambda: defaultdict(lambda: [0,0,0,0]))
bases = ["A", "C", "G", "T"]

for seq in np2_seqs:
    for i, base in enumerate(seq[:len(seq) - 1]):
        if base in bases:
            next_base = seq[i + 1]
            if next_base in bases:
                if next_base == "A":
                     transition_probs_np2_position[i][base][0] += 1
                elif next_base == "C":
                    transition_probs_np2_position[i][base][1] += 1
                elif next_base == "G":
                    transition_probs_np2_position[i][base][2] += 1
                elif next_base == "T":
                    transition_probs_np2_position[i][base][3] += 1

for length, dictionary in transition_probs_np2_position.items():
    for key in dictionary:
        count = sum(dictionary[key])
        for i, base in enumerate(dictionary[key]):
            dictionary[key][i] = base/count


transition_probs_np2_position_df = pd.DataFrame(list(unroll(transition_probs_np2_position)))
transition_probs_np2_position_df.columns = ["Length", "Base", "A", "C", "G", "T"]
transition_probs_np2_position_df = transition_probs_np2_position_df.sort_values(by=['Length'])
transition_probs_np2_position_df.to_csv("np2_transition_probs_per_position.csv", index = False)
