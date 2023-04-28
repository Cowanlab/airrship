#!/usr/bin/env python3

"""
Produces SHM reference files from AIRR tsv.
"""

import pandas as pd
import argparse
from collections import defaultdict

parser = argparse.ArgumentParser(allow_abbrev=False)

parser.add_argument('-i', '--input_file',
                    action="store",
                    required=True)

args = parser.parse_args()

input_df = pd.read_csv(args.input_file, sep = "\t")

if set(['v_call','germline_alignment','sequence_alignment',"fwr1", "fwr2","fwr3","cdr1", "cdr2", "cdr3"]).issubset(input_df.columns) == False:
    raise ValueError("One or more of the required columns are missing. v_call, germline_alignment, sequence_alignment, fwr1, fwr2, fwr3, cdr1, cdr2, cdr3 are required.")

def default_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def get_mut_freq(sequence, germline):   
    if len(sequence) != len(germline):
        return(np.nan)
    mut_count = 0
    length = 0
    for i, base in enumerate(sequence):
        if base != "." and base != "-" and base != "N" and germline[i] != "N" and germline[i] != "-" and germline[i] !=  "." :
            length += 1
            if base != germline[i]:
                mut_count += 1
    return(mut_count / length)


input_df = input_df.dropna(subset = ["germline_alignment", "sequence_alignment"])
input_df["v_family"] = input_df["v_call"].str.split("-").str[0]


fwr1_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
fwr2_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
fwr3_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
fwr4_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
cdr1_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
cdr2_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
cdr3_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
fwr_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})
cdr_kmer_base_usage = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})


for row in input_df.itertuples():
    
    seq = row.sequence_alignment
    gaps = 0
    for x in seq:
        if x == ".":
            gaps += 1
        else:
            break

    ungapped_seq = row.sequence_alignment.replace(".", "")
    germline = row.germline_alignment
    germline = germline[gaps: ]
    ungapped_germline = germline.replace(".", "")
    
    fwr1 = row.fwr1
    cdr1 = row.cdr1
    fwr2 = row.fwr2
    cdr2 = row.cdr2
    fwr3 = row.fwr3
    cdr3 = row.cdr3
    
    fwr1_len = len(fwr1.replace(".", ""))
    cdr1_len = len(cdr1.replace(".", ""))
    fwr2_len = len(fwr2.replace(".", ""))
    cdr2_len = len(cdr2.replace(".", ""))
    fwr3_len = len(fwr3.replace(".", ""))
    cdr3_len = len(cdr3.replace(".", ""))
    # fwr4_len = len(fwr4.replace(".", ""))

    cdr1_start = fwr1_len
    cdr1_end = fwr1_len + cdr1_len
    fwr2_end = cdr1_end + fwr2_len
    cdr2_end = fwr2_end + cdr2_len
    fwr3_end = cdr2_end + fwr3_len
    cdr3_end = fwr3_end + cdr3_len
         
    
    if len(ungapped_seq) != len(ungapped_germline):
        continue
       
    for x in range(0, len(ungapped_germline) - 4):
        seq_kmer = ungapped_seq[x:x+5]
        if "-" in seq_kmer or "N" in seq_kmer or "." in seq_kmer:
            continue
        germ_kmer = ungapped_germline[x:x+5]
        if "-" in germ_kmer or "N" in germ_kmer or "." in germ_kmer:
            continue
        base = seq_kmer[2]
        if cdr1_start <= x + 2 < cdr1_end:
            cdr1_kmer_base_usage[germ_kmer][base] += 1
            cdr_kmer_base_usage[germ_kmer][base] += 1
        elif fwr2_end <= x+2 < cdr2_end:
            cdr2_kmer_base_usage[germ_kmer][base] += 1
            cdr_kmer_base_usage[germ_kmer][base] += 1
        elif fwr3_end <= x+2 < cdr3_end:
            cdr3_kmer_base_usage[germ_kmer][base] += 1
            cdr_kmer_base_usage[germ_kmer][base] += 1       
        elif x + 2 < cdr1_start:
            fwr1_kmer_base_usage[germ_kmer][base] += 1
            fwr_kmer_base_usage[germ_kmer][base] += 1
        elif cdr1_end <= x + 2 < fwr2_end:
            fwr2_kmer_base_usage[germ_kmer][base] += 1
            fwr_kmer_base_usage[germ_kmer][base] += 1
        elif cdr2_end <= x + 2 < fwr3_end:
            fwr3_kmer_base_usage[germ_kmer][base] += 1
            fwr_kmer_base_usage[germ_kmer][base] += 1
        else:
            fwr4_kmer_base_usage[germ_kmer][base] += 1
            fwr_kmer_base_usage[germ_kmer][base] += 1


                 

for kmer_dict,region in zip([cdr1_kmer_base_usage, cdr2_kmer_base_usage, cdr3_kmer_base_usage, fwr1_kmer_base_usage, fwr2_kmer_base_usage, fwr3_kmer_base_usage, fwr4_kmer_base_usage, cdr_kmer_base_usage, fwr_kmer_base_usage], ["cdr1","cdr2","cdr3","fwr1", "fwr2","fwr3","fwr4","cdr","fwr"]):
    for kmer in kmer_dict.values():
        total = sum(kmer.values())
        for k,v in kmer.items():
            kmer[k] = v/total

    kmer_df = pd.DataFrame.from_dict(kmer_dict, orient = "index")
    kmer_df["kmer"] = kmer_df.index
    kmer_df = kmer_df[["kmer", "A", "C", "G", "T"]]
    kmer_df.to_csv(f"{region}_kmer_base_usage.csv", index = False)
    
    
input_df['mut_freq'] = input_df.apply(lambda row : get_mut_freq(row['sequence_alignment'],
                 row['germline_alignment']), axis = 1)

mut_per_fam_df = input_df.groupby("v_family")["mut_freq"].value_counts(normalize = True).reset_index(name = "proportion").copy()

mut_per_fam_df.to_csv("mut_freq_per_seq_per_family.csv", index = False)
