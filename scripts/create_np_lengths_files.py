#!/usr/bin/env python3

"""
Produces NP length files from AIRR format TSV.
"""

import pandas as pd
from collections import defaultdict
import argparse


parser = argparse.ArgumentParser(allow_abbrev=False)

parser.add_argument('-i', '--input_file',
                    action="store",
                    required=True)

parser.add_argument('--group',
                   action="store",
                   help="Name of column containing dataset metadata to group by.",
                   default = None)

args = parser.parse_args()

input_df = pd.read_csv(args.input_file, sep = "\t")

if args.group == None:
    group_column = "group"
    input_df["group"] = "dummy"
else:
    group_column = args.group

if set(['np1_length','np2_length']).issubset(input_df.columns) == False:
    raise ValueError("One or more of the required columns are missing. np1_length, np2_length are required.")


np1 = input_df.copy()
np1 = np1[np1["np1_length"] <= 40] #only take regions up to 40 nuc in length

np1 = np1.dropna(subset = ["np1_length"])

#Get the np1 sequences 

np2 = input_df.copy()
np2 = np2[np2["np2_length"] <= 40]

np2 = np2.dropna(subset = ["np2_length"])

#Get the np lengths

group_np1 = np1.groupby(group_column).np1_length.value_counts(normalize = True).reset_index(name = "prop")
group_np1 = group_np1.groupby("np1_length").mean().reset_index()
group_np1.to_csv("np1_lengths_proportions.csv", index = False)


group_np2 = np2.groupby(group_column).np2_length.value_counts(normalize = True).reset_index(name = "prop")
group_np2 = group_np2.groupby("np2_length").mean().reset_index()
group_np2.to_csv("np2_lengths_proportions.csv", index = False)
