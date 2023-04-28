#!/usr/bin/env python3
"""
Produces VDJ usage files from an AIRR format TSV.
"""

import pandas as pd
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

if set(['v_call','d_call','j_call']).issubset(input_df.columns) == False:
    raise ValueError("One or more of the required columns are missing. v_call, d_call and j_call are required.")

input_df["v_family"] = input_df["v_call"].str.split('-').str[0]
input_df["v_gene"] = input_df["v_call"].str.split('*').str[0]

input_df["d_family"] = input_df["d_call"].str.split('-').str[0]
input_df["d_gene"] = input_df["d_call"].str.split('*').str[0]

input_df["j_family"] = input_df["j_call"].str.split('*').str[0]
input_df["j_gene"] = input_df["j_call"].str.split('*').str[0]

if args.group == None:
    group_column = "group"
    input_df["group"] = "dummy"
else:
    group_column = args.group
    
v_na = input_df.dropna(subset = ["v_family"]).copy()
v_distrib = v_na.groupby(group_column).v_family.value_counts(normalize = True).reset_index(name = "prop").copy()
v_distrib = v_distrib.groupby("v_family").mean().reset_index()
v_distrib["prop"] = v_distrib["prop"]/v_distrib.prop.sum()
v_distrib.to_csv("IGHV_usage.csv", index = False)

v_distrib_gene = v_na.groupby(group_column).v_gene.value_counts(normalize = True).reset_index(name = "prop").copy()
v_distrib_gene = v_distrib_gene.groupby("v_gene").mean().reset_index()
v_distrib_gene["prop"] = v_distrib_gene["prop"]/v_distrib_gene.prop.sum()
v_distrib_gene.to_csv("IGHV_usage_gene.csv", index = False)

d_na = input_df.dropna(subset = ["d_family"]).copy()
d_distrib = d_na.groupby(group_column).d_family.value_counts(normalize = True).reset_index(name = "prop").copy()
d_distrib = d_distrib.groupby("d_family").mean().reset_index()
d_distrib["prop"] = d_distrib["prop"]/d_distrib.prop.sum()
d_distrib.to_csv("IGHD_usage.csv", index = False)

d_distrib_gene = d_na.groupby(group_column).d_gene.value_counts(normalize = True).reset_index(name = "prop").copy()
d_distrib_gene = d_distrib_gene.groupby("d_gene").mean().reset_index()
d_distrib_gene["prop"] = d_distrib_gene["prop"]/d_distrib_gene.prop.sum()
d_distrib_gene.to_csv("IGHD_usage_gene.csv", index = False)

j_na = input_df.dropna(subset = ["j_family"]).copy()
j_distrib = j_na.groupby(group_column).j_family.value_counts(normalize = True).reset_index(name = "prop").copy()
j_distrib = j_distrib.groupby("j_family").mean().reset_index()
j_distrib["prop"] = j_distrib["prop"]/j_distrib.prop.sum()
j_distrib.to_csv("IGHJ_usage.csv", index = False)

j_distrib_gene = j_na.groupby(group_column).j_gene.value_counts(normalize = True).reset_index(name = "prop").copy()
j_distrib_gene = j_distrib_gene.groupby("j_gene").mean().reset_index()
j_distrib_gene["prop"] = j_distrib_gene["prop"]/j_distrib_gene.prop.sum()
j_distrib_gene.to_csv("IGHJ_usage_gene.csv", index = False)

