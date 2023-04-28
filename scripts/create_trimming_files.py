#!/usr/bin/env python3

"""
Produces gene end trimming files from AIRR tsv.
"""

import pandas as pd
import argparse

parser = argparse.ArgumentParser(allow_abbrev=False)

parser.add_argument('-i', '--input_file',
                    action="store",
                    required=True)

parser.add_argument('--imgt_v',
                    action="store",
                    default="imgt_human_IGHV.fasta",
                    help = "Path to IMGT V allele file." )


parser.add_argument('--imgt_d',
                    action="store",
                    default="imgt_human_IGHD.fasta",
                    help = "Path to IMGT D allele file." )

args = parser.parse_args()

input_df = pd.read_csv(args.input_file, sep = "\t")

v_fasta = args.imgt_v
d_fasta = args.imgt_d

if set(['v_call','d_call','v_germline_end',"j_germline_start", "d_germline_start","d_sequence_end","d_sequence_start"]).issubset(input_df.columns) == False:
    raise ValueError("One or more of the required columns are missing. v_call,d_call,v_germline_end,j_germline_start,d_germline_start,d_sequence_end and d_sequence_start are required.")

input_df["v_family"] = input_df["v_call"].str.split('-').str[0]

input_df["d_family"] = input_df["d_call"].str.split('-').str[0]

input_df["j_family"] = input_df["j_call"].str.split('*').str[0]

v_lengths_gapped=dict()

with open(f'{v_fasta}') as germline_file: #imgt gapped IGHV reference file
    allele='' #make empty val
    for line in germline_file:
        if line.startswith('>'):
            allele = line.rstrip('\n').split('|')[1]
            length = line.rstrip('\n').split("|", 13)[12].split("=")[1]
            v_lengths_gapped[allele]=length 
            
v_lengths_gapped_df = pd.DataFrame.from_dict(v_lengths_gapped, orient = "index").reset_index()
v_lengths_gapped_df.columns = ["v_call", "v_call_germ_length_gapped"]

d_lengths=dict()


with open(f'{d_fasta}') as germline_file: 
    allele='' #make empty val
    for line in germline_file:
        if line.startswith('>'):
            allele=line.rstrip('\n').split('|')[1]
            length=line.rstrip('\n').split("|", 7)[6:7][0].split(" ")[0]
            d_lengths[allele]=length 

d_lengths_df = pd.DataFrame.from_dict(d_lengths, orient = "index").reset_index()
d_lengths_df.columns = ["d_call", "d_call_germ_length"]

input_df = pd.merge(input_df, v_lengths_gapped_df, on="v_call", how = "left")
input_df = pd.merge(input_df, d_lengths_df, on="d_call", how = "left")

v_trim = input_df.dropna(subset=["v_call_germ_length_gapped", "v_germline_end"]).copy()
d_trim = input_df.dropna(subset=["j_germline_start", "d_call_germ_length", "d_sequence_end", "d_sequence_start"]).copy()
j_trim = input_df.dropna(subset=["d_germline_start"]).copy()

# Manually work out trimming lengths

v_trim["v_3p_del"] = v_trim["v_call_germ_length_gapped"].astype(int) - v_trim["v_germline_end"].astype(int) 
j_trim["j_5p_del"] = j_trim["j_germline_start"] - 1
d_trim["d_5p_del"] = d_trim["d_germline_start"] - 1
d_trim.d_5p_del = d_trim.d_5p_del.astype(int)
d_trim["d_3p_del"] = d_trim["d_call_germ_length"].astype(int) - (d_trim["d_sequence_end"].astype(int) - d_trim["d_sequence_start"].astype(int) + 1)

grouped_v = v_trim.query("v_3p_del >= 0").groupby(['v_family', 'v_3p_del']).size()
prop_v_df = grouped_v.groupby(level=0).apply(lambda x: x / float(x.sum())).reset_index(name='proportions')
prop_v_df.to_csv("V_family_trimming_proportions.csv", index = False)

grouped_d_3 = d_trim.query("d_3p_del >= 0").groupby(['d_family', 'd_3p_del']).size()
prop_d_3_df = grouped_d_3.groupby(level=0).apply(lambda x: x / float(x.sum())).reset_index(name='proportions')
prop_d_3_df.to_csv("D_3_family_trimming_proportions.csv", index = False)

grouped_d_5 = d_trim.groupby(['d_family', 'd_5p_del']).size()
prop_d_5_df = grouped_d_5.groupby(level=0).apply(lambda x: x / float(x.sum())).reset_index(name='proportions')
prop_d_5_df.to_csv("D_5_family_trimming_proportions.csv", index = False)

grouped_j = j_trim.groupby(['j_family', 'j_5p_del']).size()
prop_j_df = grouped_j.groupby(level=0).apply(lambda x: x / float(x.sum())).reset_index(name='proportions')
prop_j_df.to_csv("J_family_trimming_proportions.csv", index = False)
