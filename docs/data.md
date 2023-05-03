# Input files

AIRRSHIP uses various input data in order to accurately replicate real repertoires. 
If desired, an alternative data directory can be specified. Files in this directory must follow the format of the included data, including exactly matching headers and column order, and be named in the same manner.

## Scripts for producing input files

A selection of scripts are in included in the scripts folder on Github and can be used to process AIRR TSV files into the required formats for AIRRSHIP input.

Each reference file requires that certain AIRR-C format columns are present - to produce all reference files, the following columns must be included in your input: sequence,sequence_alignment, germline_alignment, fwr1, fwr2, fwr3, cdr1, cdr2, cdr3, v_call, d_call, j_call, np1_length, np2_length, v_sequence_end, d_sequence_start, j_sequence_start. You must also provide paths to IMGT V and D allele FASTA files when calculating trimming metrics. Help for each script can be accessed using script_name.py -h.

These scripts apply very limited pre-processing to the provided sequences to allow for flexibility and will only drop sequences that contain missing values for required columns. The user may wish to filter to only productive sequences in advance. Some scripts allow the specification of a group column. In these cases this column will be used to group sequences by each unique level of this column and the resulting metrics will be averaged. Otherwise metrics are calculated from the provided file as a whole, treating every sequence equally.

Details of the input files and the scripts that can be used to process them are included below. These scripts do require the installation of pandas (tested with v1.3.5)

## IMGT alleles

Germline VDJ alleles in FASTA format with IMGT headers. Only alleles that are marked as F (functional) or ORF (open reading frame) will be used for sequence generation:

* IMGT_human_IGHV.fasta
* IMGT_human_IGHD.fasta
* IMGT_human_IGHJ.fasta

AIRRSHIP includes the IMGT dataset as of 18.02.2022. If an updated dataset is used, or the user wishes to use an alternate database, all other files must contain information to support all the genes present in the new dataset (i.e. if IMGT adds new genes or gene families, and you try to substitute in these files, AIRRSHIP may fail as the currently included distributions do not account for them). AIRRSHIP currently only supports alleles with names in the IMGT format,

See Non-human Species for allele files for information on simulating sequences from species other than humans.

!!! note
    To avoid confusion when benchmarking alignments, AIRRSHIP ignores the IMGT allele IGHV3-30-3\*03 as it has the same sequence as allele IGHV3-30\*04. AIRRSHIP will also not utilise V or J alleles that do not have the expected junctional anchors.

## VDJ usage

VDJ gene and family usage, giving the proportion of sequences to include each gene or family:

* IGHV_usage.csv
* IGHV_usage_gene.csv
* IGHD_usage.csv
* IGHD_usage_gene.csv
* IGHJ_usage.csv
* IGHJ_usage_gene.csv

###### Example: 

| v_family   | prop        |
| :--------- | :---------- |
| IGHV1      | 0.210282513 |
| IGHV1/OR15 | 0.000783292 |
| ...        | ...         |


#### Scripts:

```bash
create_vdj_files.py  [-h] -i INPUT_FILE [--group GROUP]

-i INPUT_FILE, --input_file INPUT_FILE
--group GROUP   Name of column containing dataset metadata to group by.
```

Required columns: v_call, d_call, j_call


## Trimming

Distributions of number of nucleotides to be removed from gene ends at the junctions, per gene family:

* V_family_trimming_proportions.csv
* D_5_family_trimming_proportions.csv
* D_3_family_trimming_proportions.csv
* J_family_trimming_proportions.csv

###### Example: 

| v_family | v_3p_del | proportions |
| :------- | :------- | :---------- |
| IGHV1    | 0        | 0.406392257 |
| IGHV1    | 1        | 0.290851348 |
|  ...     | ...      | ...         |

#### Scripts:

```bash
create_trimming_files.py [-h] -i INPUT_FILE [--imgt_v IMGT_V] [--imgt_d IMGT_D] 

-i INPUT_FILE, --input_file INPUT_FILE
--imgt_v IMGT_V Path to IMGT V allele file.
--imgt_d IMGT_D Path to IMGT D allele file.
```

Required columns:  v_call, d_call, v_germline_end, j_germline_start, d_germline_start, d_sequence_end and d_sequence_start.

## NP additions

Proportion of NP regions starting with each nucleotide base:

* np1_first_base_probs.csv
* np2_first_base_probs.csv

###### Example:

|     | proportion  |
| :-- | :---------- |
| A   | 0.111702543 |
| C   | 0.246122379 |
| G   | 0.284580985 |
| T   | 0.357594093 |

Proportion of NP regions of each length:

* np1_lengths_proportions.csv
* np2_lengths_proportions.csv

###### Example:

| np1_length | prop        |
| :--------- | :---------- |
| 0          | 0.052651438 |
| 1          | 0.040751106 |
| ...        | ...         |

Transition probabilities for each position in NP region (i.e. the likelihood of the next base being A,C,G or T when the current base is A,C,G or T). For both mutated and unmutated sequences:

* np1_transition_probs_per_position_igdm.csv
* np1_transition_probs_per_position_igag.csv
* np2_transition_probs_per_position_igdm.csv
* np2_transition_probs_per_position_igag.csv

###### Example:

| Length | Base | A           | C           | G           | T           |
| :----- | :--- | :---------- | :---------- | :---------- | :---------- |
| 0      | T    | 0.169400291 | 0.390251379 | 0.232798504 | 0.207549826 |
| 0      | A    | 0.24858057  | 0.238914924 | 0.34526599  | 0.167238516 |
| ...    | ...  | ...         | ...         | ...         | ...         |

#### Scripts:

```bash
create_np_lengths_files.py [-h] -i INPUT_FILE [--group GROUP]  

-i INPUT_FILE, --input_file INPUT_FILE

--group GROUP   Name of column containing dataset metadata to group by.
```

Required columns: np1_length, np2_length

```bash
create_np_first_bases.py [-h] -i INPUT_FILE

-i INPUT_FILE, --input_file INPUT_FILE
```

Required columns: np1_length,np2_length,sequence,v_sequence_end,d_sequence_start and j_sequence_start

```bash
create_np_transitions_files.py [-h] -i INPUT_FILE

-i INPUT_FILE, --input_file INPUT_FILE
```

Required columns: np1_length,np2_length,sequence,v_sequence_end,d_sequence_start and j_sequence_start

When producing NP transition matrices, only a single file for each NP region is produced. AIRRSHIP requires two, np[1|2]_transition_probs_per_position_igdm.csv and np[1|2]_transition_probs_per_position_igag.csv, to account for differences in inserted nucleotides following SHM. Therefore this script must be run once on unmutated sequences and once on mutated sequences and the files manually renamed. If required, the same file with different names can be used although this may affect insertion accuracy.


## Somatic Hypermutation

Rates of mutation frequency per sequence, per V family:

* mut_freq_per_seq_per_family.csv

###### Example:

| v_family | mut_freq    | proportion  |
| :------- | :---------- | :---------- |
| IGHV1    | 0           | 0.005959742 |
| IGHV1    | 0.1         | 0.003065621 |
| IGHV1    | 0.111111111 | 0.002811285 |
| ...       | ...        | ...         |

Per base mutation rates for the centre base of each unique 5mer (the likelihood of the centre base being one of A,C,G,T for that kmer when pulled from mutated sequences). For each IMGT region and for the CDR and FWR regions in aggregate:

* cdr_kmer_base_usage.csv
* cdr1_kmer_base_usage.csv
* cdr2_kmer_base_usage.csv
* cdr3_kmer_base_usage.csv
* fwr_kmer_base_usage.csv
* fwr1_kmer_base_usage.csv
* fwr2_kmer_base_usage.csv
* fwr3_kmer_base_usage.csv
* fwr4_kmer_base_usage.csv

##### Example:

| kmer  | A           | C           | G           | T           |
| :---- | :---------- | :---------- | :---------- | :---------- |
| GTGCA | 0.054353444 | 0.001915083 | 0.927897253 | 0.015834221 |
| TGCAC | 0.001906821 | 0.901529425 | 0.017524232 | 0.079039522 |
| GCACA | 0.847932735 | 0.043780581 | 0.074638153 | 0.033648531 |
| ...   | ...         | ...         | ...         | ...         |

#### Scripts:

```bash
create_shm_files.py [-h] -i INPUT_FILE

-i INPUT_FILE, --input_file INPUT_FILE
```

Required columns: v_call, germline_alignment, sequence_alignment, fwr1, fwr2, fwr3, cdr1, cdr2, cdr3

## Experimental data used

The input files included with AIRRSHIP were generated using repertoire data from the following public datasets:

1. Cowan, G et al., unpublished data.
2.	Ghraichy,M. et al. (2020) Maturation of the Human Immunoglobulin Heavy Chain Repertoire With Age. Front. Immunol., 11, 1734.
    * Preprocessed data downloaded from https://zenodo.org/record/3585046#.Y2upy4LP2JF.
    * Only individuals aged > 9 years included.
3.	Gidoni,M. et al. (2019) Mosaic deletion patterns of the human antibody heavy chain gene locus shown by Bayesian haplotyping. Nat. Commun., 10, 628.
    * Raw sequence data downloaded from the ENA, accession PRJEB26509.
    * Healthy controls only included.
4.	Waltari,E. et al. (2018) 5′ Rapid Amplification of cDNA Ends and Illumina MiSeq Reveals B Cell Receptor Features in Healthy Adults, Adults With Chronic HIV-1 Infection, Cord Blood, and Humanized Mice. Front. Immunol., 9, 628.
    * Raw sequence data downloaded from the SRA, accession PRJNA393446.
    * Healthy adults only included.
5.	Yang,X. et al. (2021) Large-scale analysis of 2,152 Ig-seq datasets reveals key features of B cell biology and the antibody repertoire. Cell Rep., 35, 109110.
    * Raw sequence data downloaded from the SRA, accession PRJNA564936.

Sequences were processed using pRESTO [1] and Change-O [2] from the Immcantation pipeline, with VDJ assignment carried out using IgBLAST [3] v1.18.0.

[1]	Vander Heiden, J. A. et al. pRESTO: a toolkit for processing high-throughput sequencing raw reads of lymphocyte receptor repertoires. Bioinformatics 30, 1930–1932 (2014). 

[2]	Gupta, N. T. et al. Change-O: a toolkit for analyzing large-scale B cell immunoglobulin repertoire sequencing data. Bioinformatics 31, 3356–3358 (2015).

[3]	Ye, J., Ma, N., Madden, T. L. & Ostell, J. M. IgBLAST: an immunoglobulin variable domain sequence analysis tool. Nucleic Acids Res. 41, W34–W40 (2013).




