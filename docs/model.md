# Simulation Model

## 1. Locus creation <a name="haplotype"></a>

AIRRSHIP takes IMGT FASTA files of VDJ alleles as input. These files are used to populate artificial heavy chain loci, each including single alleles representing every V, D and J gene. If more than two alleles exist for that gene, the allele included on each chromosome may differ (resembling heterozygous carriage) and the user can control the proportion of genes for which this is desired. By default, this means no one repertoire will contain sequences formed from more than two alleles of a single gene. However, the user can override this requirement and all alleles present in the input dataset will be used (using ```--all_alleles```). The proportion of heterozygous positions can be controlled by the ```--het``` argument. 

It is possible to override this and use all alleles present in the input dataset (specify with ```--all_alleles```).

## 2. VDJ recombination <a name="vdj"></a>

To generate the initial recombined sequence, VDJ alleles are chosen based on gene usage distributions established from published datasets produced by 5’ RACE amplification. Gene segments are recombined from the same synthetic loci only. 

VDJ usage may also be requested to be flat at either the gene or family level (```--flat_vdj gene``` or ```--flat_vdj family```).

## 3. Trimming of gene ends <a name="trim"></a>

The 3’ end of the V gene, 5’ and 3’ end of the D gene and 5’ end of the J gene are then trimmed. For each gene end, the number of nucleotides to be removed is sampled from distributions of trimming lengths for each IMGT gene family. To maintain sequence productivity, trimming lengths that would extend past the C-104 anchor in the V gene or the W/F - 118 anchor in the J gene, or that would result in complete deletion of the D gene, are resampled. 

Trimming may be turned off at any or all of these gene ends (```--no_trim```, ```--no_trim_v3```, ```--no_trim_d5```, ```--no_trim_d3``` or ```--no_trim_j5```).

## 4. Addition of NP nucleotides <a name="np"></a>

AIRRSHIP does not distinguish between N and P nucleotides when insertions are modelled at the VD (NP1) and DJ (NP2) junctions. For each insertion, the number of nucleotides to be added is sampled from distributions of NP1 or NP2 lengths. The first nucleotide to be inserted is then randomly selected with probabilities for each base determined from public data for NP1 and NP2. Addition of further nucleotides follows a position dependent Markov process. For each position in the NP region, the next base is chosen with likelihood determined by a transition matrix that considers the current position in the sequence and the nucleotide base present at that position. When simulating non-mutated sequences and hypermutated sequences, a transition matrix determined from IgD/IgM sequences and from IgA/IgG sequences respectively, is used. This compensates for the inability to determine which positions within an NP region have been hypermutated.

NP addition may be turned off at either or both of the junctions (```--no_np```, ```--no_np1``` or ```--no_np2```).

## 5. Somatic hypermutation (SHM) <a name="shm"></a>

Somatic hypermutation (SHM) is replicated at both the per sequence and per nucleotide position level. For each sequence, the overall mutation frequency is chosen from a distribution established from published datasets produced by 5’ RACE amplification methods. Each 5mer within the sequence, excluding those where the centre base is an NP nucleotide, is then randomly iterated over. A new base is chosen to replace the centre base of this 5mer, sampled from a set of distributions that give the frequency of each nucleotide occurring at the centre position of each unique kmer for each sequence region (e.g. FWR1, CDR1). If the base chosen differs from the germline base, a mutation is introduced at that position and the process repeats until the number of mutations required to give the desired mutation frequency is reached. 

A flat per sequence mutation rate may also be requested (using ```--shm_flat``` and specified by either ```--mut_num``` or ```--mut_rate```). Each sequence in the repertoire will have the same number or frequency of mutations.

When ```--shm_random``` is specified, per position mutation is not explicitly modelled. Instead, each base is mutated randomly, irrespective of it's regional or kmer context. This may be combined with ```--shm_flat```.

## 6. Removal of non-productive sequences <a name="productive"></a>

AIRRSHIP aims to produce only productive BCR sequences, defined here as those where the V and J segments are in-frame, no stop codon is present and the correct junction anchor residues (C-104 and W/F-118) are present. Alleles which do not have the correct anchor residues are excluded when processing the input data. Checks for productivity occur following trimming and NP addition. As certain VDJ combinations may be more likely to result in non-productive rearrangements, multiple attempts are made at trimming and nucleotide insertion using the same allele set to maintain VDJ usage distributions. SHM also commonly renders sequences non-productive, so multiple attempts at introducing hypermutation at the same rate will be made to ensure the per-sequence hypermutation reflects published datasets.