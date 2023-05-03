# Simulating Non-human Species

AIRRSHIP has not specifically been designed to allow the simulation of sequences from non-human sequences and does not come with in-built reference files for other species. However, if the species of interest has a recombination process that follows the same steps as are present in humans, it should be possible for a user to simulate from a reference directory containing files generated from this species.

The ```--species``` flag must be used to find the VDJ allele files if simulating non-human sequences. File names are expected to be in the format imgt_[species]_IGH[V|D|J].fasta. 

AIRRSHIP is primarily intended for simulation of human sequences and has been tested on this basis. If you are simulating non-human sequences, ensure the output looks as you would expect. Any unusual behaviour or errors by can be reported by raising an issue on Github. 


## Mouse Reference Datasets

We have tested this premise by simulating sequences from C57BL/6 mice and reference files are available at the AIRRSHIP [Github](https://github.com/Cowanlab/airrship/c57bl6_reference). These reference files were produced from the datasets detailed below. The germline alleles included are those that are present in both the IMGT mouse reference files and the OGRDB C57BL/6 reference germline set [1,2]. The number of individuals and sequences is considerably smaller than that used to generate the human references. We have also not fully optimised AIRRSHIP for non-human data, and as a result the simulations are slightly less realistic than for human sequences, especially across the junction region. 

To simulate using these files, the following command would be used after the c57bl6 reference directory had been downloaded:

```bash
airrship -o c57bl6_repertoire \
            --datadir path/to/c57bl6_reference \
            --species mouse
```

## Experimental data used

1. Collins, A.M. et al. (2015) The mouse antibody heavy chain repertoire is germline-focused and highly variable between inbred strains. Philos. Trans. R. Soc. B Biol. Sci., 370.
    * Raw sequence data downloaded from ENA accession PRJEB8745 
    * Only C57BL/6 samples used
2. Corcoran, M.M. et al. (2016) Production of individualized V gene databases reveals high levels of immunoglobulin genetic diversity. Nat. Commun. 2016 71, 7, 1–14.
    * Raw sequence data downloaded from ENA accession PRJEB15295
    * C57Bl/6 M3 sample used
3. Greiff, V. et al. (2017) Systems Analysis Reveals High Genetic and Antigen-Driven Predetermination of Antibody Repertoires throughout B Cell Development. Cell Rep., 19, 1467–1478.
    * Raw sequence data downloaded from ArrayExpress accession E-MTAB-5349
    * Naive splenocyte samples from uninfected C57BL/6 used
4. Kräutler, N.J. et al. (2020) Quantitative and Qualitative Analysis of Humoral Immunity Reveals Continued and Personalized Evolution in Chronic Viral Infection. Cell Rep., 30, 997-1012.e6.
    * Raw sequence data downloaded from ArrayExpress accession E-MTAB-8585
5. Mouat, I. et. al., unpublished data
    * Splenocytes from eight control C57BL/6 mice


[1]	Lees, W. et al. OGRDB: a reference database of inferred immune receptor genes. Nucleic Acids Res. 48, D964–D970 (2020).

[2]	Jackson, K. J. L. et al. A BALB/c IGHV Reference Set, Defined by Haplotype Analysis of Long-Read VDJ-C Sequences From F1 (BALB/c x C57BL/6) Mice. Front. Immunol. 13, 2490 (2022).
