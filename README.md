# AIRRSHIP - Adaptive Immune Receptor Repertoire Simulation of Human Immunoglobulin Production

AIRRSHIP simulates B cell receptor (BCR) sequences for use in benchmarking applications where BCR sequences of known origin are required.

AIRRSHIP replicates the VDJ recombination process from haplotype through to somatic hypermutation. Recombination metrics are derived from a range of experimental sequences allowing faithful replication of true repertoires. Users may also control a wide range of parameters that influence allele usage, junctional diversity and somatic hypermutation rates. The current model extends to human heavy chain BCR sequences only. 

## Installation

```
pip install airrship
```

## Documentation

Full documentation is available [here](https://airrship.readthedocs.io). 

## Example datasets

Four pre-generated datasets, each with 500,000 sequences, are available at [GitHub](https://github.com/Cowanlab/airrship). These were produced using default parameters, two with and two without SHM, and may be used if you do not wish to install and run AIRRSHIP yourself or want to explore the output yourself first.