# Quickstart

## Download

The latest version of AIRRSHIP can be downloaded from PyPi or [GitHub](https://github.com/Cowanlab/airrship).

## Installation

The easiest way to install is using pip, either directly:

```bash
pip install airrship
```

Or after downloading the latest release:

```bash
pip install airrship-x.y.z.tar.gz
```
## Requirements

AIRRSHIP intentionally uses only Python standard libraries and requires only the installation of base Python (version 3.7 or above).

## Examples

A very small example repertoire is held at the AIRRSHIP [GitHub repository](https://github.com/Cowanlab/airrship) to provide an example of the expected output. Larger example repertoire files are available at [Zenodo](https://doi.org/10.5281/zenodo.7568251).

## Running from the command line <a name="command line"></a>


The most basic call to AIRRSHIP requires only an output name.

```bash
airrship -o my_repertoire
```
This will create a repertoire of 1000 unmutated human, heavy chain BCR sequences with metrics derived from experimental distributions.

Four output files will be generated: 

* my_repertoire.fasta - final sequences in FASTA format
* my_repertoire.tsv - information regarding sequence generation
* my_repertoire_locus.csv - the simulated locus
* my_repertoire_summary.txt - summary of input commands 

Please see [Output Files](output.md) for further details on output file format.

### Customising repertoire generation

By default, AIRRSHIP attempts to replicate real experimental repertoires as closely as possible. However, there a large number of command line options that can be used to produce repertoires with specific desired features. 

For example, we could create a repertoire with:

* 16,000 sequences 
* a locus where every gene is homozygous (only one allele present per gene)
* balanced usage of the gene families
* no trimming of the 5' end of the D gene
* no insertion of nucleotides between the D and J gene (no NP2 regions)
* include non-productive sequences in output and limit them to 10% of the repertoire

```bash
airrship -o complex_repertoire \
                     -n 16000 \
                     --het 0 0 0 \
                     --flat_vdj family \
                     --no_trim_d5 \
                     --no_np2 \
                     --non_productive \
                     --prop_non_productive 0.1
```

Full details can be found in [Command line Usage](parameters.md).


!!! note
    Occasionally AIRRSHIP may fail to generate a productive sequence from a specific combination of alleles and will print a warning. This should not be of concern unless it happens with high frequency. In this case you may need to check your chosen parameters or input data.


### Adding somatic hypermutation


The ```--shm``` flag will generate SHM according to observed distributions (see [Simulation Model](model.md) for more information). 

```bash
airrship -o shm_repertoire --shm
```

Mutation rates can be controlled by passing a multiplication factor with  ```--shm_multiplier```. For example, the below command will create a repertoire with sequences mutated to rates half that as specified in the reference files.

```bash
airrship -o shm_repertoire --shm --shm_multiplier 0.5
```

To request a constant mutation frequency across all sequences, the ```--shm_flat ``` option can be used. The desired mutation rate or number can be specified with either ```--mut_rate``` or ```mut_num```.

The below command will create 1000 sequences, each of which with a mutation rate of 0.08 (i.e. number of mutations in sequence / length of sequence = 0.08). The distribution will be as close to flat as is possible but may fluctuate slightly.

```bash
airrship -o shm_flat_repertoire --shm_flat --mutation_rate 0.08
```

The default SHM algorithm treats each base in the sequence differently, depending on the 5mer context of the base and the region of the sequence it is found in. To make per base mutation independent of sequence context, ```--shm_random``` can be used. 

```bash
airrship -o shm_random_repertoire --shm_random
```

The per sequence mutation rate will still follow the observed experimental distribution unless ```--shm_flat``` is also specified. 

!!! note
    Setting mutation rates higher than 0.2 will result in a warning and, depending on the other options specified, may result in very slow performance or a failure to generate sufficient sequences. Other distributions may also be skewed. Mutation rates above 0.5 are not supported.


## Using the package in Python <a name="python"></a>

If desired, instead of running from the command line, the package can be imported and a call to main() made within Python, specifying the same parameters as discussed above.

```python
from airrship import create_repertoire

create_repertoire.main(['-o', 'my_repertoire', '--outdir', 'output'])

create_repertoire.main(('-o my_repertoire --outdir output').split())
```

It is also possible to use individual functions from the package. A simple three step workflow to generate sequences is described below.


#### Read in data to be used

```python
from airrship import create_repertoire

data_dict = create_repertoire.load_data()
```

#### Create a locus from which to generate sequences

```python
locus = create_repertoire.get_genotype()
```

#### Generate sequences

```python
sequence = create_repertoire.generate_sequence(locus, data_dict, mutate = True)
```

```generate_sequence``` returns an individual Sequence class object with all information about its generation stored as attributes. For example, the final mutated sequence can be accessed using sequence.mutated_seq. For full details, see [Python](python_use.md).

