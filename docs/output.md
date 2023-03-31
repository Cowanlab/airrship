# Output files

When used at the command line, AIRRSHIP produces four output files per run:

### 1. Sequence FASTA <a name="fasta"></a>

**outname.fasta**

A FASTA file containing the final simulated sequences. In the case of SHM,
these will be the mutated sequences. The FASTA headers correspond to the sequence_id column in outname.csv.

### 2. Sequence information TSV <a name="tsv"></a>

**outname.tsv**

A tab-delimited file containing information about each sequence and its formation. The file format follows the [AIRR-C Rearrangement Schema](https://docs.airr-community.org/en/stable/datarep/rearrangements.html) where possible. 

The below columns are present regardless of simulation criteria:

| **Name**         | **Description**                                                                 |
| :--------------- | :------------------------------------------------------------------------------ |
| sequence_id      | Unique sequence identifier.                                                     |
| sequence         | Final simulated nucleotide sequence.                                            |
| productive       | True if sequence is predicted to be productive.                                 |
| stop_codon       | True if the sequence contains a stop codon.                                     |
| vj_in_frame      | True if the V and J segments are in frame.
| v_call           | V gene with allele.                                                             |
| d_call           | D gene with allele.                                                             |
| j_call           | J gene with allele.                                                             |
| junction         | Junction region nucleotide sequence. CDR3 plus two conserved codons.            |
| junction_aa      | Junction region amino acid translation.                                         |
| junction_length  | Length of the junction region.                                                  |
| np1_length       | Length of the combined N/P region between the V and D gene.                     |
| np1              | Nucleotide sequence of the combined N/P region between the  <br/> V and D gene. |
| np2_length       | Length of the combined N/P region between the D and J gene.                     |
| np2              | Length of the combined N/P region between the D and J gene.                     |
| v_3_trim         | Number of nucleotides trimmed from the 3' end of the V gene.                    |
| d_5_trim         | Number of nucleotides trimmed from the 5' end of the D gene.                    |
| d_3_trim         | Number of nucleotides trimmed from the 3' end of the D gene.                    |
| j_5_trim         | Number of nucleotides trimmed from the 5' end of the J gene.                    |
| v_sequence       | Part of the sequence originating from the V gene.                               |
| d_sequence       | Part of the sequence originating from the D gene.                               |
| j_sequence       | Part of the sequence originating from the J gene.                               |
| v_sequence_start | Start position of the V gene in the sequence (1-based closed interval).         |
| v_sequence_end   | End position of the V gene in the sequence (1-based closed interval).           |
| d_sequence_start | Start position of the D gene in the sequence (1-based closed interval).         |
| d_sequence_end   | End position of the D gene in the sequence (1-based closed interval).           |
| j_sequence_start | Start position of the J gene in the sequence (1-based closed interval).         |
| j_sequence_end   | End position of the J gene in the sequence (1-based closed interval).           |

Some columns are present only when SHM is not simulated:

| **Name**        | **Description**                                                                   |
| :-------------- | :-------------------------------------------------------------------------------- |
| gapped_sequence | Simulated nucleotide sequence with gaps inserted according to IMGT  <br/> schema. |

Some columns are present only when SHM is simulated:

| **Name**                  | **Description**                                                                            |
| :------------------------ | :----------------------------------------------------------------------------------------- |
| shm_events                | Comma-delimited list of mutation events. In the format <br/> position:base>mutated_base    |
| shm_count                 | Number of mutations in the sequence.                                                       |
| shm_freq                  | Mutation frequency (number of mutations divided by  <br/> length of sequence)              |
| unmutated_sequence        | Unmutated simulated nucleotide sequence.                                                   |
| gapped_unmutated_sequence | Unmutated simulated nucleotide sequence with gaps inserted <br/> according to IMGT schema. |
| gapped_mutated_sequence   | Mutated simulated nucleotide sequence with gaps inserted <br/> according to IMGT schema.   |


### 3. Locus file <a name="locus"></a>

**outname_locus.csv**

A two column CSV file containing the alleles chosen for each simulated "chromosome". Can be used in subsequent runs to simulate sequences from the same genetic background.

### 4. Summary file <a name="summary"></a>

**outname_summary.txt**

A text file listing the arguments provided to the AIRRSHIP call.
