#!/usr/bin/env python

'''AIRRSHIP - Simulation of B cell receptor repertoires

Copyright (C) 2022  Catherine Sutherland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.'''


import random
from collections import defaultdict
import csv
import time
import warnings
import re
import argparse
import sys
import importlib.resources

# Classes


class Allele:
    """Class that represents a V, D or J allele.

    Attributes:
        name (str): The IMGT name of the allele
        gapped_seq (str): The IMGT gapped germline nucleotide sequence
        length (str): IMGT defined length of the allele
        ungapped_sq (str): Ungapped germline nucleotide sequence
        trim_5 (int): Number of nucleotides to be trimmed from 5' end
        trim_3 (int): Number of nucleotides to be trimmed from 3' end

    """

    def __init__(self, name, gapped_seq, length):
        """Initialises an Allele class instance.

        Args:
            name (str): The IMGT name of the allele
            gapped_seq (str): The IMGT gapped nucleotide sequence
            length (str): IMGT defined length of the allele
        """
        self.name = name
        self.segment = name.split("IGH", 1)[1][0]
        self.gapped_seq = gapped_seq
        self.length = int(length)
        self.ungapped_seq = gapped_seq.replace(".", "")
        self.ungapped_len = len(self.ungapped_seq)
        if "-" in self.name:
            self.family = self.name.split("-")[0]
        else:
            self.family = self.name.split("*")[0]
        self.gene = self.name.split("*")[0]
        if self.segment == "V":
            # look for 2nd Cys at position 104
            cys = self.gapped_seq[309:312]
            cys_wider = self.gapped_seq[306:315]
            self.anchor = self.ungapped_seq.rfind(cys_wider) + 3
        if self.segment == "J":
            # look for F/W - G - X - G motif
            motif = re.compile(
                '(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
            match = motif.search(self.ungapped_seq)
            self.anchor = match.span()[0]

    def __repr__(self):  # defines how the class object is represented
        return f"{self.__class__.__name__}({self.name}, {self.gapped_seq[:10]}, {self.length})"

    def get_trim_length(self, no_trim_list, trim_dicts):  # use per family trimming
        """Chooses trimming lengths for allele.

        Adds two class attributes - trim_3, 3' prime trimming value and trim_5,
        5' prime trimming value.

        Args:
            no_trim_list (list): List of 5 Booleans, specifying whether to not
                trim [all_ends, v_3_end, d_5_end, d_3_end, j_5_end]. 
            trim_dicts (dict): A dictionary of dictionaries of trimming length 
                proportions by gene family for each segment (V, D or J).
        """

        no_all, no_v3, no_d5, no_d3, no_j5 = no_trim_list

        trim_3 = 0  # set to 0 - J will never be trimmed at 3'
        trim_5 = 0  # set to 0 - V will never be trimmed at 5'

        if no_all == False:
            if self.segment == "V":
                if no_v3 == False:
                    trim_3_dict = trim_dicts["V_3"]
                    # choose trim length/prob dict by gene family
                    if self.family in trim_3_dict:
                        prob_dict = trim_3_dict[self.family]
                    else:
                        prob_dict = random.choice(list(trim_3_dict.values()))
                    lengths, probs = zip(*prob_dict.items())
                    trim_3 = choice(lengths, probs)
                    # trim_3 = weighted_choice(prob_dict.items())
                    # prevent entire allele or anchor from being removed
                    this_loop = 0
                    while (trim_3 >= self.length) or (trim_3 >= (self.length - self.anchor - 1)):
                        if this_loop > 50:
                            print("failed to trim", self.name)
                            break
                        this_loop += 1
                        trim_3 = choice(lengths, probs)

            elif self.segment == "D":
                if no_d5 == False:
                    trim_5_dict = trim_dicts["D_5"]
                    if self.family in trim_5_dict:
                        prob_5_dict = trim_5_dict[self.family]
                    else:
                        prob_5_dict = random.choice(list(trim_5_dict.values()))
                    lengths_5, probs_5 = zip(*prob_5_dict.items())
                    trim_5 = choice(lengths_5, probs_5)

                if no_d3 == False:
                    trim_3_dict = trim_dicts["D_3"]
                    if self.family in trim_3_dict:
                        prob_3_dict = trim_3_dict[self.family]
                    else:
                        prob_3_dict = random.choice(list(trim_3_dict.values()))
                    lengths_3, probs_3 = zip(*prob_3_dict.items())
                    trim_3 = choice(lengths_3, probs_3)

                while trim_3 + trim_5 >= self.length:
                    if no_d5 == False:
                        trim_5 = choice(lengths_5, probs_5)

                    if no_d3 == False:
                        trim_3 = choice(lengths_3, probs_3)

            elif self.segment == "J":
                if no_j5 == False:
                    trim_5_dict = trim_dicts["J_5"]
                    if self.family in trim_5_dict:
                        prob_dict = trim_5_dict[self.family]
                    else:
                        prob_dict = random.choice(list(trim_5_dict.values()))
                    lengths, probs = zip(*prob_dict.items())
                    trim_5 = choice(lengths, probs)
                    # trim_5 = weighted_choice(prob_dict.items())

                    while (trim_5 >= self.length) or (trim_5 >= self.anchor):
                        trim_5 = choice(lengths, probs)
                        # trim_5 = weighted_choice(prob_dict.items())

        self.trim_5 = trim_5
        self.trim_3 = trim_3


class Sequence:
    """
    Represents a recombined Ig sequence consisting of V, D and J segments.

    Attributes:
        v_allele (Allele): IMGT V gene allele.
        d_allele (Allele): IMGT D gene allele.
        j_allele (Allele): IMGT J gene allele.
        alleles (list): List of IMGT alleles.
        NP1_region (str): NP1 region - between V and D gene.
        NP1_length (int): Length of NP1 region.
        NP2_region (str): NP2 region - between V and D gene.
        NP2_length (int): Length of NP2 region.
        ungapped_seq (str): Ungapped nucleotide sequence.
        gapped_seq (str): Gapped nucleotide sequence.
        mutated_seq (str): Ungapped mutated nucleotide sequence.
        gapped_mutated_seq (str): Ungapped mutated nucleotide sequence.
        mutated_seq (str): Ungapped mutated nucleotide sequence.
        junction (str): Nucleotide sequence of junction region.
        v_seq (str):  Nucleotide sequence of V region.
        d_seq (str): Nucleotide sequence of D region.
        j_seq (str): Nucleotide sequence of J region.
        v_seq_start (int): Start position of V region.
        d_seq_start (int): Start position of D region.
        j_seq_start (int): Start position of J region.
        v_seq_end (int): End position of V region.
        d_seq_end (int): End position of D region.
        j_seq_end (int): End position of J region.
        mutations (str): Mutation events.
        mut_count (int): Mutation count.
        mut_freq (int): Mutation frequency.
    """

    def __init__(self, v_allele, d_allele, j_allele):
        """Initialises a Sequence class instance.

        Args:
            v_allele (Allele): IMGT V gene allele, required.
            d_allele (Allele): IMGT D gene allele, required.
            j_allele (Allele): IMGT J gene allele, required.
        """
        self.v_allele = v_allele
        self.d_allele = d_allele
        self.j_allele = j_allele
        self.alleles = [v_allele, d_allele, j_allele]
        self.NP1_region = ""
        self.NP2_region = ""
        self.junction = ""
        self.v_seq = ""
        self.d_seq = ""
        self.j_seq = ""
        self.v_seq_start = 0
        self.d_seq_start = 0
        self.j_seq_start = 0
        self.v_seq_end = 0
        self.d_seq_end = 0
        self.j_seq_end = 0
        self.mutations = ""
        self.mut_count = 0
        self.mut_freq = 0
        self.ungapped_seq = ""
        self.mutated_seq = None
        self.gapped_seq = ""
        self.gapped_mutated_seq = None

    def get_nuc_seq(self, no_trim_list, trim_dicts, no_np_list, NP_lengths, NP_transitions, NP_first_bases, gapped=False):
        """Creates the recombined nucleotide sequence with trimming and np addition.

        Args:
            no_trim_list (list): List of 5 Booleans, specifying whether to not
                trim [all_ends, v_3_end, d_5_end, d_3_end, j_5_end]. 
            trim_dicts (dict): A dictionary of dictionaries of trimming length 
                proportions by gene family for each segment (V, D or J).
            no_np_list (list): List of 3 Booleans, specifying whether to not
                add [both_np, np1, np2].
            NP_lengths (dict): Dictionary of possible NP region lengths and the
                proportion of sequences to use them. In the format
                {NP region length: proportion}.
            NP_transitions (dict): Nested dictionary containing transition matrix of
                probabilities of moving from one nucleotide (A, C, G, T) to any other 
                for each position in the NP region.
            NP_first_bases (dict): Nested dictionary of the proportion of NP 
                sequences starting with each base for NP1 and NP2.
                gapped (bool): Specify whether to return sequence with IMGT gaps
                or not.

        Returns:
            nuc_seq (str): The recombined nucleotide sequence.
        """

        v_allele = self.v_allele
        d_allele = self.d_allele
        j_allele = self.j_allele

        no_NP, no_NP1, no_NP2 = no_np_list

        if no_NP == False:
            if no_NP1 == False:
                self.NP1_region = get_NP_regions_markov_pos(
                    NP_lengths, NP_transitions, "NP1", NP_first_bases)
            if no_NP2 == False:
                self.NP2_region = get_NP_regions_markov_pos(
                    NP_lengths, NP_transitions, "NP2", NP_first_bases)

        self.NP1_length = len(self.NP1_region)
        self.NP2_length = len(self.NP2_region)

        if gapped == True:

            nuc_seq = (
                get_gapped_trimmed_seq(v_allele, no_trim_list, trim_dicts)
                + self.NP1_region
                + get_gapped_trimmed_seq(d_allele, no_trim_list, trim_dicts)
                + self.NP2_region
                + get_gapped_trimmed_seq(j_allele, no_trim_list, trim_dicts)
            )

        if gapped == False:

            nuc_seq = (
                get_trimmed_seq(v_allele, no_trim_list, trim_dicts)
                + self.NP1_region
                + get_trimmed_seq(d_allele, no_trim_list, trim_dicts)
                + self.NP2_region
                + get_trimmed_seq(j_allele, no_trim_list, trim_dicts)
            )

        return nuc_seq

    def get_junction_length(self):
        """Calculates the junction length of the sequence (CDR3 region plus both
        anchor residues).

        Returns:
            junction_length (int): Number of nucleotides in junction (CDR3 + anchors)
        """

        junction_length = self.v_allele.length - (self.v_allele.anchor - 1) - self.v_allele.trim_3 + self.NP1_length + self.d_allele.length - \
            self.d_allele.trim_5 - self.d_allele.trim_3 + \
            self.NP2_length + (self.j_allele.anchor + 2) - self.j_allele.trim_5
        return junction_length


class NP_Region:
    """
    Class that represents an NP region where base usage is determined by a first
    order Markov chain for which the transition matrix varies by position.

    Attributes:
        transition_probs (dict) : Dictionary of transition matrices per position
            in sequence.
        first_base (str) : The first base in the NP region.
        length (int) : The chosen length of the NP region.
    """

    def __init__(self, transition_probs, first_base, length):
        """Initialises an NP_region class instance.
        """

        self.transition_probs = transition_probs
        self.bases = ["A", "C", "G", "T"]
        self.first_base = first_base
        self.length = length

    def next_base(self, current_base, position):
        """Get the next base in NP region from transition matrix.

        Args:
            current_base (str): Current base in NP region sequence.
            position (int): Current position in NP region sequence.

        Returns:
            base (str): Next base in the NP region sequence.
        """

        base = choice(
            self.bases, [self.transition_probs[position][current_base][next_base]
                         for next_base in self.bases])
        return base

    def generate_np_seq(self):
        """Creates an NP region sequence using a first order Markov chain.

        Returns:
            sequence (str): Final NP region sequence.
        """
        sequence = ""
        current_base = self.first_base
        sequence += current_base
        for i in range(self.length - 1):
            next_base = self.next_base(current_base, i)
            sequence += next_base
            current_base = next_base
        return sequence


# Helper functions

def choice(options, probs):
    """Choose an element from a list according to a list of corresponding
    (cumulative) probabilities.

    Args:
        options (list): List of options to be selected from.
        probs (list): List of probabilities of each option.

    Returns:
        List element chosen from options, according to the probabilities given.
    """
    x, c = random.random(), 0
    for i, p in enumerate(probs):
        c += p
        if x < c:
            break
    return options[i]


def translate(seq):
    """Translates a nucleotide sequence to an amino acid sequence.

    Args:
        seq (str): Nucleotide sequence to be translated.

    Returns:
        protein (str): Amino acid sequence.
    """

    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }
    protein = ""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3].upper()
        if len(codon) < 3:
            protein += '.'
        elif "." in codon:
            protein += '_'
        elif "-" in codon:
            protein += 'X'
        elif "N" in codon:
            protein += 'X'
        else:
            protein += table[codon]
    return protein


def check_stops(sequence):
    """Check for stop codons in a nucleotide sequence.

    Args:
        sequence (str): Nucleotide sequence.

    Returns:
        stop (bool): True if stop codon is present.
    """
    stops = ["TAG", "TAA", "TGA"]
    stop = False
    for x in range(0, len(sequence), 3):
        if sequence[x:x+3].upper() in stops:
            stop = True
            break
    return stop


def choose_base(seq_kmer):
    """Choose a random base that is different to current centre base of kmer of
    size 5.

    Args:
        seq_kmer (str): 5 base kmer.

    Returns:
        base (str): Randomly chosen base ("A", "C", "G" or "T").
    """
    if seq_kmer[2].upper() == "A":
        base = random.choice(["C", "G", "T"])
    if seq_kmer[2].upper() == "C":
        base = random.choice(["A", "G", "T"])
    if seq_kmer[2].upper() == "G":
        base = random.choice(["A", "C", "T"])
    if seq_kmer[2].upper() == "T":
        base = random.choice(["A", "C", "G"])

    return base


def float_range(min, max):
    """ Used with argparse to check that value is within a float range - 
    min <= arg <= max

    Args:
        min (float): Minimum float value.
        max (float): Maximum float value.
    """

    def check_float_range(arg):
        """New Type function for argparse - a float within predefined range."""
        try:
            f = float(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("Expect a float value")
        if f < min or f > max:
            raise argparse.ArgumentTypeError(
                f"Must be in range {str(min)} -  {str(max)}. Mutation rates above {str(max)} are unachievable.")
        return f

    return check_float_range


def int_range(min, max):
    """ Used with argparse to check that value is within a float range - 
    min <= arg <= max

    Args:
        min (int): Minimum integer value.
        max (int): Maximum integer value.
    """

    def check_int_range(arg):
        """New Type function for argparse - a integer within predefined range."""
        try:
            f = int(arg)
        except ValueError:
            raise argparse.ArgumentTypeError("Expect an integer value")
        if f < min or f > max:
            raise argparse.ArgumentTypeError(
                f"Must be in range {str(min)} -  {str(max)}. Mutation rates above {str(max)} are unachievable.")
        return f

    # Return function handle
    return check_int_range

# Data processing functions


def parse_fasta(f):
    """Parses a FASTA format file.

    Args:
        f (file): A FASTA formatted file, each sequence is expected to have a
            header line starting with ">".

    Yields:
        generator : A generator object that produces each header and sequence
    """

    header, seq = None, []
    for line in f:
        line = line.rstrip("\n")  # remove newline character
        if line.startswith(">"):  # if line is a FASTA header line
            if header:
                yield (header, "".join(seq))  # yield header and sequence
            header, seq = line, []  # set header to contents of line and blank sequence
        else:
            seq.append(line)  # if sequence line, append contents to sequence
    if header:
        yield (header, "".join(seq))  # yield header and sequence


def create_allele_dict(fasta):
    """Creates an allele dictionary from an IMGT formatted reference FASTA.

    Args:
        fasta (file): A FASTA file consisting of reference Ig V, D or J alleles,
            with IMGT-gapped sequences. Header is expected to follow IMGT format.

    Returns:
        allele_dict (dict): Allele dictionary with V, D or J gene as keys and 
            lists of the corresponding alleles, in the form of Allele class 
            instances, as values.
            Only returns full length alleles denoted as functional or open reading
            frame by the IMGT (marked by "F" or "ORF" in the sequence header).
            {Gene : [Allele, Allele ...]}
    """

    # TODO - add exception handling here to check expected allele format etc

    with open(fasta) as f:
        allele_dict = defaultdict(list)
        for header, seq in parse_fasta(f):
            allele = header.split("|")[1]
            if allele == "IGHV3-30-3*03" or allele == "IGHV5-10-1*02":
                # IGHV3-30-3*03 has identical sequence to IGHV3-30*04 and IGHV5-10-1*02 has anchor Cys that makes it out of frame
                continue
            if "-" in allele:
                family = allele.split("-")[0]
            else:
                family = allele.split("*")[0]
            if "IGHD" in family:
                seq = seq.replace(".", "").lower()
            if "IGHV" in family:
                cys = seq[309:312]
                if cys != "tgt" and cys != "tgc":  # if allele doesn't have expected V anchor then don't use it
                    continue
            if "IGHJ" in family:  # if allele doesn't have expected J anchor then don't use it
                motif = re.compile(
                    '(ttt|ttc|tgg)(ggt|ggc|gga|ggg)[a-z]{3}(ggt|ggc|gga|ggg)')
                match = motif.search(seq)
                if match is None:
                    continue

            gene = allele.split("*")[0]
            coding = header.split("|")[3]
            length = header.split("|")[6].split(" ")[0]
            if "partial" not in header and (coding == "F" or coding == "ORF"):
                allele_dict[gene].append(Allele(allele, seq, length))

    return dict(allele_dict)


def create_family_use_dict(usage_csv):
    """Creates a dictionary for V, D or J family usage.

    Args:
        usage_csv (file): A CSV file consisting of gene family and proportion
            of sequences to use alleles from that gene family.

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        family_dict (dict): Dictionary of gene families and the proportion of 
            sequences to use their alleles. In the format {gene family: proportion}
    """
    with open(usage_csv, "r") as f:
        family_dict = {}
        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 2:
                raise ValueError(
                    f"Invalid csv format. Expected line length 2. Line {csvlines.line_num} has length {len(line)}.")
            else:
                family, prop = line
                family_dict[family] = float(prop)
        return dict(family_dict)


def create_trimming_dict(trim_csv):
    """Creates a dictionary of V,D or J trimming proportions per gene family.

    Args:
        trim_csv (file): A CSV file consisting of IGH gene family (e.g. IGHV1),
            number of trimmed nucleotides, and proportion of repertoire trimmed in
            such a way.

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        trim_dict (dict): A nested dictionary containing trimming distributions for each gene
            family. In the format {Gene family: {Trim length : proportion of sequences}}.
    """
    with open(trim_csv, "r") as f:

        trim_dict = defaultdict(lambda: defaultdict(float))

        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 3:
                raise ValueError(
                    f"Invalid csv format. Expected line length 3. Line {csvlines.line_num} has length {len(line)}.")
            else:
                family, trim, prop = line  # ! Assumes column order in file
                trim_dict[family][int(trim)] = float(prop)
        return dict(trim_dict)


def create_NP_length_dict(lengths_csv):
    """Creates a dictionary that contains distributions of NP region lengths.

    Args:
        lengths_csv (file): A CSV file consisting of NP region length and the
            proportion of sequences to use this length.

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        lengths_dict (dict): Dictionary of possible NP region lengths and the proportion of
            sequences to use them. In the format {NP region length: proportion}
    """
    with open(lengths_csv, "r") as f:
        lengths_dict = {}
        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 2:
                raise ValueError(
                    f"Invalid csv format. Expected line length 2. Line {csvlines.line_num} has length {len(line)}.")
            else:
                length, prop = line
                lengths_dict[int(float(length))] = float(prop)
        return dict(lengths_dict)


def create_first_base_dict(first_base_csv):
    """Creates a dictionary for choosing the first base of an NP region.

    Args:
        first_base_csv (file): A CSV file consisting of nucleotide base 
            (A,C,G,T) and the probability of NP regions using this nucleotide 
            as the first base.

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        first_base_dict (dict): Dictionary of nucleotide bases and the 
            probability of NP regions using this nucleotide as the first base. 
            In the format {Nucleotide base: probability, ...}
    """

    with open(first_base_csv, "r") as f:
        first_base_dict = {}
        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 2:
                raise ValueError(
                    f"Invalid csv format. Expected line length 2. Line {csvlines.line_num} has length {len(line)}.")
            else:
                base, prop = line
                first_base_dict[base] = float(prop)
        return dict(first_base_dict)


def create_NP_position_transition_dict(transitions_csv):
    """Creates a dictionary representing a transition matrix for NP region.

    Creates a dictionary representing a transition matrix with the 
    probabilities of transitioning from one nucleotide base (A,C,G,T) to another
    within the NP region based on the current position in the sequence.

    Args:
        transitions_csv (file): A CSV file consisting of nucleotide base 
            (A,C,G,T) and the probability of transitioning to each other base for
            each position in the NP region

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        NP_len_transition_dict (dict): Dictionary of nucleotide bases and 
            the probability of transitioning to each other base at each position 
            in NP region. In the format 
            {Position: {Nucleotide base: 
                    {Nucleotide base: probability, ...}, ...}, ...}
    """

    with open(transitions_csv, "r") as f:

        NP_len_transition_dict = defaultdict(lambda: defaultdict(dict))

        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 6:
                raise ValueError(
                    f"Invalid csv format. Expected line length 6. Line {csvlines.line_num} has length {len(line)}.")
            else:
                position, base, A, C, G, T = line
                NP_len_transition_dict[int(float(position))][base] = {
                    "A": float(A),
                    "C": float(C),
                    "G": float(G),
                    "T": float(T),
                }
        return dict(NP_len_transition_dict)


def create_kmer_base_dict(kmer_csv):
    """
    Creates a dictionary of kmer mutation probabilities.

    Args:
        kmer_csv (file): A CSV file consisting of unique kmer and the 
            proportion of sequences using each nucleotide base (A, T, C, G)
            at the centre position of this kmer.

    Returns:
        kmer_dict (dict): Nested dictionary containing nucleotide distributions for the 
            centre position of each kmer. In the format
            {kmer: {A : Proportion, T : Proportion,
            C : Proportion, G : Proportion }}
    """

    with open(kmer_csv, "r") as f:
        kmer_dict = defaultdict(lambda: defaultdict(dict))

        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 5:
                raise ValueError(
                    f"Invalid csv format. Expected line length 5. Line {csvlines.line_num} has length {len(line)}.")
            else:
                kmer, A, C, G, T = line
                kmer_dict[kmer.lower()] = {
                    "A": float(A),
                    "C": float(C),
                    "G": float(G),
                    "T": float(T),
                }
        return dict(kmer_dict)


def load_data(data_folder=None, mutate=False):
    """Loads and processes required data files from data folder.

    Args:
        data_folder (path, optional): Path to data folder with required data.
            If not specified then uses inbuilt package data. Defaults to None.
        mutate (bool, optional): Whether to read in data for mutated sequences 
            or not. Defaults to False.

    Returns:
        data_dict (dict): Dictionary containing all required data for generating 
            sequences. Includes family_use_dict, gene_use_dict, trim_dicts,
            NP_transitions, NP_first_bases, NP_lengths, mut_rate_per_seq and kmer_dicts.
    """

    data_dict = {}

    if data_folder != None:
        path_to_data = data_folder

        # GENE USAGE

        v_usage = create_family_use_dict(
            f"{path_to_data}/IGHV_usage.csv")
        d_usage = create_family_use_dict(
            f"{path_to_data}/IGHD_usage.csv")
        j_usage = create_family_use_dict(
            f"{path_to_data}/IGHJ_usage.csv")

        data_dict["family_use_dict"] = {"V": v_usage,
                                        "D": d_usage,
                                        "J": j_usage}

        v_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHV_usage_gene.csv")

        d_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHD_usage_gene.csv")

        j_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHJ_usage_gene.csv")

        data_dict["gene_use_dict"] = {"V": v_gene_usage,
                                      "D": d_gene_usage,
                                      "J": j_gene_usage}

        v_trim = create_trimming_dict(
            f"{path_to_data}/V_family_trimming_proportions.csv")
        d_3_trim = create_trimming_dict(
            f"{path_to_data}/D_3_family_trimming_proportions.csv")
        d_5_trim = create_trimming_dict(
            f"{path_to_data}/D_5_family_trimming_proportions.csv")
        j_trim = create_trimming_dict(
            f"{path_to_data}/J_family_trimming_proportions.csv")

        data_dict["trim_dicts"] = {"V_3": v_trim, "D_5": d_5_trim,
                                   "D_3": d_3_trim, "J_5": j_trim}

        # N REGIONS

        # LENGTHS
        NP1_lengths = create_NP_length_dict(
            f"{path_to_data}/np1_lengths_proportions.csv")
        NP2_lengths = create_NP_length_dict(
            f"{path_to_data}/np2_lengths_proportions.csv")

        # FIRST BASE USAGE
        NP1_first_base_use = create_first_base_dict(
            f"{path_to_data}/np1_first_base_probs.csv"
        )
        NP2_first_base_use = create_first_base_dict(
            f"{path_to_data}/np2_first_base_probs.csv"
        )

        # TRANSITION MATRICES
        if mutate == False:
            NP1_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np1_transition_probs_per_position_igdm.csv")

            NP2_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np2_transition_probs_per_position_igdm.csv")

        if mutate == True:
            NP1_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np1_transition_probs_per_position_igag.csv")

            NP2_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np2_transition_probs_per_position_igag.csv")

        data_dict["NP_transitions"] = {"NP1": NP1_transitions,
                                       "NP2": NP2_transitions}
        data_dict["NP_first_bases"] = {"NP1": NP1_first_base_use,
                                       "NP2": NP2_first_base_use}
        data_dict["NP_lengths"] = {"NP1": NP1_lengths,
                                   "NP2": NP2_lengths}

        with open(f"{path_to_data}/mut_freq_per_seq_per_family.csv", newline='') as f:
            mut_rate_per_seq = defaultdict(dict)
            csvlines = csv.reader(f, delimiter=",")
            next(csvlines)  # skip header
            for line in csvlines:
                if not len(line) == 3:
                    raise ValueError(
                        f"Invalid csv format. Expected line length 3. Line {csvlines.line_num} has length {len(line)}.")
                else:
                    family, mut_rate, prop = line
                    mut_rate_per_seq[family][float(mut_rate)] = float(prop)

        data_dict["mut_rate_per_seq"] = mut_rate_per_seq

        # PER BASE MUTATION

        cdr1_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr1_kmer_base_usage.csv")
        cdr2_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr2_kmer_base_usage.csv")
        cdr3_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr3_kmer_base_usage.csv")
        cdr_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr_kmer_base_usage.csv")
        fwr1_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr1_kmer_base_usage.csv")
        fwr2_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr2_kmer_base_usage.csv")
        fwr3_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr3_kmer_base_usage.csv")
        fwr4_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr4_kmer_base_usage.csv")
        fwr_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr_kmer_base_usage.csv")

        data_dict["kmer_dicts"] = {"fwr1": fwr1_kmers, "fwr2": fwr2_kmers, "fwr3": fwr3_kmers, "fwr4": fwr4_kmers,
                                   "cdr1": cdr1_kmers, "cdr2": cdr2_kmers, "cdr3": cdr3_kmers, "cdr": cdr_kmers, "fwr": fwr_kmers}

    else:
        try:
            # account for differing versions of importlib.resources
            path_to_data = importlib.resources.files(
                'airrship').joinpath("data")
        except AttributeError:
            with importlib.resources.path('airrship', 'data') as p:
                path_to_data = p
            # GENE USAGE

        v_usage = create_family_use_dict(
            f"{path_to_data}/IGHV_usage.csv")
        d_usage = create_family_use_dict(
            f"{path_to_data}/IGHD_usage.csv")
        j_usage = create_family_use_dict(
            f"{path_to_data}/IGHJ_usage.csv")

        data_dict["family_use_dict"] = {"V": v_usage,
                                        "D": d_usage,
                                        "J": j_usage}

        v_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHV_usage_gene.csv")

        d_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHD_usage_gene.csv")

        j_gene_usage = create_family_use_dict(
            f"{path_to_data}/IGHJ_usage_gene.csv")

        data_dict["gene_use_dict"] = {"V": v_gene_usage,
                                      "D": d_gene_usage,
                                      "J": j_gene_usage}

        v_trim = create_trimming_dict(
            f"{path_to_data}/V_family_trimming_proportions.csv")
        d_3_trim = create_trimming_dict(
            f"{path_to_data}/D_3_family_trimming_proportions.csv")
        d_5_trim = create_trimming_dict(
            f"{path_to_data}/D_5_family_trimming_proportions.csv")
        j_trim = create_trimming_dict(
            f"{path_to_data}/J_family_trimming_proportions.csv")

        data_dict["trim_dicts"] = {"V_3": v_trim, "D_5": d_5_trim,
                                   "D_3": d_3_trim, "J_5": j_trim}

        # N REGIONS

        # LENGTHS
        NP1_lengths = create_NP_length_dict(
            f"{path_to_data}/np1_lengths_proportions.csv")
        NP2_lengths = create_NP_length_dict(
            f"{path_to_data}/np2_lengths_proportions.csv")

        # FIRST BASE USAGE
        NP1_first_base_use = create_first_base_dict(
            f"{path_to_data}/np1_first_base_probs.csv"
        )
        NP2_first_base_use = create_first_base_dict(
            f"{path_to_data}/np2_first_base_probs.csv"
        )

        # TRANSITION MATRICES
        if mutate == False:
            NP1_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np1_transition_probs_per_position_igdm.csv")

            NP2_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np2_transition_probs_per_position_igdm.csv")

        if mutate == True:
            NP1_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np1_transition_probs_per_position_igag.csv")

            NP2_transitions = create_NP_position_transition_dict(
                f"{path_to_data}/np2_transition_probs_per_position_igag.csv")

        data_dict["NP_transitions"] = {"NP1": NP1_transitions,
                                       "NP2": NP2_transitions}
        data_dict["NP_first_bases"] = {"NP1": NP1_first_base_use,
                                       "NP2": NP2_first_base_use}
        data_dict["NP_lengths"] = {"NP1": NP1_lengths,
                                   "NP2": NP2_lengths}

        # PER SEQUENCE MUTATION

        with open(f"{path_to_data}/mut_freq_per_seq_per_family.csv") as f:
            mut_rate_per_seq = defaultdict(dict)
            csvlines = csv.reader(f, delimiter=",")
            next(csvlines)  # skip header
            for line in csvlines:
                if not len(line) == 3:
                    raise ValueError(
                        f"Invalid csv format. Expected line length 3. Line {csvlines.line_num} has length {len(line)}.")
                else:
                    family, mut_rate, prop = line
                    mut_rate_per_seq[family][float(mut_rate)] = float(prop)

        data_dict["mut_rate_per_seq"] = mut_rate_per_seq

        cdr1_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr1_kmer_base_usage.csv")
        cdr2_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr2_kmer_base_usage.csv")
        cdr3_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr3_kmer_base_usage.csv")
        cdr_kmers = create_kmer_base_dict(
            f"{path_to_data}/cdr_kmer_base_usage.csv")
        fwr1_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr1_kmer_base_usage.csv")
        fwr2_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr2_kmer_base_usage.csv")
        fwr3_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr3_kmer_base_usage.csv")
        fwr4_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr4_kmer_base_usage.csv")
        fwr_kmers = create_kmer_base_dict(
            f"{path_to_data}/fwr_kmer_base_usage.csv")

        data_dict["kmer_dicts"] = {"fwr1": fwr1_kmers, "fwr2": fwr2_kmers, "fwr3": fwr3_kmers, "fwr4": fwr4_kmers,
                                   "cdr1": cdr1_kmers, "cdr2": cdr2_kmers, "cdr3": cdr3_kmers, "cdr": cdr_kmers, "fwr": fwr_kmers}

    return data_dict


# Sequence creation functions

def create_locus_segment(segment, allele_dict, hetero_prop=0):
    """Chooses alleles present in gene segment of locus.

    Randomly selects the V, D or J alleles to be included on each chromosome.
    The proportion of positions at which an individual is heterozygous (i.e.
    where each chromosome carries a different allele can be determined)

    Args:
        segment (str): The V, D, or J segment for which alleles are to be chosen.
        allele_dict (dict): Allele dictionary with V/D/J gene as keys and lists of the
            corresponding alleles, in the form of Allele class instances, as values.
            {Gene : [Allele, Allele ...], ...}
        hetero_prop (int, optional): Proportion of gene positions for which each
            chromosome carries a different allele. Defaults to 0.

    Returns:
        locus_segment (list): Returns list of two lists of Alleles, each representing the
            V, D or J alleles present on one chromosome.
    """

    chromosome1_segment, chromosome2_segment = [], []

    # get number of genes which will be heterozygous
    ref_het = int(round(hetero_prop * len(allele_dict)))
    het = int(round(hetero_prop * len(allele_dict)))
    # randomly sort allele dictionary values
    for i, alleles in enumerate(
        sorted(allele_dict.values(), key=lambda x: random.random())
    ):
        allele = random.choice(alleles)  # choose first allele
        chromosome1_segment.append(allele)
        # if only one allele of that gene exists, has to be homozygous
        if len(alleles) > 1:
            if het > 0:  # if count for heterozygous genes still > 0
                # remove the chr1 chosen allele from list
                alleles.remove(allele)
                # re-choose to get different allele
                het_allele = random.choice(alleles)
                chromosome2_segment.append(het_allele)
                het -= 1  # decrease counter to account for new heterozygous loci
            else:
                chromosome2_segment.append(allele)
        else:
            chromosome2_segment.append(allele)
        if i == len(allele_dict) - 1:
            if het > 0:
                print(
                    f"Requested proportion heterozygous genes: {hetero_prop} for {segment} genes, maximum proportion possible given: {round((ref_het-het)/len(allele_dict), 2)}"
                )
    locus_segment = [chromosome1_segment, chromosome2_segment]
    return locus_segment


def create_locus(segments, allele_dicts, hetero_props, haplotype=True):
    """Creates a two chromosome gene locus containing V, D and J alleles.

    Args:
        segments (str): Which segments to include, expect V, D and J.
            allele_dicts (dict): Dictionary of dictionaries, in the format of
            {Segment : {Gene : [Allele, Allele ...]}}
        hetero_props (list): List of integer proportion of positions to be
            heterozygous for each position

    Returns:
        locus (list): List of two dictionaries. Each is a dictionary containing the gene
            segment as keys and the chosen alleles as values. Format is
            {Segment : [Allele, Allele ...], ...}
    """
    chromosome1, chromosome2 = defaultdict(list), defaultdict(list)

    if haplotype == True:
        for segment, hetero_prop in zip(segments, hetero_props):
            chromosome1[segment], chromosome2[segment] = create_locus_segment(
                segment, allele_dicts[segment], hetero_prop
            )
    else:
        for segment in segments:
            allele_dict = allele_dicts[segment]
            for gene in allele_dict.values():
                for allele in gene:
                    chromosome1[segment].append(allele)
                    chromosome2[segment].append(allele)

    locus = [chromosome1, chromosome2]
    return locus


def get_locus_from_file(locus_csv, allele_dicts):
    """Return a a two chromosome gene locus defined by a file.

    Args:
        locus_csv (file): CSV file containing alleles present.
        allele_dicts (dict): Nested dictionary of dictionaries. Outer dictionary
            has V,D,J as keys and dictionary, with V, D or J gene as keys and lists 
            of the corresponding alleles, in the form of Allele class instances, 
            as values, as values.
            {"V" : {Gene : [Allele, Allele ...]}, ...}

    Raises:
        ValueError: Warns of unexpected file format - checks line length.

    Returns:
        locus (list): List of two dictionaries. Each is a dictionary containing 
            the gene segment as keys and the chosen alleles as values. Format is
            {Segment : [Allele, Allele ...], ...}
    """

    with open(locus_csv, "r") as f:

        chromosome1, chromosome2 = defaultdict(list), defaultdict(list)

        csvlines = csv.reader(f, delimiter=",")
        next(csvlines)  # skip header
        for line in csvlines:
            if not len(line) == 2:
                raise ValueError(
                    f"Invalid csv format. Expected line length 2. Line {csvlines.line_num} has length {len(line)}.")
            else:
                allele_1_name, allele_2_name = line
                segment = allele_1_name.split("IGH", 1)[1][0]
                gene = allele_1_name.split("*")[0]
                for allele in allele_dicts[segment][gene]:
                    if allele.name == allele_1_name:
                        chromosome1[segment].append(allele)
                segment = allele_2_name.split("IGH", 1)[1][0]
                gene = allele_2_name.split("*")[0]
                for allele in allele_dicts[segment][gene]:
                    if allele.name == allele_2_name:
                        chromosome2[segment].append(allele)

        locus = [chromosome1, chromosome2]
        return locus


def write_locus_file(out_name, locus):
    """Write out locus to file.

    Args:
        out_name (str): Prefix to add to filename.
        locus (list): List of two dictionaries. Each is a dictionary containing 
            the gene segment as keys and the chosen alleles as values. Format is
            {Segment : [Allele, Allele ...], ...}
    """
    save_locus = {}
    for i, chromosome in enumerate(locus):
        save_locus[i] = []
        for segment in chromosome.values():
            for allele in segment:
                save_locus[i].append(allele.name)

    with open(out_name, 'w') as f:
        writer = csv.writer(f, delimiter=",", lineterminator="\n")
        writer.writerow(
            ["Chr 1"]
            + ["Chr 2"])
        for allele_1, allele_2 in zip(save_locus[0], save_locus[1]):
            row = [allele_1] + [allele_2]
            writer.writerow(row)


def write_summary_file(out_name, args):
    """Write out summary of the arguments passed to commandline call.

    Args:
        out_name (str): Prefix to add to filename.
        args (Namespace): Output from argparse.ArgumentParser.parse_args(), 
            commandline arguments specified in run.
    """
    with open(out_name, 'w') as f:
        f.write("AIRRSHIP was run with the following arguments: \n")
        output = str(args).rstrip(")").lstrip("Namespace(")
        f.write(str(output))


def get_genotype(data_folder=None, het_list=[1, 1, 1], haplotype=True, locus=None):
    """Wrapper that generates a locus for use in sequence generations

    Args:
        data_folder (path, optional): Path to data folder with required data. 
            When not specified will use package data. Defaults to None.
        het_list (list, optional): Proportion of genes [V, D, J] to be heterozygous.
            Defaults to [1, 1, 1].
        haplotype (bool, optional): True when only two alleles per gene are 
            to be used. Defaults to True.
        locus (path, optional): Path to file with predefined locus. 
            Defaults to None.

    Returns:
        locus (list): List of two dictionaries. Each is a dictionary containing 
            the gene segment as keys and the chosen alleles as values. Format is
            {Segment : [Allele, Allele ...], ...}
    """
    if data_folder != None:
        path_to_data = data_folder

        v_alleles = create_allele_dict(f"{path_to_data}/imgt_human_IGHV.fasta")
        d_alleles = create_allele_dict(f"{path_to_data}/imgt_human_IGHD.fasta")
        j_alleles = create_allele_dict(f"{path_to_data}/imgt_human_IGHJ.fasta")

    else:
        try:
            path_to_data = importlib.resources.files(
                'airrship').joinpath("data")
        except AttributeError:
            with importlib.resources.path('airrship', 'data') as p:
                path_to_data = p
        v_alleles = create_allele_dict(
            f"{path_to_data}/imgt_human_IGHV.fasta")
        d_alleles = create_allele_dict(
            f"{path_to_data}/imgt_human_IGHD.fasta")
        j_alleles = create_allele_dict(
            f"{path_to_data}/imgt_human_IGHJ.fasta")

    vdj_allele_dicts = {"V": v_alleles,
                        "D": d_alleles,
                        "J": j_alleles}

    if locus == None:
        locus = create_locus(
            ["V", "D", "J"], vdj_allele_dicts, het_list, haplotype)
    else:
        locus = get_locus_from_file(locus, vdj_allele_dicts)

    return locus


def random_sequence(locus, gene_use_dict, family_use_dict, segments=["V", "D", "J"], flat=False):
    """Creates a random sequence.

    Initialises a Sequence object from randomly chosen alleles from a previously 
    created locus. The gene to pick alleles from is randomly chose based on the
    desired gene usage distribution.

    Args:
        locus (list): List of two dictionaries. Each is a dictionary containing 
            the gene segment as keys and the chosen alleles as values. Format is
            [{Segment : [Allele, Allele ...], ...}, {Segment : [Allele, Allele ...], ...}]
        gene_use_dict (_type_): Nested dictionary of gene segment and
            genes and the proportion of sequences to use their alleles. 
            In the format {segment: {gene: proportion, ...}, ...}
        family_use_dict (_type_): Nested dictionary of gene segment and
            gene families and the proportion of sequences to use their alleles. 
            In the format {segment: {gene family: proportion, ...}, ...}
        segments (list, optional): Gene segments to include in sequence. 
            Defaults to ["V", "D", "J"].
        flat (optional): gene, family or False. Gene or family specify that 
            sequences should use all genes or gene families evenly. If false, 
            usage follows experimental distributions. Defaults to False.

    Returns:
        out_sequence (Sequence): Randomly generated Sequence class object.
    """

    alleles = []
    chromosome = random.choice(locus)
    for segment in segments:
        if flat == False:  # if family usage is to be biased
            counter = 1
            while counter > 0:
                gene_dict = gene_use_dict[segment]
                genes, probs = zip(*gene_dict.items())
                gene = choice(genes, probs)
                gene_alleles = [
                    x for x in chromosome[segment] if x.gene == gene]
                if len(gene_alleles) > 0:
                    allele = random.choice(gene_alleles)
                    counter = - 1
                    alleles.append(allele)
        elif flat == "gene":
            counter = 1
            while counter > 0:
                gene_dict = gene_use_dict[segment]
                gene = random.choice(list(gene_dict.keys()))
                gene_alleles = [
                    x for x in chromosome[segment] if x.gene == gene]
                if len(gene_alleles) > 0:
                    allele = random.choice(gene_alleles)
                    counter = - 1
                    alleles.append(allele)
        elif flat == "family":  # if family usage is to be even
            family_dict = family_use_dict[segment]
            families = list(family_dict.keys())
            family = random.choice(families)
            allele = random.choice(
                [x for x in chromosome[segment] if x.family == family])
            alleles.append(allele)
    out_sequence = Sequence(alleles[0], alleles[1], alleles[2])
    return out_sequence


def get_NP_regions_markov_pos(NP_lengths, NP_transitions, which_NP, first_base_dict):
    """Creates a NP region using a first order markov chain in which the 
    transition matrix is chosen based on the current position in the NP region.

    Args:
        NP_lengths (dict): Dictionary of possible NP region lengths and the
            proportion of sequences to use them. In the format
            {NP region length: proportion}
        NP_transitions (dict): Nested dictionary containing transition matrix of
            probabilities of moving from one nucleotide (A, C, G, T) to any other 
            for each position in the NP region.
        which_NP (str): Which NP region to create - NP1 or NP2
        first_base_dict (dict): Dictionary of the proportion of sequences starting
            with each base.

    Returns:
        NP_seq (str): An NP region of chosen length. Uppercase.
    """

    NP_seq = ""
    length_dict = NP_lengths[which_NP]
    lengths, probs = zip(*length_dict.items())
    length = choice(lengths, probs)
    if length != 0:
        nucs, probs = zip(*first_base_dict[which_NP].items())
        first_base = choice(nucs, probs)
        np_region = NP_Region(NP_transitions[which_NP], first_base, length)
        NP_seq = np_region.generate_np_seq()
    return NP_seq


def get_trimmed_seq(allele, no_trim_list, trim_dicts):
    """Generate trimmed ungapped allele nucleotide sequence.

    Args:
        allele (Allele): Allele to be trimmed.
        no_trim_list (list): List of 5 Booleans, specifying whether to not
            trim [all_ends, v_3_end, d_5_end, d_3_end, j_5_end].         
        trim_dicts (dict): A dictionary of dictionaries of trimming length 
            proportions by gene family for each segment (V, D or J).

    Returns:
        trimmed_seq (str): Trimmed nucleotide sequence, lower case, ungapped.
    """

    allele.get_trim_length(no_trim_list, trim_dicts)
    trim_5 = allele.trim_5
    trim_3 = allele.trim_3
    trimmed_seq = ""
    if allele.segment == "D":
        if trim_5 > 0 and trim_3 > 0:
            trimmed_seq = allele.ungapped_seq[trim_5:-trim_3]
        elif trim_5 > 0 and trim_3 == 0:
            trimmed_seq = allele.ungapped_seq[trim_5:]
        elif trim_5 == 0 and trim_3 > 0:
            trimmed_seq = allele.ungapped_seq[:-trim_3]
        else:
            trimmed_seq = allele.ungapped_seq
    elif allele.segment == "V":
        if trim_3 > 0:
            trimmed_seq = allele.ungapped_seq[:-trim_3]
        else:
            trimmed_seq = allele.ungapped_seq
    elif allele.segment == "J":
        if trim_5 > 0:
            trimmed_seq = allele.ungapped_seq[trim_5:]
        else:
            trimmed_seq = allele.ungapped_seq
    else:
        print("error in gene name")
    return trimmed_seq


def get_gapped_trimmed_seq(allele, no_trim_list, trim_dicts):
    """Generate trimmed ungapped allele nucleotide sequence.

    Args:
        allele (Allele): Allele to be trimmed.
        no_trim_list (list): List of 5 Booleans, specifying whether to not
            trim [all_ends, v_3_end, d_5_end, d_3_end, j_5_end].         
        trim_dicts (dict): A dictionary of dictionaries of trimming length 
            proportions by gene family for each segment (V, D or J).

    Returns:
        trimmed_seq (str): Trimmed nucleotide sequence, lower case, gapped.
    """

    allele.get_trim_length(no_trim_list, trim_dicts)
    trim_5 = allele.trim_5
    trim_3 = allele.trim_3
    trimmed_seq = ""
    if allele.segment == "D":
        if trim_5 > 0 and trim_3 > 0:
            trimmed_seq = allele.gapped_seq[trim_5:-trim_3]
        elif trim_5 > 0 and trim_3 == 0:
            trimmed_seq = allele.gapped_seq[trim_5:]
        elif trim_5 == 0 and trim_3 > 0:
            trimmed_seq = allele.gapped_seq[:-trim_3]
        else:
            trimmed_seq = allele.gapped_seq
    elif allele.segment == "V":
        if trim_3 > 0:
            trimmed_seq = allele.gapped_seq[:-trim_3]
        else:
            trimmed_seq = allele.gapped_seq
    elif allele.segment == "J":
        if trim_5 > 0:
            trimmed_seq = allele.gapped_seq[trim_5:]
        else:
            trimmed_seq = allele.gapped_seq
    else:
        print("error in gene name")
    return trimmed_seq


def mutate_by_kmer_random_region(sequence, v_family, junction_length, kmer_mut_info, mut_per_seq, mut_rate=None, number_muts=None):
    """Adds hypermutation to a sequence.

    Hypermutates sequence based on proportion of mutation at central position
    of each 5-mer in sequence, varies by region of sequence.

    Args:
        sequence (str): Sequence to be mutated. Expect gapped, lowercase str
            with NP regions in uppercase.
        junction_length (int): Number of nucleotides in junction (CDR3 + anchors)
        kmer_mut_info (dict): Dict of dicts, with the output for create_kmer_base_dict
            for each IMGT region (FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, FWR4) and overall
            FWR and CDR measures. {region: {kmer: {A : Proportion, T : Proportion,
            C : Proportion, G : Proportion }, ...}, ...}
        mut_per_seq (dict): Mutation frequencies and probabilities for their use.
        mut_rate (float, optional): Mutation rate to be used rather than choosing
            from distribution. Defaults to None.
        number_muts (int, optional): Number of mutations to be added rather than
            choosing from distribution. Defaults to None.

    Returns:
        mut_seq (str): Mutated sequence. Ungapped, all uppercase.
        mutations (dict): Mutation events, {position: "old_base>new_base", ...}.
        mut_rate (int): Mutation rate used. 

    """

    sequence = sequence
    mut_seq = sequence[:2]
    ungapped_seq = sequence.replace(".", "")
    mutations = {}
    v_family = v_family

    if mut_rate == None and number_muts == None:
        mut_dict = mut_per_seq[v_family]
        rates, probs = zip(*mut_dict.items())
        mut_rate = choice(rates, probs)
        # mut_rate = random.choice(mut_per_seq)
        number_muts = int(round(len(sequence) * mut_rate))
    elif number_muts == None:
        number_muts = int(round(len(sequence) * mut_rate))
    else:
        number_muts = number_muts
    if mut_rate == 0 or number_muts == 0:
        return ungapped_seq, mutations, mut_rate

    unmut_kmers = defaultdict()
    final_bases = defaultdict()

    fwr1 = sequence[:78]  # this has to be the gapped sequence
    cdr1 = sequence[78:114]
    fwr2 = sequence[114:165]
    cdr2 = sequence[165:195]
    fwr3 = sequence[195:312]
    cdr3 = sequence[312: 312 + junction_length - 6]
    # fwr4 = sequence[312 + junction_length - 6: len(sequence)]

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

    for x in range(0, len(ungapped_seq) - 4):
        seq_kmer = ungapped_seq[x:x+5]
        unmut_kmers[x+2] = seq_kmer

    items = list(unmut_kmers.items())  # List of tuples
    random.shuffle(items)  # random order

    # whilst still need to mutate bases and bases left to be mutated
    while number_muts > 0 and len(items) > 0:
        for item in items:  # loop through position, kmer tuples
            position = item[0]
            kmer = item[1]
            if number_muts > 0:  # double check that still want to mutate
                if cdr1_start <= position < cdr1_end:
                    region = 'cdr1'
                elif fwr2_end <= position < cdr2_end:
                    region = 'cdr2'
                elif fwr3_end <= position < cdr3_end:
                    region = 'cdr3'
                elif position < cdr1_start:
                    region = 'fwr1'
                elif cdr1_end <= position < fwr2_end:
                    region = 'fwr2'
                elif cdr2_end <= position < fwr3_end:
                    region = 'fwr3'
                else:
                    region = 'fwr4'
                # ignore N regions, reliant on them being added in upper case (alternative way?)
                if kmer[2] in ("A", "C", "G", "T"):
                    base = kmer[2]
                    final_bases[position] = base.upper()
                    items.remove(item)
                elif kmer in kmer_mut_info[region]:
                    prob_dict = kmer_mut_info[region][kmer]
                    # possible that this base is never mutated in data
                    if prob_dict[kmer[2].upper()] == 1:
                        base = kmer[2]
                        final_bases[position] = base.upper()
                        items.remove(item)
                    nucs, probs = zip(*prob_dict.items())
                    base = choice(nucs, probs)
                    final_bases[position] = base.lower()
                    if base.lower() != kmer[2]:
                        number_muts -= 1
                        items.remove(item)
                        mutations[position+1] = f"{kmer[2].upper()}>{base}"
                elif kmer in kmer_mut_info[region[:-1]]:
                    prob_dict = kmer_mut_info[region[:-1]][kmer]
                    if prob_dict[kmer[2].upper()] == 1:
                        base = kmer[2]
                        final_bases[position] = base.upper()
                        items.remove(item)
                    nucs, probs = zip(*prob_dict.items())
                    base = choice(nucs, probs)
                    final_bases[position] = base.lower()
                    if base.lower() != kmer[2]:
                        number_muts -= 1
                        items.remove(item)
                        mutations[position+1] = f"{kmer[2].upper()}>{base}"
                else:
                    base = kmer[2]
                    final_bases[position] = base.upper()
                    items.remove(item)
            else:
                base = kmer[2]
                final_bases[position] = base.upper()

    for item in items:
        position = item[0]
        kmer = item[1]
        base = kmer[2]
        final_bases[position] = base.upper()

    for x in range(0, len(ungapped_seq) - 4):
        mut_seq += final_bases[x+2]

    mut_seq += ungapped_seq[-2:]
    if number_muts > 0:
        warnings.warn("failed to achieve desired mutations")
    return mut_seq, mutations, mut_rate


def mutate_randomly(sequence, mut_per_seq, mut_rate=False, number_muts=False):
    """Adds random hypermutation to a sequence.

    Hypermutates sequence randomly without taking into account kmer context of
    bases - each position has the same likelihood of mutation. 

    Args:
        sequence (str): Sequence to be mutated. Expect gapped, lowercase str
            with NP regions in uppercase.
        mut_per_seq (dict): Mutation frequencies and probabilities for their use.
        mut_rate (float, optional): Mutation rate to be used rather than choosing
            from distribution. Defaults to None.
        number_muts (int, optional): Number of mutations to be added rather than
            choosing from distribution. Defaults to None.

    Returns:
        mut_seq (str): Mutated sequence. Ungapped, all uppercase.
        mutations (dict): Mutation events, {position: "old_base>new_base", ...}.
        mut_rate (int): Mutation rate used. 
    """

    sequence = sequence.replace(".", "")

    if mut_rate == False:
        if number_muts == False:
            rates, probs = zip(*mut_per_seq.items())
            mut_rate = choice(rates, probs)
            # mut_rate = random.choice(mut_per_seq)
            number_muts = int(round(len(sequence) * mut_rate))
    else:
        number_muts = int(round(len(sequence) * mut_rate))
    unmut_bases = defaultdict()
    final_bases = defaultdict()

    for x in range(0, len(sequence)):
        base = sequence[x]
        unmut_bases[x] = base
    mut_seq = ""

    mutations = {}

    while len(final_bases) < len(sequence):
        # loop through position, kmer tuples:
        for item in sorted(unmut_bases.items(), key=lambda x: random.random()):
            position = item[0]
            base = item[1]
            if number_muts > 0:  # check that still want to mutate
                fate = random.choice(["yes", "no"])
                if fate == "yes":
                    if base == "a":
                        new_base = random.choice(["C", "G", "T"])
                    elif base == "c":
                        new_base = random.choice(["A", "G", "T"])
                    elif base == "g":
                        new_base = random.choice(["A", "C", "T"])
                    elif base == "t":
                        new_base = random.choice(["A", "C", "G"])
                    else:  # should only happen for the NP regions, don't mutate them still
                        new_base = base
                        final_bases[position] = new_base
                        del unmut_bases[position]
                    if new_base.upper() != base.upper():
                        mutations[position+1] = f"{base.upper()}>{new_base}"
                else:
                    new_base = base
                if new_base.upper() != base.upper():
                    number_muts -= 1
                    final_bases[position] = new_base.upper()
                    del unmut_bases[position]
            else:
                final_bases[position] = base
                del unmut_bases[position]

    for x in range(0, len(sequence)):
        mut_seq += final_bases[x]

    if number_muts > 0:
        print("failed to achieve desired mutations")
    return mut_seq, mutations, mut_rate


def overarch_mutation(sequence, v_family, junction_length, kmer_mut_info, shm_flat, shm_random, mut_per_seq, mut_rate=None, mut_num=None, repeat=False):
    """Wrapper function for mutating sequences.

    Args:
        sequence (str): Sequence to be mutated. Expect gapped, lowercase str
            with NP regions in uppercase.
        junction_length (int): Number of nucleotides in junction (CDR3 + anchors)
        kmer_mut_info (dict): Dict of dicts, with the output for create_kmer_base_dict
            for each IMGT region (FWR1, CDR1, FWR2, CDR2, FWR3, CDR3, FWR4) and overall
            FWR and CDR measures. {region: {kmer: {A : Proportion, T : Proportion,
            C : Proportion, G : Proportion }, ...}, ...}        
        shm_flat (bool): True if SHM is to be even across all sequences.
        shm_random (bool): True if per base mutation is to be random.
        mut_per_seq (list): List of mutation frequencies to choose from.
        mut_rate (float, optional): Mutation rate to be used rather than choosing
            from distribution. Defaults to None.        
        mut_num (int, optional): Number of mutations to be added rather than
            choosing from distribution. Defaults to None.
        repeat (bool, optional): True if a repeated attempt at SHM after 
            previous failed to generate productive sequence (ensures preservation 
            of mutation rate). Defaults to False.

    Returns:
        mut_seq (str): Mutated sequence. Ungapped, all uppercase.
        mutations (dict): Mutation events, {position: "old_base>new_base", ...}.
        mut_rate (int): Mutation rate used. 
    """

    if shm_flat == True:
        if shm_random == False:
            if mut_num is not None:
                mutated_seq, mutations, mutation_rate = mutate_by_kmer_random_region(
                    sequence, v_family, junction_length, kmer_mut_info, mut_per_seq, number_muts=mut_num)
            else:
                mutated_seq, mutations, mutation_rate = mutate_by_kmer_random_region(
                    sequence, v_family, junction_length, kmer_mut_info, mut_per_seq, mut_rate=mut_rate)
        elif shm_random == True:
            if mut_num is not None:
                mutated_seq, mutations, mutation_rate = mutate_randomly(
                    sequence, mut_per_seq, number_muts=mut_num)
            else:
                mutated_seq, mutations, mutation_rate = mutate_randomly(
                    sequence, mut_per_seq, mut_rate=mut_rate)
    elif shm_random == True:
        if repeat == False:
            mutated_seq, mutations, mutation_rate = mutate_randomly(
                sequence, mut_per_seq)
        else:
            mutated_seq, mutations, mutation_rate = mutate_randomly(
                sequence, mut_per_seq, mut_rate=mut_rate)
    else:
        if repeat == False:
            mutated_seq, mutations, mutation_rate = mutate_by_kmer_random_region(
                sequence, v_family, junction_length, kmer_mut_info, mut_per_seq)
        else:
            mutated_seq, mutations, mutation_rate = mutate_by_kmer_random_region(
                sequence, v_family, junction_length, kmer_mut_info, mut_per_seq, mut_rate=mut_rate)

    return mutated_seq, mutations, mutation_rate


def gap_seq(gapped_seq, mutated_seq):
    """Adds IMGT gaps back into a mutated sequence based on its unmutated form.

    Args:
        gapped_seq (str): Unmutated gapped sequence.
        mutated_seq (str): Mutated ungapped sequence.

    Returns:
        out_seq (str): Gapped mutated sequence.
    """

    gaps = []
    for i, base in enumerate(gapped_seq):
        if base == ".":
            gaps.append(i)
    gap_add = 0
    out_seq = ""
    for i in range(len(gapped_seq)):
        if i in gaps:
            out_seq += "."
            gap_add += 1
        else:
            out_seq += mutated_seq[i - gap_add]

    return out_seq

# Sequence generation


def generate_sequence(locus, data_dict, mutate=False, flat_usage=False, no_trim_list=(False, False, False, False, False), no_np_list=(False, False, False), shm_flat=False, shm_random=False, mutation_rate=None, mutation_number=None):
    """Wrapper to bring together entire sequence generation process.

    Recombines, trims and mutates. Will only produce functional sequences 
    (sequences with an in-frame V and J gene, no stop codons and the expected 
    junction anchor residues). 

    Args:
        locus (list): List of two dictionaries. Each is a dictionary containing 
            the gene segment as keys and the chosen alleles as values. Format is
            {Segment : [Allele, Allele ...], ...}        
        data_dict (dict): Output of load_data(). Includes family_use_dict, 
            gene_use_dict, trim_dicts, NP_transitions, NP_first_bases, NP_lengths, 
            mut_rate_per_seq and kmer_dicts.
        mutate (bool, optional): True if SHM to be introduced. Defaults to False.
        flat_usage (optional): gene, family or False. Gene or family specify that 
            sequences should use all genes or gene families evenly. If false, 
            usage follows experimental distributions. Defaults to False.
        no_trim_list (tuple, optional): List of 5 Booleans, specifying whether to not
            trim [all_ends, v_3_end, d_5_end, d_3_end, j_5_end].  Defaults to 
            (False, False, False, False, False).
        no_np_list (tuple, optional): List of 3 Booleans, specifying whether to not
            add [both_np, np1, np2]. Defaults to (False, False, False).
        shm_flat (bool): True if SHM is to be even across all sequences. Defaults to False.
        shm_random (bool): True if per base mutation is to be random. Defaults to False.
        mutation_rate (float, optional): Mutation rate to be used rather than choosing
            from distribution. Defaults to None.
        mutation_number (int, optional): Number of mutations to be added rather than
            choosing from distribution. Defaults to None.

    Returns:
        sequence (Sequence): Final recombined sequence, with trimming, 
            NP region addition and SHM if requested. 
    """

    gene_use_dict = data_dict["gene_use_dict"]
    family_use_dict = data_dict["family_use_dict"]
    trim_dicts = data_dict["trim_dicts"]
    NP_lengths = data_dict["NP_lengths"]
    NP_transitions = data_dict["NP_transitions"]
    NP_first_bases = data_dict["NP_first_bases"]
    kmer_dicts = data_dict["kmer_dicts"]
    mut_rate_per_seq = data_dict["mut_rate_per_seq"]

    functional = False
    if mutate == True:
        while functional == False:
            sequence = random_sequence(
                locus, gene_use_dict, family_use_dict, ["V", "D", "J"], flat_usage)
            attempts = 0
            while functional == False:
                attempts += 1
                if attempts > 40:
                    break
                sequence.gapped_seq = sequence.get_nuc_seq(no_trim_list,
                                                           trim_dicts,
                                                           no_np_list,
                                                           NP_lengths,
                                                           NP_transitions,
                                                           NP_first_bases,
                                                           gapped=True)
                sequence.ungapped_seq = sequence.gapped_seq.replace(".", "")
                sequence.junction_length = sequence.get_junction_length()
                this_loop = 0
                while (sequence.junction_length % 3) != 0 or check_stops(sequence.ungapped_seq) == True:
                    this_loop += 1
                    if this_loop > 80:
                        print(
                            f"couldn't make productive (before mutation): {sequence.v_allele.name}, {sequence.d_allele.name}, {sequence.j_allele.name}")
                        break
                    sequence.gapped_seq = sequence.get_nuc_seq(no_trim_list,
                                                               trim_dicts,
                                                               no_np_list,
                                                               NP_lengths,
                                                               NP_transitions,
                                                               NP_first_bases,
                                                               gapped=True)
                    sequence.ungapped_seq = sequence.gapped_seq.replace(
                        ".", "")
                    sequence.junction_length = sequence.get_junction_length()

                sequence.mutations = ""

                sequence.mutated_seq, sequence.mutations, mutation_rate = overarch_mutation(
                    sequence.gapped_seq, sequence.v_allele.family, sequence.junction_length, kmer_dicts, shm_flat, shm_random, mut_rate_per_seq, mut_rate=mutation_rate, mut_num=mutation_number)

                this_loop = 0
                while check_stops(sequence.mutated_seq) == True:
                    this_loop += 1
                    if this_loop > 50:
                        # print(
                        #     f"couldn't make productive: {sequence.v_allele.name}, {sequence.d_allele.name}, {sequence.j_allele.name}")
                        break

                    sequence.mutated_seq, sequence.mutations, mutation_rate = overarch_mutation(
                        sequence.gapped_seq, sequence.v_allele.family, sequence.junction_length, kmer_dicts, shm_flat, shm_random, mut_rate_per_seq, mut_rate=mutation_rate, mut_num=mutation_number, repeat=True)

                sequence.junction = sequence.mutated_seq[sequence.v_allele.anchor:
                                                         sequence.v_allele.anchor + sequence.junction_length].upper()

                if (sequence.junction_length % 3) == 0 and check_stops(sequence.ungapped_seq) == False:
                    sequence.junction_aa = translate(sequence.junction)
                    if sequence.junction_aa.startswith("C") and (sequence.junction_aa.endswith("F") or sequence.junction_aa.endswith("W")):
                        sequence.ungapped_seq = sequence.ungapped_seq.upper()
                        sequence.mutated_seq = sequence.mutated_seq.upper()
                        sequence.gapped_seq = sequence.gapped_seq.upper()

                        sequence.v_seq_start = 1
                        sequence.v_seq_end = sequence.v_allele.ungapped_len - sequence.v_allele.trim_3
                        sequence.d_seq_start = sequence.v_seq_end + sequence.NP1_length + 1
                        sequence.d_seq_end = sequence.d_seq_start + sequence.d_allele.ungapped_len - \
                            sequence.d_allele.trim_3 - sequence.d_allele.trim_5 - 1
                        sequence.j_seq_start = sequence.d_seq_end + sequence.NP2_length + 1
                        sequence.j_seq_end = sequence.j_seq_start + sequence.j_allele.ungapped_len - \
                            sequence.j_allele.trim_5 - 1
                        sequence.v_seq = sequence.mutated_seq[:sequence.v_seq_end]
                        sequence.d_seq = sequence.mutated_seq[sequence.d_seq_start -
                                                              1:sequence.d_seq_end]
                        sequence.j_seq = sequence.mutated_seq[sequence.j_seq_start - 1:]

                        sequence.gapped_mutated_seq = gap_seq(
                            sequence.gapped_seq, sequence.mutated_seq).upper()

                        sequence.mut_count = len(sequence.mutations)
                        sequence.mut_freq = sequence.mut_count / \
                            len(sequence.mutated_seq)

                        functional = True
                        return sequence

    if mutate == False:
        while functional == False:
            sequence = random_sequence(
                locus, gene_use_dict, family_use_dict, ["V", "D", "J"], flat_usage)
            attempts = 0
            while functional == False:
                attempts += 1
                if attempts > 40:
                    break
                sequence.gapped_seq = sequence.get_nuc_seq(no_trim_list,
                                                           trim_dicts,
                                                           no_np_list,
                                                           NP_lengths,
                                                           NP_transitions,
                                                           NP_first_bases,
                                                           gapped=True)
                sequence.ungapped_seq = sequence.gapped_seq.replace(
                    ".", "").upper()
                sequence.junction_length = sequence.get_junction_length()

                this_loop = 0
                while (sequence.junction_length % 3) != 0 or check_stops(sequence.ungapped_seq) == True:
                    this_loop += 1
                    if this_loop > 50:
                        print(
                            f"couldn't make productive: {sequence.v_allele.name}, {sequence.d_allele.name}, {sequence.j_allele.name}")
                        break
                    sequence.gapped_seq = sequence.get_nuc_seq(no_trim_list,
                                                               trim_dicts,
                                                               no_np_list,
                                                               NP_lengths,
                                                               NP_transitions,
                                                               NP_first_bases,
                                                               gapped=True)
                    sequence.ungapped_seq = sequence.gapped_seq.replace(
                        ".", "").upper()
                    sequence.junction_length = sequence.get_junction_length()

                sequence.junction = sequence.ungapped_seq[sequence.v_allele.anchor:
                                                          sequence.v_allele.anchor + sequence.junction_length].upper()

                if (sequence.junction_length % 3) == 0 and check_stops(sequence.ungapped_seq) == False:
                    sequence.junction_aa = translate(sequence.junction)
                    if sequence.junction_aa.startswith("C") and (sequence.junction_aa.endswith("F") or sequence.junction_aa.endswith("W")):
                        sequence.gapped_seq = sequence.gapped_seq.upper()
                        sequence.v_seq_start = 1
                        sequence.v_seq_end = sequence.v_allele.ungapped_len - sequence.v_allele.trim_3
                        sequence.d_seq_start = sequence.v_seq_end + sequence.NP1_length + 1
                        sequence.d_seq_end = sequence.d_seq_start + sequence.d_allele.ungapped_len - \
                            sequence.d_allele.trim_3 - sequence.d_allele.trim_5 - 1
                        sequence.j_seq_start = sequence.d_seq_end + sequence.NP2_length + 1
                        sequence.j_seq_end = sequence.j_seq_start + sequence.j_allele.ungapped_len - \
                            sequence.j_allele.trim_5 - 1
                        sequence.v_seq = sequence.ungapped_seq[:sequence.v_seq_end].upper(
                        )
                        sequence.d_seq = sequence.ungapped_seq[sequence.d_seq_start -
                                                               1:sequence.d_seq_end].upper()
                        sequence.j_seq = sequence.ungapped_seq[sequence.j_seq_start - 1:].upper(
                        )
                        functional = True
                        return sequence


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    parser = argparse.ArgumentParser(allow_abbrev=False)

    parser.add_argument("-o",
                        "--outfname",  # specify prefix to name files
                        action="store",
                        help="Output file name.",
                        dest="out_name",
                        required=True)
    parser.add_argument("--outdir",
                        action="store",
                        help="Specify the output directory.",
                        default=".")
    parser.add_argument("--datadir",
                        action="store",
                        help="Specify alternate location for input repertoire metrics.",
                        default=None)
    parser.add_argument("-n",
                        "--number_seqs",
                        type=int,
                        help="Number of sequences to simulate.",
                        default=1000
                        )
    parser.add_argument("--het",
                        # if using haplotype, decide how many of the genes should have two differing alleles
                        action="store",
                        metavar="PROP",
                        help="Proportion of genes to be heterozygous, specify as  V D J, values must be between 0 and 1. Not compatible with --all_alleles.",
                        nargs=3,
                        type=float,
                        default=[0.74, 0.16, 0.67])
    parser.add_argument("--shm",
                        action="store_true",  # SHM is false unless "-shm" is specified
                        help="Apply somatic hypermutation.")
    parser.add_argument("--shm_flat",
                        action="store_true",
                        help="Use the same mutation rate for each sequence (specify using --mut_rate or --mut_num). Default is mutation rate of 0.05 for all sequences.")
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--mut_rate",
                       action="store",
                       type=float_range(0, 0.5),
                       help="Mutation frequency for flat SHM. Value between 0 and 0.5. Cannot use with --mut_num. No default.")
    group.add_argument("--mut_num",
                       action="store",
                       type=int_range(0, 180),
                       help="Number of mutations for flat SHM. Integer between 0 and 180. Cannot use with --mut_rate. No default.")
    parser.add_argument("--shm_random",
                        action="store_true",
                        help="Do not mutate according to 5mer context. Choose position to mutate and resulting base randomly.")
    parser.add_argument("--all_alleles",  # Use the haplotype method unless "-all_alleles" is specified
                        action="store_true",
                        help="Use all available alleles from all available genes, i.e. do not generate a 'haplotype'.")
    parser.add_argument("--locus",
                        metavar="LOCUS FILE",
                        action="store",
                        help="Specify a file to use as a predefined locus. Expect a csv with two columns of alleles.")
    parser.add_argument("--flat_vdj",
                        choices=['gene', 'family', False],
                        action="store",
                        default=False,  # Family usage is biased unless "-flat" is specified
                        help="Make VDJ usage across either IMGT gene or family flat.")
    parser.add_argument("--no_trim",
                        action="store_true",  # Trimming is performed as default
                        help="Don't trim any end of V,D,J genes during recombination.")
    parser.add_argument("--no_trim_v3",
                        action="store_true",  # Trimming is performed as default
                        help="Don't trim 3' end of V genes during recombination.")
    parser.add_argument("--no_trim_d3",
                        action="store_true",  # Trimming is performed as default
                        help="Don't trim 3' end of D genes during recombination.")
    parser.add_argument("--no_trim_d5",
                        action="store_true",  # Trimming is performed as default
                        help="Don't trim 5' end of D genes during recombination.")
    parser.add_argument("--no_trim_j5",
                        action="store_true",  # Trimming is performed as default
                        help="Don't trim 5' end of J genes during recombination.")
    parser.add_argument("--no_np",
                        action="store_true",  # NP regions are generate unless otherwise specified
                        help="Don't insert nucleotides at either junction, i.e. do not create NP regions. ")
    parser.add_argument("--no_np1",
                        action="store_true",  # NP regions are generate unless otherwise specified
                        help="Don't insert nucleotides at the VD junctions, i.e. do not create NP1 regions. ")
    parser.add_argument("--no_np2",
                        action="store_true",  # NP regions are generate unless otherwise specified
                        help="Don't insert nucleotides at DJ junctions, i.e. do not create NP2 regions. ")
    parser.add_argument("--seed",
                        type=int,
                        help="Set random seed.")

    args = parser.parse_args(args)

    # Print copyright notice

    print("AIRRSHIP Copyright (C) 2022  Catherine Sutherland \nThis program comes with ABSOLUTELY NO WARRANTY. \nThis is free software, and you are welcome to redistribute \nit under certain conditions. Refer to the GNU Affero General \nPublic License for further details.")

    # Set path for reference data

    path_to_data = args.datadir

    # Process arguments

    if args.all_alleles == False:
        haplotype = True
    else:
        haplotype = False

    if args.shm == True or args.shm_flat == True or args.shm_random == True:
        mutate = True
    else:
        mutate = False

    no_trim_args = (args.no_trim, args.no_trim_v3, args.no_trim_d5,
                    args.no_trim_d3, args.no_trim_j5)

    no_np_args = (args.no_np, args.no_np1, args.no_np2)

    out_fasta = args.outdir + "/" + args.out_name + ".fasta"
    out_tsv = args.outdir + "/" + args.out_name + ".tsv"
    out_locus = args.outdir + "/" + args.out_name + "_locus.csv"
    out_summary = args.outdir + "/" + args.out_name + "_summary.txt"

    if args.seed:
        random.seed(args.seed)
    else:
        random.seed(random.random())

    # Set up mutation rates

    mutation_number = None
    mutation_rate = None
    if args.shm_flat == True:
        if args.mut_num != None:
            mutation_number = args.mut_num
            if args.mut_num > 70:
                print(
                    "Warning: High mutation number specified. Performance will be slow and may fail to generate sequences.")
        else:
            if args.mut_rate != None:
                mutation_rate = args.mut_rate
                if args.mut_rate > 0.2:
                    print(
                        "Warning: High mutation rate specified. Performance will be slow and may fail to generate sequences.")

            else:
                mutation_rate = 0.05

    # Read in all the reference data

    data_dict = load_data(path_to_data, mutate)

    # Set up the genotype to be used

    locus = get_genotype(path_to_data, args.het, haplotype, args.locus)

    # Write out locus file

    write_locus_file(out_locus, locus)

    # Write out summary file

    write_summary_file(out_summary, args)

    # TIME
    start_time = time.time()

    # OPEN FILES AND WRITE CSV HEADINGS

    print("Simulating sequences:")
    with open(out_fasta, "w") as fasta:
        with open(out_tsv, "w") as out:
            writer = csv.writer(out, delimiter="\t", lineterminator="\n")
            if mutate == True:
                writer.writerow(
                    ["sequence_id"]
                    + ["sequence"]
                    + ["v_call"]
                    + ["d_call"]
                    + ["j_call"]
                    + ["junction"]
                    + ["junction_aa"]
                    + ["junction_length"]
                    + ["np1_length"]
                    + ["np1"]
                    + ["np2_length"]
                    + ["np2"]
                    + ["v_3_trim"]
                    + ["d_5_trim"]
                    + ["d_3_trim"]
                    + ["j_5_trim"]
                    + ["v_sequence"]
                    + ["d_sequence"]
                    + ["j_sequence"]
                    + ["v_sequence_start"]
                    + ["v_sequence_end"]
                    + ["d_sequence_start"]
                    + ["d_sequence_end"]
                    + ["j_sequence_start"]
                    + ["j_sequence_end"]
                    + ["shm_events"]
                    + ["shm_count"]
                    + ["shm_freq"]
                    + ["unmutated_sequence"]
                    + ["gapped_unmutated_sequence"]
                    + ["gapped_mutated_sequence"]
                )

                # GENERATE REPERTOIRE

                counter = 0

                while counter < args.number_seqs:
                    sequence = generate_sequence(locus, data_dict, mutate, args.flat_vdj, no_trim_args,
                                                 no_np_args, args.shm_flat, args.shm_random, mutation_rate, mutation_number)
                    mutations_write = []
                    for x in sorted(sequence.mutations.items()):
                        mutations_write.append(f"{x[0]}:{x[1]}")
                    mutations_write = (','.join(mutations_write))
                    fasta.write(f">{counter}" + "\n")
                    fasta.write(sequence.mutated_seq + "\n")
                    row = (
                        [counter]
                        + [sequence.mutated_seq]
                        + [sequence.v_allele.name]
                        + [sequence.d_allele.name]
                        + [sequence.j_allele.name]
                        + [sequence.junction]
                        + [sequence.junction_aa]
                        + [sequence.junction_length]
                        + [sequence.NP1_length]
                        + [sequence.NP1_region]
                        + [sequence.NP2_length]
                        + [sequence.NP2_region]
                        + [sequence.v_allele.trim_3]
                        + [sequence.d_allele.trim_5]
                        + [sequence.d_allele.trim_3]
                        + [sequence.j_allele.trim_5]
                        + [sequence.v_seq]
                        + [sequence.d_seq]
                        + [sequence.j_seq]
                        + [sequence.v_seq_start]
                        + [sequence.v_seq_end]
                        + [sequence.d_seq_start]
                        + [sequence.d_seq_end]
                        + [sequence.j_seq_start]
                        + [sequence.j_seq_end]
                        + [mutations_write]
                        + [sequence.mut_count]
                        + [sequence.mut_freq]
                        + [sequence.ungapped_seq]
                        + [sequence.gapped_seq]
                        + [sequence.gapped_mutated_seq]
                    )
                    writer.writerow(row)
                    counter += 1
                    if counter % 1000 == 0:
                        print(counter)

            # WRITE CSV HEADINGS
            if mutate == False:
                writer.writerow(
                    ["sequence_id"]
                    + ["sequence"]
                    + ["v_call"]
                    + ["d_call"]
                    + ["j_call"]
                    + ["junction"]
                    + ["junction_aa"]
                    + ["junction_length"]
                    + ["np1_length"]
                    + ["np1"]
                    + ["np2_length"]
                    + ["np2"]
                    + ["v_3_trim"]
                    + ["d_5_trim"]
                    + ["d_3_trim"]
                    + ["j_5_trim"]
                    + ["v_sequence"]
                    + ["d_sequence"]
                    + ["j_sequence"]
                    + ["v_sequence_start"]
                    + ["v_sequence_end"]
                    + ["d_sequence_start"]
                    + ["d_sequence_end"]
                    + ["j_sequence_start"]
                    + ["j_sequence_end"]
                    + ["gapped_sequence"]
                )

                # GENERATE REPERTOIRE
                counter = 0
                non_productive = 0
                while counter < args.number_seqs:
                    sequence = generate_sequence(locus, data_dict, mutate, args.flat_vdj, no_trim_args,
                                                 no_np_args, args.shm_flat, args.shm_random, mutation_rate, mutation_number)
                    fasta.write(f">{counter}" + "\n")
                    fasta.write(sequence.ungapped_seq + "\n")
                    row = (
                        [counter]
                        + [sequence.ungapped_seq]
                        + [sequence.v_allele.name]
                        + [sequence.d_allele.name]
                        + [sequence.j_allele.name]
                        + [sequence.junction]
                        + [sequence.junction_aa]
                        + [sequence.junction_length]
                        + [sequence.NP1_length]
                        + [sequence.NP1_region]
                        + [sequence.NP2_length]
                        + [sequence.NP2_region]
                        + [sequence.v_allele.trim_3]
                        + [sequence.d_allele.trim_5]
                        + [sequence.d_allele.trim_3]
                        + [sequence.j_allele.trim_5]
                        + [sequence.v_seq]
                        + [sequence.d_seq]
                        + [sequence.j_seq]
                        + [sequence.v_seq_start]
                        + [sequence.v_seq_end]
                        + [sequence.d_seq_start]
                        + [sequence.d_seq_end]
                        + [sequence.j_seq_start]
                        + [sequence.j_seq_end]
                        + [sequence.gapped_seq]
                    )
                    writer.writerow(row)
                    counter += 1
                    if counter % 1000 == 0:
                        print(counter)

    end_time = time.time()
    print("Elapsed time was %g seconds" % (end_time - start_time))
    print(f'Output FASTA: {out_fasta}')


if __name__ == "__main__":
    main()
