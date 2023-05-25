#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""util.py: Contains utility functions for genomics."""

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans

rev_comp_table = maketrans(b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR")


def reverse_complement(sequence: str) -> str:
    """
    reverse_complement retuzrns the reverse complement of given sequence.

    Parameters
    ----------
    sequence : str
        Input sequence.

    Returns
    -------
    str
        Reverse complement of input sequence.
    """
    return sequence[::-1].translate(rev_comp_table)


three_to_one = {
    "Ala": "A",  # Alanine
    "Arg": "R",  # Arginine
    "Asn": "N",  # Asparagine
    "Asp": "D",  # Aspartic acid
    "Cys": "C",  # Cysteine
    "Gln": "Q",  # Glutamine
    "Glu": "E",  # Glutamic acid
    "Gly": "G",  # Glycine
    "His": "H",  # Histidine
    "Ile": "I",  # Isoleucine
    "Leu": "L",  # Leucine
    "Lys": "K",  # Lysine
    "Met": "M",  # Methionine
    "Phe": "F",  # Phenylalanine
    "Pro": "P",  # Proline
    "Pyl": "O",  # Pyrrolysine
    "Ser": "S",  # Serine
    "Sec": "U",  # Selenocysteine
    "Thr": "T",  # Threonine
    "Trp": "W",  # Tryptophan
    "Tyr": "Y",  # Tyrosine
    "Val": "V",  # Valine
    "Asx": "B",  # Aspartic acid or Asparagine
    "Glx": "Z",  # Glutamic acid or Glutamine
    "Xaa": "X",  # Any amino acid
    "Xle": "J",  # Leucine or Isoleucine
    "TERM": "*",  # Termination
}

one_to_three = {v: k for k, v in three_to_one.items()}


def get_one_letter_amino_acid_code(three_letter_code: str) -> str:
    """
    get_one_letter_amino_acid_code returns the one letter code for the given three letter code.

    Parameters
    ----------
    three_letter_code : str
        Three letter amino acid code.

    Returns
    -------
    str
        One letter amino acid code.
    """
    return three_to_one[three_letter_code]


def get_three_letter_amino_acid_code(one_letter_code: str) -> str:
    """
    get_three_letter_amino_acid_code returns the three letter code for the given one letter code.

    Parameters
    ----------
    one_letter_code : str
        One letter amino acid code.

    Returns
    -------
    str
        Three letter amino acid code.
    """
    return one_to_three[one_letter_code]
