#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""functions.py: Contains some utility functions that go nowhere else."""

from typing import Optional, Callable, List, Dict, Tuple, Any
from IPython.display import Markdown, display
import numpy as np


__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def dm(text):
    "display helper for lazy typers"
    display(Markdown(text))


def filter_function(
    column_names: List[str],
    threshold: float,
    canonical_chromosomes: List[str],
    biotypes: List[str] = None,
) -> Callable:
    """
    Filter function for RNAseq used to filter genes instances.

    This is my regular goto filter for DEG analysis. Filters for canonical chromosomes
    and biotypes. In addition, you can set a threshold on expression value columns.
    This is used to prefilter the genes before the actual DEG analysis to limit the
    set of genes we need to consider and exclude lowly expressed genes.

    Parameters
    ----------
    column_names : List[str]
        List of expression columns to filter for.
    threshold : float
        Threshold for expression values to be considered as 'measured'.
    canonical_chromosomes : List[str]
        List of canonical chromosomes to consider.
    biotypes : List[str], optional
        Biotypes that are relevant, by default None.

    Returns
    -------
    Callable
        Filter function t be passed to genes.filter.
    """

    def __filter(df):
        keep = np.zeros(len(df), dtype=bool)
        for column_name in column_names:
            keep = keep | (df[column_name] >= threshold)
        chromosomes = [x in canonical_chromosomes for x in df["chr"].values]
        keep = keep & chromosomes
        if biotypes is not None:
            bio = [x in biotypes for x in df["biotype"].values]
            keep = keep & bio
        return keep

    return __filter


def get_label_fuction(columns_names_nice, row_names_nice):
    def label_func(label):
        if label in columns_names_nice:
            return columns_names_nice[label]
        elif label in row_names_nice:
            return row_names_nice[label]
        return label

    return label_func


def dict_to_string_of_items(dictionary: dict) -> str:
    return " ".join([str(item) for tuple in dictionary.items() for item in tuple if item!=""])
