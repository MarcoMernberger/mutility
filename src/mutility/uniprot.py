#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
uniprot.py: Contains some functions for retrieving ENSEMBL ids for certain
protein classes.
"""

from pathlib import Path
from typing import List, Dict
from pypipegraph import Job
from pandas import DataFrame
import pandas as pd
import requests
import parse
import re

__author__ = "Marco Mernberger"
__copyright__ = "Copyright (c) 2020 Marco Mernberger"
__license__ = "mit"


def get_kinase_from_url(species: str = "mouse") -> DataFrame:
    """
    Retrieves a list of kinases from the uniprot server and returns a DataFrame
    with kinase uniprot ids and corresponding classes.

    Parameters
    ----------
    species : str, optional
        Species to use, must be either "mouse" or "human", by default "mouse".

    Returns
    -------
    DataFrame
        A DataFrame with kinase uniprot ids and corresponding classes.

    Raises
    ------
    ValueError
        If the wrong species was provided.
    """
    if species not in ["mouse", "human"]:
        raise ValueError(
            "SwissProt has only information for human and mouse kinases, set species accordingly."
        )
    species = species.upper()
    pattern = species + "{uniprot})"
    url_for_kinases = "https://www.uniprot.org/docs/pkinfam.txt"
    r = requests.get(url_for_kinases)
    to_df: Dict[str, List] = {"uniprot_id": [], "kinase class": []}
    for p in parse.findall("=\n{group}\n{equals}=\n{kinases}\n\n", r.text):
        group = p["group"]
        if not re.match(r"\s", group):
            kinase_block = p["kinases"]
            for k in parse.findall(pattern, kinase_block):
                uniprot = k["uniprot"].strip()[1:]
                to_df["uniprot_id"].append(uniprot)
                to_df["kinase class"].append(group)
    return pd.DataFrame(to_df)


def uniprot_to_ensenbl(uniprot_ids: List[str]) -> DataFrame:
    """
    Converts a list of uniprot IDs to Ensembl IDs using the Uniprot conversion
    tool.

    Parameters
    ----------
    uniprot_ids : List[str]
        List of IDs to convert.

    Returns
    -------
    DataFrame
        DataFrame with uniprot IDs and corresponding Ensembl IDs.
    """
    query = " ".join(uniprot_ids)
    url_uniprot = "https://www.uniprot.org/uploadlists/"
    payload = {
        "from": "ACC+ID",
        "to": "ENSEMBL_ID",
        "format": "tab",
        "query": query,
    }
    r = requests.get(url_uniprot, params=payload)
    to_df: Dict[str, List] = {"uniprot": [], "gene_stable_id": []}
    lines = r.text.split("\n")
    for line in lines[1:]:
        try:
            u, e = line.split("\t")
            to_df["uniprot"].append(u)
            to_df["gene_stable_id"].append(e)
        except ValueError:
            pass
    return pd.DataFrame(to_df)


def get_kinases_with_ensembl_ids(species: str) -> DataFrame:
    """
    Returns a DataFrame containing all kinase for the given species from uniprot,
    their corresponding enzyme class and ensembl stable ids.

    Parameters
    ----------
    species : str
        The species to use must be either "human" or "mouse".

    Returns
    -------
    DataFrame
        DataFrame containing ensembl ids, uniprot ids and kinase class.
    """
    df_uniprot = get_kinase_from_url(species)
    df_uniprot.index = df_uniprot["uniprot_id"]
    df_ensembl = uniprot_to_ensenbl(df_uniprot["uniprot_id"].values)
    df_ensembl.index = df_ensembl["uniprot"]
    df = df_uniprot.join(df_ensembl)
    df = df[["gene_stable_id", "uniprot_id", "kinase class"]]
    df = df.drop_duplicates(
        subset=["gene_stable_id"]
    )  # we want that unique for indexing
    assert len(df["gene_stable_id"].unique()) == len(df)
    return df


def write_kinases(outfile: Path, species: str) -> Job:
    """
    Returns a job that writes a DataFrame containing all kinase for the given
    species from uniprot, their corresponding enzyme class and ensembl stable ids.

    Parameters
    ----------
    outfile : Path
        The filepath to write the DataFrame to.
    species : str
        The species to use must be either "human" or "mouse".

    Returns
    -------
    Job
        The Job that writes the DataFrame.
    """
    outfile.parent.mkdir(parents=True, exist_ok=True)

    def __write():
        df = get_kinases_with_ensembl_ids(species)
        df.to_csv(outfile, sep="\t", index=False)

    return ppg.FileGeneratingJob(outfile, __write).depends_on(
        ppg.FunctionInvariant(
            "get_kinases_with_ensembl_ids", get_kinases_with_ensembl_ids
        )
    )
