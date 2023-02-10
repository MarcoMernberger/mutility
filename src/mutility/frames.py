import pandas as pd
import re
from typing import Union, Optional, List
from pathlib import Path
from pandas import DataFrame


def trim_dna_strings(sequence: str) -> str:
    sequence = str(sequence)
    return sequence.strip().upper()


def read_excel_from_biologists(
    infile: Union[Path, str],
    sheet_name: Optional[Union[str, int]] = 0,
    trim_columns: Optional[List[str]] = None,
    merge_columns: List[str] = [],
    **kwargs
) -> DataFrame:
    converters = None
    if trim_columns is not None:
        converters = dict(zip(trim_columns, [trim_dna_strings] * len(trim_columns)))
    df = pd.read_excel(infile, sheet_name=sheet_name, converters=converters, **kwargs)
    # deal with merged cells
    for col in merge_columns:
        df[col] = df[col].fillna(method="ffill")
    return df
