import pandas as pd
import re
from typing import Union, Optional
from pathlib import Path
from pandas import DataFrame


def read_excel_from_biologists(
    infile: Union[Path, str], sheet_name: Optional[Union[str, int]] = 0, **kwargs
) -> DataFrame:
    df = pd.read_excel(infile, sheet_name=sheet_name, **kwargs)
    print(df)
    # whitespace = re.compile(r"\s+")
    for col in df.columns:
        if df[col].dtype == "object":
            df[col].str.strip()
    return df
