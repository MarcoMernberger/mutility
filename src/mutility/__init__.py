# -*- coding: utf-8 -*-
try:
    # Py3.8+
    from importlib.metadata import (
        version as _get_version,
        PackageNotFoundError as _PackageNotFoundError,
    )
except Exception:
    # Fallback for older python
    from importlib_metadata import version as _get_version, PackageNotFoundError as _PackageNotFoundError  # type: ignore

# Change here if project/distribution name differs from package name
dist_name = __name__
try:
    __version__ = _get_version(dist_name)
except _PackageNotFoundError:
    __version__ = "unknown"
finally:
    del _get_version, _PackageNotFoundError

from .functions import *  # noqa: E403
from .frames import read_excel_from_biologists  # noqa: E401
from .genomics import (
    reverse_complement,  # noqa: E401
    get_one_letter_amino_acid_code,  # noqa: E401
    get_three_letter_amino_acid_code,  # noqa: E401
)
from .mutalizer import CodonComparison, extract_codon_from_sequence  # noqa: E401
from .fastq import count_most_common_sequences


__all__ = [
    "reverse_complement",
    "get_one_letter_amino_acid_code",
    "get_three_letter_amino_acid_code",
    "CodonComparison",
    "extract_codon_from_sequence",
    "read_excel_from_biologists",
    "count_most_common_sequences",
]
