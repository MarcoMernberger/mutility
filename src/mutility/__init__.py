# -*- coding: utf-8 -*-
from pkg_resources import get_distribution, DistributionNotFound

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = __name__
    __version__ = get_distribution(dist_name).version
except DistributionNotFound:
    __version__ = "unknown"
finally:
    del get_distribution, DistributionNotFound

from .functions import *  # noqa: E403
from .frames import read_excel_from_biologists  # noqa: E401
from .genomics import (
    reverse_complement,  # noqa: E401
    get_one_letter_amino_acid_code,  # noqa: E401
    get_three_letter_amino_acid_code,  # noqa: E401
)
from .mutalizer import CodonComparison, extract_codon_from_sequence  # noqa: E401
from .fastq import count_most_common_sequences
