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

from .functions import *
from .frames import read_excel_from_biologists
from .genomics import (
    reverse_complement,
    get_one_letter_amino_acid_code,
    get_three_letter_amino_acid_code,
)
