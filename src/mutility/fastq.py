import pandas as pd
import gzip
import collections
from pathlib import Path
from gzip import GzipFile
from typing import BinaryIO, Union, Optional
from dataclasses import dataclass, replace

try:
    import string

    maketrans = string.maketrans
except (ImportError, NameError, AttributeError):
    maketrans = bytes.maketrans


rev_comp_table = maketrans(
    b"ACBDGHKMNSRUTWVYacbdghkmnsrutwvy", b"TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR"
)


AdapterMatch = collections.namedtuple(
    "AdapterMatch", ["astart", "astop", "rstart", "rstop", "matches", "errors"]
)


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


@dataclass
class Read:
    """Data class for sequencing reads"""

    Name: str
    Sequence: str
    Quality: str

    def reverse():
        return Read(
            self.Name, reverse_complement(self.Sequence[::-1]), self.Quality[::-1]
        )

    def __str__(self):
        return f"{self.Name}\n{self.Sequence}\n+\n{self.Quality}\n"


class Fragment:
    """Data class for single-end and paired-end Reads/Fragments."""

    def __init__(self, *reads: Read):
        self.Read1 = reads[0]
        if len(reads) == 2:
            self.Read2 = reads[1]

    @property
    def is_paired(self):
        return hasattr(self, "Read2")

    @property
    def reads(self):
        if self.is_paired:
            return [self.Read1, self.Read2]
        else:
            return [
                self.Read1,
            ]

    def __iter__(self):
        for read in self.reads:
            yield read

    def copy(self):
        return Fragment(replace(self.Read1), replace(self.Read2))

    def __str__(self):
        return f"{self.Read1}\n{self.Read2}\n"


def _open_auto(filename: str):
    if filename.endswith(".gz"):
        return gzip.open(filename, "rb")
    if filename.endswith(".bz2"):
        return bz2.open(filename, "rb")
    return open(filename, "rb", buffering=4 * 1024 * 1024)  # großer Buffer


def iterate_fastq(filename: str, reverse_reads: bool) -> Read:
    # op = mbf.align._common.BlockedFileAdaptor(filename)
    op = _open_auto(filename)
    while True:
        try:
            name = op.readline()[1:-1].decode()
            seq = op.readline()[:-1].decode()
            op.readline()
            qual = op.readline()[:-1].decode()
            if reverse_reads:
                seq = seq[::-1].translate(rev_comp_table)
                qual = qual[::-1]
            yield Read(name, seq, qual)
        except StopIteration:
            break


def get_fastq_iterator(filepath: Path):
    if filepath.suffix == ".gz":
        fileobj = gzip.open(filepath, "r")
    else:
        fileobj = filepath.open("rb")
    return read_fastq_iterator(fileobj)


def read_fastq_iterator(
    file_object: Union[BinaryIO, GzipFile], reverse_reads: bool = False
):
    """A very dump and simple fastq reader, mostly for testing the other more sophisticated variants

    Yield (seq, name, quality)
    """
    row1 = file_object.readline().decode()
    row2 = file_object.readline().decode()
    row3 = file_object.readline().decode()
    row4 = file_object.readline().decode()
    while row1:
        seq = row2[:-1]
        quality = row4[:-1]
        name = row1[1:-1]
        if reverse_reads:
            seq = seq[::-1].translate(rev_comp_table)
            qual = qual[::-1]
        yield (seq, name, quality)
        row1 = file_object.readline().decode()
        row2 = file_object.readline().decode()
        _ = file_object.readline().decode()
        row4 = file_object.readline().decode()


def count_most_common_sequences(
    output_file: Union[str, Path],
    r1: Union[str, Path],
    r2: Optional[Union[str, Path]] = None,
    max: int = 100000,
    index: Optional[int] = None,
):
    if isinstance(output_file, str):
        outfile = Path(output_file)
    else:
        outfile = output_file
    outfile.parent.mkdir(parents=True, exist_ok=True)

    iter1 = get_fastq_iterator(r1)
    iterlist = [iter1]
    if r2 is not None:
        iter2 = get_fastq_iterator(r2)
        iterlist.append(iter2)
    counter = collections.Counter()
    examples = {}
    count = 0
    for tup in zip(*iterlist):
        seqs = []
        names = []
        for read in tup:
            s, n, _ = read
            s = s[:index]
            seqs.append(s)
            names.append(n)
        key = tuple(seqs)
        counter[key] += 1
        if key not in examples:
            examples[key] = tuple(names)
        count += 1
        if count >= max:
            break
    if r2 is not None:
        to_df = {
            "Seq1": [],
            "Seq2": [],
            "Count": [],
            "Example": [],
        }
        for key in counter:
            to_df["Seq1"].append(key[0])
            to_df["Seq2"].append(key[1])
            to_df["Count"].append(counter[key])
            to_df["Example"].append(examples[key])
    else:
        to_df = {
            "Seq": [],
            "Count": [],
            "Example": [],
        }
        for key in counter:
            to_df["Seq"].append(key[0])
            to_df["Count"].append(counter[key])
            to_df["Example"].append(examples[key][0])

    df = pd.DataFrame(to_df)
    df = df.sort_values("Count", ascending=False)
    df.to_csv(output_file, sep="\t", index=False)
