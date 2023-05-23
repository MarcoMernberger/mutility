import pandas as pd
import gzip
from pathlib import Path
from gzip import GzipFile
from typing import BinaryIO, Union, Optional
import collections


def get_fastq_iterator(filepath: Path):
    if filepath.suffix == ".gz":
        fileobj = gzip.open(filepath, "r")
    else:
        fileobj = filepath.open("rb")
    return read_fastq_iterator(fileobj)


def read_fastq_iterator(file_object: Union[BinaryIO, GzipFile]):
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
