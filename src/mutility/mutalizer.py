import pandas as pd
import re
from pathlib import Path
from .genomics import get_one_letter_amino_acid_code
from typing import Tuple, Union


# effect pattern
pattern_effect_substitution = re.compile(r"p\.(?P<from>\D{1})(?P<codon>\d+)(?P<to>\D{1})$")
pattern_effect_deletion_genomic = re.compile(r"g\.(?P<position>\d+)del(?P<ref>[ACTG]{1,4})$")
pattern_effect_insertion = re.compile(r"g\.(?P<position>\d+)ins\d$")
# protein pattern
pattern_protein_synonym = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\(=\)"
)  # NC_000017.11(NP_000537.3):p.(=)
pattern_protein_del = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\((?P<from>\D{3})(?P<codon>\d+)del\)"
)  # NC_000017.11(NP_000537.3):p.(Tyr126del)
pattern_protein_fs = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\((?P<from>\D{3})(?P<codon>\d+)(?P<to>\D{3})fs\*\d+\)$"
)  # NC_000017.11(NP_000537.3):p.(Tyr126Serfs*44)
pattern_protein_delins = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\((?P<from1>\D{3})(?P<codon1>\d+)_(?P<from2>\D{3})(?P<codon2>\d+)delins(?P<to>\D{3})\)"
)  # NC_000017.11(NP_000537.3):p.(Met133_Phe134delinsIle)
pattern_protein_nonsense = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\((?P<from>\D{3})(?P<codon>\d+)\*\)$"
)  # NC_000017.11(NP_000537.3):p.(Tyr220*)
pattern_protein_substitution = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\((?P<from>\D{3})(?P<codon>\d+)(?P<to>[A-Z]\D{2})\)$"
)  # NC_000017.11(NP_000537.3):p.(Tyr126Ile)
pattern_protein_nonsense_2 = re.compile(
    r"NC_000017\.11\(NP_000537\.3\):p\.\?"
)  # NC_000017.11(NP_000537.3):p.?
# genomic pattern
pattern_genomic_range_del = re.compile(r"NC_000017\.11:g\.(?P<start>\d+)_(?P<stop>\d+)del$")
pattern_genomic_single_del = re.compile(r"NC_000017\.11:g\.(?P<start>\d+)del$")
pattern_genomic_ins = re.compile(
    r"NC_000017\.11:g\.(?P<start>\d+)_(?P<stop>\d+)(?P<ins>ins[ATCG])$"
)
pattern_genomic_sub_mult = re.compile(
    r"NC_000017\.11:g\.\[(?P<sublist>((\d+[ATCG]>[ATCG]\;*)+))\]$"
)
pattern_genomic_ins_dup = re.compile(r"NC_000017\.11:g\.\d+dup$")
pattern_genomic_pattern_sub1 = re.compile(
    r"NC_000017\.11:g\.(?P<start>\d+)(?P<from>[ATCG])>(?P<to>[ATCG])$"
)


class CodonComparison:
    def __init__(self, path_to_df: Path = Path("incoming/Sequences_lib_5678_with_bbs1.tsv")):
        self.path = path_to_df
        self.init_df_exons()

    def init_df_exons(self):
        self.df_wt_exons = self.read_df()
        self.assert_frame()
        self.init_codons()

    def init_codons(self):
        wt_codons = {}
        for exon_id, row in self.df_wt_exons.iterrows():
            wt_codons[exon_id] = re.findall("...", row["Coding"])
        self.wt_codons = wt_codons

    def assert_frame(self):
        assert all(self.df_wt_exons["Coding"].str.len() % 3 == 0)

    def read_df(self):
        df_wt_exons = pd.read_csv(self.path, sep="\t")
        df_wt_exons = df_wt_exons.fillna("")
        df_wt_exons["Coding"] = df_wt_exons[
            ["beginning_of_truncated_start_codon", "Exon", "end_of_truncated_end_codon"]
        ].apply("".join, axis=1)
        df_wt_exons = df_wt_exons.set_index("ID")
        return df_wt_exons

    def get_first_codon_difference(self, row: pd.Series) -> Tuple[str, str]:
        exon_id = row["ID"]
        print(self.df_wt_exons)
        wt = self.df_wt_exons.loc[exon_id]
        sequence = "".join(
            [
                wt["beginning_of_truncated_start_codon"],
                row["Sequence"][25:-25],
                wt["end_of_truncated_end_codon"],
            ]
        )
        alt_codons = re.findall("...", sequence)
        found = False
        for codon_ref, codon_alt in zip(self.wt_codons[exon_id], alt_codons):
            if codon_ref != codon_alt:
                found = True
                break
        if found:
            return codon_ref, codon_alt
        else:
            return "", ""

    def add_codon_columns_from_sequence(self, df: pd.DataFrame) -> pd.DataFrame:
        df[["RefCodon", "AltCodon"]] = df.apply(
            self.get_first_codon_difference, axis=1, result_type="expand"
        )
        return df


def extract_codon_from_sequence(df: pd.DataFrame) -> pd.DataFrame:
    comparator = CodonComparison()
    df = comparator.add_codon_columns_from_sequence(df)
    df[
        [
            "type_g",
            "type_g_fine",
            "effect",
            "type_p",
            "codon",
            "codon_ref",
            "codon_alt",
            "aa_ref",
            "aa_alt",
        ]
    ] = df.apply(extract_mutation_details, axis=1, result_type="expand")
    return df


def extract_mutation_details(row: pd.Series) -> Tuple[str, str, str, str, str, str, str, str, str]:
    """
    Extracts mutation details from a row of a dataframe.
    """
    type_g, type_g_fine = extract_genomic_info(row)
    effect, type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt = extract_protein_info(row)
    return type_g, type_g_fine, effect, type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt


def extract_genomic_info(row: pd.Series) -> Tuple[str, str]:
    for matcher in [
        match_pattern_genomic_single_del,
        match_pattern_genomic_range_del,
        match_pattern_genomic_ins,
        match_pattern_genomic_sub_mult,
        match_pattern_genomic_sub1,
        match_pattern_genomic_ins_dup,
    ]:
        match = matcher(row["hg38_genomic"])
        if match:
            return match
    raise ValueError(f'Could not match {row["hg38_genomic"]}.')


def extract_protein_info(row: pd.Series) -> Tuple[str, str, str, str, str, str, str]:
    for matcher in [
        match_pattern_protein_missense,
        match_pattern_protein_frameshift,
        match_pattern_protein_deletion,
        match_pattern_protein_delins,
        match_pattern_protein_nonsense,
        match_pattern_protein_nonsense_2,
        match_pattern_protein_synonym,
    ]:
        match = matcher(row)
        if match:
            return match
    raise ValueError(f'Could not match {row["hg38 protein"]}, {row["Effect New"]}.')


def match_pattern_genomic_range_del(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_range_del, genomic)
    if m is not None:
        d = int(m.group("stop")) - int(m.group("start")) + 1
        result = "del", f"del{d}"
        return result
    return None


def match_pattern_genomic_single_del(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_single_del, genomic)
    if m is not None:
        result = "del", "del1"
        return result
    return None


def match_pattern_genomic_ins(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_ins, genomic)
    if m is not None:
        result = "ins", m.group("ins")
        return result
    return None


def match_pattern_genomic_sub_mult(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_sub_mult, genomic)
    if m is not None:
        result = "sub", f"sub{m.group('sublist').count(';')+1}"
        print(result)
        return result
    return None


def match_pattern_genomic_sub1(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_pattern_sub1, genomic)
    if m is not None:
        result = "sub", "sub1"
        return result
    return None


def match_pattern_genomic_ins_dup(genomic: str) -> Union[Tuple[str, str], None]:
    m = re.match(pattern_genomic_ins_dup, genomic)
    if m is not None:
        result = "dup", "ins1"
        return result
    return None


def match_pattern_protein_frameshift(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein = row["hg38 protein"]
    m = re.match(pattern_protein_fs, protein)
    if m is not None:
        codon = m.group("codon")
        aa_ref = get_one_letter_amino_acid_code(m.group("from"))
        aa_alt = get_one_letter_amino_acid_code(m.group("to"))
        type_p = "fs"
        codon_ref, codon_alt = "", ""
        result = f"p.{aa_ref}{codon}{aa_alt}fs", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
        return result
    return None


def match_pattern_protein_deletion(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein, effect = row["hg38 protein"], row["Effect New"]
    m = re.match(pattern_protein_del, protein)
    if m is not None:
        codon = m.group("codon")
        aa_ref = get_one_letter_amino_acid_code(m.group("from"))
        aa_alt = ""
        match_effect = match_pattern_effect_deletion_genomic(effect)
        if not match_effect:
            raise ValueError(f"Could not match {effect}.")
        type_p = "del"
        codon_ref, codon_alt = "", ""
        result = f"p.{aa_ref}{codon}del", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
        return result
    return None


def match_pattern_protein_delins(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein = row["hg38 protein"]
    m = re.match(pattern_protein_delins, protein)
    if m is not None:
        codon1 = m.group("codon1")
        aa_ref1 = get_one_letter_amino_acid_code(m.group("from1"))
        aa_alt = get_one_letter_amino_acid_code(m.group("to"))
        codon2 = m.group("codon2")
        aa_ref2 = get_one_letter_amino_acid_code(m.group("from2"))
        type_p = "delins"
        codon_ref, codon_alt = "", ""
        result = (
            f"p.{aa_ref1}{codon1}_{aa_ref2}{codon2}delins{aa_alt}",
            type_p,
            codon1,
            codon_ref,
            codon_alt,
            aa_ref1,
            aa_alt,
        )
        return result
    return None


def match_pattern_protein_nonsense(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein = row["hg38 protein"]
    m = re.match(pattern_protein_nonsense, protein)
    if m is not None:
        codon = m.group("codon")
        aa_ref = get_one_letter_amino_acid_code(m.group("from"))
        aa_alt = "X"
        type_p = "non"
        codon_ref, codon_alt = "", ""  # if deletion
        result = f"p.{aa_ref}{codon}{aa_alt}", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
        return result
    return None


def match_pattern_protein_nonsense_2(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein, effect = row["hg38 protein"], row["Effect New"]
    m = re.match(pattern_protein_nonsense_2, protein)
    if m is not None:
        match_effect = match_pattern_effect_substitution(effect)
        if match_effect is not None:
            aa_ref, codon, aa_alt = match_effect
        else:
            raise ValueError("Could not match {effect}.")
        aa_alt = "X"
        type_p = "non"
        codon_ref, codon_alt = "", ""
        result = f"p.{aa_ref}{codon}{aa_alt}", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
        return result
    return None


def match_pattern_protein_synonym(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein, effect = row["hg38 protein"], row["Effect New"]
    m = re.match(pattern_protein_synonym, protein)
    effect_matchers = {
        "sub": match_pattern_effect_substitution,
        "ins": match_pattern_effect_insertion_genomic,
        "del": match_pattern_effect_deletion_genomic,
    }
    if m is not None:
        type_p = "syn"
        codon = ""
        codon_ref = row["RefCodon"] if "RefCodon" in row else ""
        codon_alt = row["AltCodon"] if "AltCodon" in row else ""
        if effect != "":
            for g_type in effect_matchers:
                match_effect = effect_matchers[g_type](effect)
                if match_effect is not None:
                    aa_ref, codon, aa_alt = match_effect
                    result = "p.(=)", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
                    return result
        else:
            return "p.(=)", type_p, codon, codon_ref, codon_alt, "", ""
        raise ValueError(f"Coud not match {effect}.")
    return None


def match_pattern_protein_missense(
    row: pd.Series,
) -> Union[Tuple[str, str, str, str, str, str, str], None]:
    protein = row["hg38 protein"]
    print(protein)
    m = re.match(pattern_protein_substitution, protein)
    if m is not None:
        type_p = "mis"
        codon = m.group("codon")
        aa_ref = get_one_letter_amino_acid_code(m.group("from"))
        aa_alt = get_one_letter_amino_acid_code(m.group("to"))
        codon_ref, codon_alt = row["RefCodon"], row["AltCodon"]
        result = f"p.{aa_ref}{codon}{aa_alt}", type_p, codon, codon_ref, codon_alt, aa_ref, aa_alt
        return result
    return None


def match_pattern_effect_substitution(effect: str) -> Union[Tuple[str, str, str], None]:
    m = re.match(pattern_effect_substitution, effect)
    if m is not None:
        codon = m.group("codon")
        aa_ref = m.group("from")
        aa_alt = m.group("to")
        return (aa_ref, codon, aa_alt)
    return None


def match_pattern_effect_insertion_genomic(effect: str) -> Union[Tuple[str, str, str], None]:
    m = re.match(pattern_effect_insertion, effect)
    if m is not None:
        return ("", "", "")
    return None


def match_pattern_effect_deletion_genomic(effect: str) -> Union[Tuple[str, str, str], None]:
    m = re.match(pattern_effect_deletion_genomic, effect)
    if m is not None:
        # codon_ref = m.group("ref")
        return ("", "", "")
    return None
