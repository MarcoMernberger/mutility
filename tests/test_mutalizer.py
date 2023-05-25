import pandas as pd
from pathlib import Path
from mutility.mutalizer import (
    match_pattern_genomic_range_del,
    match_pattern_genomic_single_del,
    match_pattern_genomic_ins,
    match_pattern_genomic_sub1,
    match_pattern_genomic_ins_dup,
    match_pattern_genomic_sub_mult,
    match_pattern_protein_nonsense,
    match_pattern_protein_nonsense_2,
    match_pattern_protein_deletion,
    match_pattern_protein_frameshift,
    match_pattern_protein_delins,
    match_pattern_protein_missense,
    match_pattern_protein_synonym,
    match_pattern_effect_deletion_genomic,
    match_pattern_effect_insertion_genomic,
    match_pattern_effect_substitution,
    CodonComparison,
)

tests = [
    "NC_000017.11:g.7675248_7675250del",  # multidel
    "NC_000017.11:g.7675247del",  # single del
    "NC_000017.11:g.7675248_7675249insA",  # ins
    "NC_000017.11:g.7675248G>A",  # sub1
    "NC_000017.11:g.[7675204T>C;7675206G>A]",  # sub2
    "NC_000017.11:g.[7675204T>C;7675205G>A;7675206G>A]",  # sub3
    "NC_000017.11:g.7675249dup",  # dup or sub1
]

tests_protein = {
    "Test1": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "g.7675239delC",
    },  # syn, del
    "Test2": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "g.7675239delCA",
    },  # syn, del
    "Test3": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "g.7675239delCAG",
    },  # syn, del
    "Test4": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126Thrfs*44)",
        "Effect New": "g.7675237delGT",
    },  # fs, del
    "Test5": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126Serfs*44)",
        "Effect New": "g.7675237delGTA",
    },  # fs, del
    "Test6": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126del)",
        "Effect New": "g.7675236delTAC",
    },  # del
    "Test7": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Ala129_Leu130delinsVal)",
        "Effect New": "g.7675226delCCC",
    },  # delins
    "Test8": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr205*)",
        "Effect New": "g.7674916delTTT",
    },
    "Test9": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Lys132*)",
        "Effect New": "g.7675218ins1",
    },  # ins nonsense
    "Test10": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Lys132Argfs*17)",
        "Effect New": "g.7675217ins1",
    },  # ins frameshift
    "Test11": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "g.7675050ins1",
    },  # ins frameshift
    "Test12": {"hg38 protein": "NC_000017.11(NP_000537.3):p.(=)", "Effect New": ""},  # snp syn
    "Test13": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "p.G187G",
        "RefCodon": "GGT",
        "AltCodon": "GGA",
    },  # snp syn
    "Test14": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Gly187Glufs*22)",
        "Effect New": "",
    },  # snp fs
    "Test15": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Ala189Val)",
        "Effect New": "p.A189V",
        "RefCodon": "GCC",
        "AltCodon": "GTC",
    },  # snp missense
    "Test16": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Cys182*)",
        "Effect New": "p.C182*",
    },  # snp nonsense
    "Test17": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126Leu)",
        "Effect New": "p.Y126L",
        "RefCodon": "TAC",
        "AltCodon": "CTT",
    },  # sub missense
    "Test18": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126Arg)",
        "Effect New": "p.Y126R",
        "RefCodon": "TAC",
        "AltCodon": "CGT",
    },  # sub3 missense
    "Test19": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Tyr126Ser)",
        "Effect New": "p.Y126S",
        "RefCodon": "TAC",
        "AltCodon": "AGT",
    },  # sub2 missense
    "Test20": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Ser127Pro)",
        "Effect New": "p.S127P",
        "RefCodon": "TCC",
        "AltCodon": "CCT",
    },  # sub2 missense
    "Test21": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Ser127Tyr)",
        "Effect New": "p.S127Y",
        "RefCodon": "TCC",
        "AltCodon": "TAT",
    },  # sub1 missense
    "Test22": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "p.A129A",
        "RefCodon": "GCC",
        "AltCodon": "GCA",
    },  # sub syn
    "Test23": {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.?",
        "Effect New": "p.N131*",
    },  # sub nonsense
}


results_protein = {
    "Test1": ("p.(=)", "syn", "", "", "", "", ""),
    "Test2": ("p.(=)", "syn", "", "", "", "", ""),
    "Test3": ("p.(=)", "syn", "", "", "", "", ""),
    "Test4": ("p.Y126Tfs", "fs", "126", "", "", "Y", "T"),
    "Test5": ("p.Y126Sfs", "fs", "126", "", "", "Y", "S"),
    "Test6": ("p.Y126del", "del", "126", "", "", "Y", ""),
    "Test7": ("p.A129_L130delinsV", "delins", "129", "", "", "A", "V"),
    "Test8": ("p.Y205X", "non", "205", "", "", "Y", "X"),
    "Test9": ("p.K132X", "non", "132", "", "", "K", "X"),
    "Test10": ("p.K132Rfs", "fs", "132", "", "", "K", "R"),
    "Test11": ("p.(=)", "syn", "", "", "", "", ""),
    "Test12": ("p.(=)", "syn", "", "", "", "", ""),
    "Test13": ("p.(=)", "syn", "187", "GGT", "GGA", "G", "G"),
    "Test14": ("p.G187Efs", "fs", "187", "", "", "G", "E"),
    "Test15": ("p.A189V", "mis", "189", "GCC", "GTC", "A", "V"),
    "Test16": ("p.C182X", "non", "182", "", "", "C", "X"),
    "Test17": ("p.Y126L", "mis", "126", "TAC", "CTT", "Y", "L"),
    "Test18": ("p.Y126R", "mis", "126", "TAC", "CGT", "Y", "R"),
    "Test19": ("p.Y126S", "mis", "126", "TAC", "AGT", "Y", "S"),
    "Test20": ("p.S127P", "mis", "127", "TCC", "CCT", "S", "P"),
    "Test21": ("p.S127Y", "mis", "127", "TCC", "TAT", "S", "Y"),
    "Test22": ("p.(=)", "syn", "129", "GCC", "GCA", "A", "A"),
    "Test23": ("p.N131X", "non", "131", "", "", "N", "X"),
}


# test cases for CodonComparison
alt_seqs = [
    {
        "ID": "Ex7",
        "Sequence": "ATCTTGGGCCTGTGTTATCTCCTAGGTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAGGACTCCAAGTCAGGAGCCACTTGCCACCCTGCA",
    },
    {
        "ID": "Ex7",
        "Sequence": "ATCTTGGGCCTGTGTTATCTCCTAGCTTGGCTCTGACTGTACCACCATCCACTACAACTACATGTGTAACAGTTCCTGCATGGGCGGCATGAACCGGAGGCCCATCCTCACCATCATCACACTGGAGGACTCCAGGTCAGGAGCCACTTGCCACCCTGCA",
    },
    {
        "ID": "Ex5",
        "Sequence": "CTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATAGTGAGCAGCTGGGGCTGGAGAGACG",
    },
    {
        "ID": "Ex5",
        "Sequence": "CTCTGTCTCCTTCCTCTTCCTACAGTACTCCCCTGCCCTCAACAAGATGTTTTGCCAACTGGCCAAGACCTGCCCTGTGCAGCTGTGGGTTGATTCCACACCCCCGCCCGGCACCCGCGTCCGCGCCATGGCCATCTACAAGCAGTCACAGCACATGACGGAGGTTGTGAGGCGCTGCCCCCACCATGAGCGCTGCTCAGATAGCGATGCTGAGCAGCTGGGGCTGGAGAGACG",
    },
]


def test_match_pattern_protein_nonsense():
    row = {"hg38 protein": "NC_000017.11(NP_000537.3):p.(Leu130*)", "Effect New": "p.L130*"}
    match = match_pattern_protein_nonsense(row)
    assert match == ("p.L130X", "non", "130", "", "", "L", "X")
    matching_tests = ["Test8", "Test9", "Test16"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_nonsense(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_nonsense(tests_protein[test])
        assert match is None


def test_match_pattern_protein_nonsense_2():
    row = {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.?",
        "Effect New": "p.T23*",
    }
    match = match_pattern_protein_nonsense_2(row)
    assert match == ("p.T23X", "non", "23", "", "", "T", "X")
    matching_tests = ["Test23"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_nonsense_2(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_nonsense_2(tests_protein[test])
        assert match is None


def test_match_pattern_protein_frameshift():
    row = {"hg38 protein": "NC_000017.11(NP_000537.3):p.(Leu130Thrfs*44)"}
    match = match_pattern_protein_frameshift(row)
    assert match == ("p.L130Tfs", "fs", "130", "", "", "L", "T")
    matching_tests = ["Test4", "Test5", "Test10", "Test14"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_frameshift(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_frameshift(tests_protein[test])
        assert match is None


def test_match_pattern_protein_delins():
    row = {"hg38 protein": "NC_000017.11(NP_000537.3):p.(Ala129_Lys130delinsVal)"}
    match = match_pattern_protein_delins(row)
    assert match == ("p.A129_K130delinsV", "delins", "129", "", "", "A", "V")
    matching_tests = ["Test7"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_delins(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_delins(tests_protein[test])
        assert match is None


def test_match_pattern_protein_missense():
    row = {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(Ala129Val)",
        "Effect New": "p.A129V",
        "RefCodon": "GCC",
        "AltCodon": "GTC",
    }
    match = match_pattern_protein_missense(row)
    assert match == ("p.A129V", "mis", "129", "GCC", "GTC", "A", "V")
    matching_tests = ["Test15", "Test17", "Test18", "Test19", "Test20", "Test21"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_missense(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_missense(tests_protein[test])
        assert match is None


def test_match_pattern_protein_synonymous():
    row = {
        "hg38 protein": "NC_000017.11(NP_000537.3):p.(=)",
        "Effect New": "p.A129A",
        "RefCodon": "GCC",
        "AltCodon": "GCA",
    }
    match = match_pattern_protein_synonym(row)
    assert match == ("p.(=)", "syn", "129", "GCC", "GCA", "A", "A")
    matching_tests = ["Test1", "Test2", "Test3", "Test11", "Test12", "Test13", "Test22"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_synonym(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_synonym(tests_protein[test])
        assert match is None


def test_match_pattern_protein_deletion():
    row = {"hg38 protein": "NC_000017.11(NP_000537.3):p.(Ala129del)", "Effect New": "g.7675delGCA"}
    match = match_pattern_protein_deletion(row)
    assert match == ("p.A129del", "del", "129", "", "", "A", "")
    matching_tests = ["Test6"]
    non_matching_tests = tests_protein.keys() - matching_tests
    for test in matching_tests:
        match = match_pattern_protein_deletion(tests_protein[test])
        assert match == results_protein[test]
    for test in non_matching_tests:
        match = match_pattern_protein_deletion(tests_protein[test])
        assert match is None


def test_match_pattern_genomic_range_del():
    genomic = "NC_000017.11:g.41244988_41244989del"
    match_result = match_pattern_genomic_range_del(genomic)
    assert match_result == ("del", "del2")
    genomic = tests[0]
    match_result = match_pattern_genomic_range_del(genomic)
    assert match_result == ("del", "del3")
    for other in tests[1:]:
        match_result = match_pattern_genomic_range_del(other)
        assert match_result is None


def test_match_pattern_genomic_single_del():
    genomic = "NC_000017.11:g.41244987del"
    match = match_pattern_genomic_single_del(genomic)
    assert match == ("del", "del1")
    genomic = tests[1]
    match = match_pattern_genomic_single_del(genomic)
    assert match == ("del", "del1")
    for other in tests[2:] + [tests[0]]:
        match = match_pattern_genomic_single_del(other)
        assert match is None


def test_match_pattern_genomic_ins():
    genomic = "NC_000017.11:g.41244987_41244989insT"
    match = match_pattern_genomic_ins(genomic)
    assert match == ("ins", "insT")
    genomic = tests[2]
    match = match_pattern_genomic_ins(genomic)
    assert match == ("ins", "insA")
    for other in tests[3:] + tests[:2]:
        match = match_pattern_genomic_ins(other)
        assert match is None


def test_match_pattern_genomic_sub1():
    genomic = "NC_000017.11:g.41244987G>A"
    match = match_pattern_genomic_sub1(genomic)
    assert match == ("sub", "sub1")
    genomic = tests[3]
    match = match_pattern_genomic_sub1(genomic)
    assert match == ("sub", "sub1")
    for other in tests[4:] + tests[:3]:
        match = match_pattern_genomic_sub1(other)
        assert match is None


def test_match_pattern_genomic_sub_dup():
    genomic = "NC_000017.11:g.41244987dup"
    match = match_pattern_genomic_ins_dup(genomic)
    assert match == ("dup", "ins1")
    genomic = tests[6]
    match = match_pattern_genomic_ins_dup(genomic)
    assert match == ("dup", "ins1")
    for other in tests[:6]:
        match = match_pattern_genomic_ins_dup(other)
        assert match is None


def test_match_pattern_genomic_sub_mult():
    genomic = "NC_000017.11:g.[41244987T>C;41244989G>A]"
    match = match_pattern_genomic_sub_mult(genomic)
    assert match == ("sub", "sub2")
    genomic = tests[4]
    match = match_pattern_genomic_sub_mult(genomic)
    assert match == ("sub", "sub2")
    genomic = tests[5]
    match = match_pattern_genomic_sub_mult(genomic)
    assert match == ("sub", "sub3")
    for other in tests[6:] + tests[:4]:
        print(other)
        match = match_pattern_genomic_sub_mult(other)
        assert match is None


def test_match_pattern_effect_deletion_genomic():
    effects = ["g.7674967delG", "g.7674967delGG", "g.7674967delGGC"]
    for effect in effects:
        match = match_pattern_effect_deletion_genomic(effect)
        assert match == ("", "", "")


def test_match_pattern_effect_insertion_genomic():
    effect = "g.7674868ins1"
    match = match_pattern_effect_insertion_genomic(effect)
    assert match == ("", "", "")


def test_match_pattern_effect_substitution():
    effects = ["p.Y220C", "p.Y220*"]
    for effect in effects:
        match = match_pattern_effect_substitution(effect)
        assert match == ("Y", "220", effect[-1])


def test_codon_comparison():
    correct_results = [("AGT", "AAT"), ("GTT", "CTT"), ("GGT", "AGT"), ("", "")]
    comparator = CodonComparison(
        Path(__file__).parent / "data" / "Sequences_lib_5678_with_bbs1.tsv"
    )
    for row, result in zip(alt_seqs, correct_results):
        ref, alt = comparator.get_first_codon_difference(row)
        assert ref == result[0]
        assert alt == result[1]


def test_add_codon_columns_from_sequence():
    comparator = CodonComparison(
        Path(__file__).parent / "data" / "Sequences_lib_5678_with_bbs1.tsv"
    )
    test_df = pd.DataFrame(alt_seqs)
    df = comparator.add_codon_columns_from_sequence(test_df)
    assert df["RefCodon"].tolist() == ["AGT", "GTT", "GGT", ""]
    assert df["AltCodon"].tolist() == ["AAT", "CTT", "AGT", ""]
