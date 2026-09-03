import pytest

from gfviewer.errors import InputValidationError
from gfviewer.genome import (
    chromosome_id_list,
    classify_sequences,
    load_genome,
)


def test_load_txt_table(genome_txt):
    g = load_genome(genome_txt)
    assert list(g) == ["chr1", "chr2", "chr3"]
    assert g["chr1"] == 4143664
    assert all(isinstance(v, int) for v in g.values())


def test_load_fasta_lengths(genome_fasta):
    g = load_genome(genome_fasta)
    assert set(g) == {"chr1", "chr2", "chr3", "chr4", "chr5"}
    assert all(v > 0 for v in g.values())


def test_order_is_preserved(tmp_path):
    p = tmp_path / "g.txt"
    p.write_text("zeta,100\nalpha,200\nmu,50\n")
    g = load_genome(str(p))
    assert list(g) == ["zeta", "alpha", "mu"]


def test_tab_and_header_row(tmp_path):
    p = tmp_path / "g.tsv"
    p.write_text("name\tlength\nchrA\t1000\nchrB\t2000\n")
    g = load_genome(str(p))
    assert g == {"chrA": 1000, "chrB": 2000}


def test_missing_file():
    with pytest.raises(InputValidationError):
        load_genome("/no/such/file.txt")


def test_bad_length(tmp_path):
    p = tmp_path / "g.txt"
    p.write_text("chrA,notanumber\n")
    with pytest.raises(InputValidationError):
        load_genome(str(p))


def test_classify_sequences_name_based():
    g = {
        "chr1": 4_000_000, "chr2": 3_000_000, "chrX": 2_000_000,
        "scaffold_12": 40_000, "contig00007": 8_000,
        "chrUn_random": 12_000, "MT": 16_000,
    }
    flags = classify_sequences(g)
    assert flags["chr1"] and flags["chr2"] and flags["chrX"]
    assert not flags["scaffold_12"]
    assert not flags["contig00007"]
    assert not flags["chrUn_random"]
    assert not flags["MT"]
    assert chromosome_id_list(g) == ["chr1", "chr2", "chrX"]


def test_classify_roman_and_arm_names():
    g = {"I": 900_000, "II": 800_000, "2L": 700_000, "X": 600_000, "scf718": 2_000}
    flags = classify_sequences(g)
    assert all(flags[s] for s in ("I", "II", "2L", "X"))
    assert not flags["scf718"]


def test_classify_explicit_override():
    g = {"a": 100, "b": 200, "c": 300}
    flags = classify_sequences(g, chromosome_ids=["a", "c"])
    assert flags == {"a": True, "b": False, "c": True}


def test_classify_never_empties_everything():
    g = {"weirdname1": 100, "weirdname2": 200}
    assert all(classify_sequences(g).values())
