import os

import pytest

ROOT = os.path.dirname(os.path.dirname(__file__))
TESTS = os.path.join(ROOT, "static", "tests")
DATA = os.path.join(TESTS, "formats")
SYNTH = os.path.join(TESTS, "synthetic")


@pytest.fixture(scope="session", autouse=True)
def _ensure_fixtures():
    """Build the BED/GFF3/GTF + synthetic datasets if they are not present."""
    if not (os.path.exists(os.path.join(DATA, "genes.gff3"))
            and os.path.exists(os.path.join(SYNTH, "celegans_20", "genes.tsv"))):
        import runpy

        runpy.run_path(os.path.join(os.path.dirname(__file__), "make_fixtures.py"),
                       run_name="__main__")


@pytest.fixture(scope="session")
def genome_txt():
    return os.path.join(TESTS, "chrs_test_2.txt")


@pytest.fixture(scope="session")
def genome_fasta():
    return os.path.join(TESTS, "chrs_test_3.fasta")


@pytest.fixture(scope="session")
def genome():
    from gfviewer.genome import load_genome

    return load_genome(os.path.join(TESTS, "chrs_test_2.txt"))


@pytest.fixture(scope="session")
def data_csv():
    return os.path.join(TESTS, "data_test_2.csv")


@pytest.fixture(scope="session")
def data_xlsx():
    return os.path.join(TESTS, "data_test_1.xlsx")


@pytest.fixture(scope="session")
def data_tsv():
    return os.path.join(TESTS, "data_test_3.tsv")


@pytest.fixture(scope="session")
def bed_file():
    return os.path.join(DATA, "genes.bed")


@pytest.fixture(scope="session")
def gff3_file():
    return os.path.join(DATA, "genes.gff3")


@pytest.fixture(scope="session")
def gtf_file():
    return os.path.join(DATA, "genes.gtf")


@pytest.fixture(scope="session")
def per_family_beds():
    d = os.path.join(DATA, "per_family")
    return [os.path.join(d, f) for f in sorted(os.listdir(d))]


@pytest.fixture(scope="session")
def named_bed():
    return os.path.join(DATA, "genes_named.bed")


@pytest.fixture(scope="session")
def mapping_file():
    return os.path.join(DATA, "mapping.tsv")


@pytest.fixture(scope="session", params=["arabidopsis_10", "celegans_20"])
def synthetic_dataset(request):
    d = os.path.join(SYNTH, request.param)
    return os.path.join(d, "genes.tsv"), os.path.join(d, "genome.txt")
