"""Build every bundled example dataset into a single directory: static/tests/.

Run once (also invoked automatically by the test-suite conftest):

    python tests/make_fixtures.py

Outputs
-------
static/tests/formats/            the 6-family Babesia MO1 set in every input
                                 format (BED, per-family BED, BED + mapping,
                                 GFF3, GTF)
static/tests/synthetic/          two synthetically generated sets:
    arabidopsis_10/              10 random gene families on the 5 A. thaliana
                                 (TAIR10) chromosomes
    celegans_20/                 20 random gene families on the 6 C. elegans
                                 (WBcel235) chromosomes
"""

import os
import random

import pandas as pd

HERE = os.path.dirname(__file__)
TESTS = os.path.join(HERE, "..", "static", "tests")
SRC = os.path.join(TESTS, "data_test_2.csv")
FORMATS = os.path.join(TESTS, "formats")
SYNTH = os.path.join(TESTS, "synthetic")


# --------------------------------------------------------------------------- #
# 1. the 6-family Babesia set re-expressed as BED / GFF3 / GTF
# --------------------------------------------------------------------------- #
def build_formats():
    os.makedirs(FORMATS, exist_ok=True)
    df = pd.read_csv(SRC)
    df = df[df["gene_family"].str.lower() != "centromere"]

    with open(os.path.join(FORMATS, "genes.bed"), "w") as fh:
        fh.write("track name=mgf\n")
        for _, r in df.iterrows():
            fh.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                r["chromosome"], int(r["start"]) - 1, int(r["end"]),
                r["gene_family"], r["strand"]))

    fam_dir = os.path.join(FORMATS, "per_family")
    os.makedirs(fam_dir, exist_ok=True)
    for fam, sub in df.groupby("gene_family"):
        with open(os.path.join(fam_dir, "{}.bed".format(fam)), "w") as fh:
            for _, r in sub.iterrows():
                fh.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                    r["chromosome"], int(r["start"]) - 1, int(r["end"]),
                    r["gene_id"], r["strand"]))

    with open(os.path.join(FORMATS, "genes_named.bed"), "w") as fh:
        for _, r in df.iterrows():
            fh.write("{}\t{}\t{}\t{}\t0\t{}\n".format(
                r["chromosome"], int(r["start"]) - 1, int(r["end"]),
                r["gene_id"], r["strand"]))
    with open(os.path.join(FORMATS, "mapping.tsv"), "w") as fh:
        fh.write("gene_id\tgene_family\n")
        for _, r in df.iterrows():
            fh.write("{}\t{}\n".format(r["gene_id"], r["gene_family"]))

    with open(os.path.join(FORMATS, "genes.gff3"), "w") as fh:
        fh.write("##gff-version 3\n")
        for _, r in df.iterrows():
            fh.write("{}\tgfviewer\tgene\t{}\t{}\t.\t{}\t.\tID={};gene_family={}\n".format(
                r["chromosome"], int(r["start"]), int(r["end"]), r["strand"],
                r["gene_id"], r["gene_family"]))

    with open(os.path.join(FORMATS, "genes.gtf"), "w") as fh:
        for _, r in df.iterrows():
            fh.write('{}\tgfviewer\tgene\t{}\t{}\t.\t{}\t.\tgene_id "{}"; gene_family "{}";\n'.format(
                r["chromosome"], int(r["start"]), int(r["end"]), r["strand"],
                r["gene_id"], r["gene_family"]))

    with open(os.path.join(FORMATS, "README.txt"), "w") as fh:
        fh.write(
            "The 6 Babesia sp. MO1 gene families (see data_test_2.csv) written in "
            "every input format GFViewer accepts.\n\n"
            "  genes.bed         one BED file  -> column 4 (name) is the family\n"
            "  per_family/*.bed  several BED files -> each file name is a family\n"
            "  genes_named.bed   BED with gene ids in column 4 ...\n"
            "  mapping.tsv       ... + this gene_id -> gene_family map\n"
            "  genes.gff3        GFF3, family in the gene_family attribute\n"
            "  genes.gtf         GTF, family in the gene_family attribute\n\n"
            "Genome for all of the above: ../chrs_test_2.txt\n"
        )


# --------------------------------------------------------------------------- #
# 2. synthetic datasets
# --------------------------------------------------------------------------- #
# reference chromosome lengths (bp)
ARABIDOPSIS = {          # Arabidopsis thaliana, TAIR10 nuclear chromosomes
    "Chr1": 30_427_671, "Chr2": 19_698_289, "Chr3": 23_459_830,
    "Chr4": 18_585_056, "Chr5": 26_975_502,
}
CELEGANS = {             # Caenorhabditis elegans, WBcel235
    "I": 15_072_434, "II": 15_279_421, "III": 13_783_801,
    "IV": 17_493_829, "V": 20_924_180, "X": 17_718_942,
}

_MODES = ("tandem", "subtelomeric", "chromosome-enriched", "dispersed")


def _family_genes(rng, genome, fam, prefix, start_id):
    """Return a list of (gene_id, family, chrom, start, end, strand, mode)."""
    chroms = list(genome)
    mode = rng.choices(_MODES, weights=(3, 2, 2, 3))[0]
    n = rng.randint(4, 46)
    lo, hi = 800, 4500
    out = []
    gid = start_id

    def add(c, s, e, strand):
        nonlocal gid
        gid += 1
        s = max(1, min(int(s), genome[c] - 2))
        e = max(s + 1, min(int(e), genome[c] - 1))
        out.append(("{}g{:05d}".format(prefix, gid), fam, c, s, e, strand, mode))

    if mode == "tandem":
        for _ in range(rng.randint(1, 2)):
            c = rng.choice(chroms)
            L = genome[c]
            anchor = rng.choice([
                rng.randint(30_000, int(L * 0.06)),
                rng.randint(int(L * 0.40), int(L * 0.60)),
                rng.randint(int(L * 0.93), L - 120_000),
            ])
            pos = anchor
            for _ in range(max(2, n // 2)):
                glen = rng.randint(lo, hi)
                add(c, pos, pos + glen, rng.choice("+-"))
                pos += glen + rng.randint(150, 6_000)
    elif mode == "subtelomeric":
        for _ in range(n):
            c = rng.choice(chroms)
            L = genome[c]
            if rng.random() < 0.5:
                s = rng.randint(5_000, int(L * 0.12))
            else:
                s = rng.randint(int(L * 0.88), L - 8_000)
            glen = rng.randint(lo, hi)
            add(c, s, s + glen, rng.choice("+-"))
    elif mode == "chromosome-enriched":
        home = rng.choice(chroms)
        for _ in range(n):
            c = home if rng.random() < 0.82 else rng.choice(chroms)
            L = genome[c]
            s = rng.randint(20_000, L - 8_000)
            glen = rng.randint(lo, hi)
            add(c, s, s + glen, rng.choice("+-"))
    else:  # dispersed
        for _ in range(n):
            c = rng.choice(chroms)
            L = genome[c]
            s = rng.randint(20_000, L - 8_000)
            glen = rng.randint(lo, hi)
            add(c, s, s + glen, rng.choice("+-"))
    return out, gid


def build_synthetic_one(key, genome, n_families, prefix, seed, blurb):
    out_dir = os.path.join(SYNTH, key)
    os.makedirs(out_dir, exist_ok=True)
    rng = random.Random(seed)

    with open(os.path.join(out_dir, "genome.txt"), "w") as fh:
        for c, L in genome.items():
            fh.write("{},{}\n".format(c, L))

    rows = []
    modes = {}
    gid = 0
    for fi in range(1, n_families + 1):
        fam = "{}_GF{:02d}".format(prefix, fi)
        recs, gid = _family_genes(rng, genome, fam, prefix.lower(), gid)
        rows.extend(recs)
        modes[fam] = recs[0][6] if recs else "dispersed"

    df = pd.DataFrame(
        [(r[0], r[1], r[2], r[3], r[4], r[5]) for r in rows],
        columns=["gene_id", "gene_family", "chromosome", "start", "end", "strand"],
    ).sort_values(["chromosome", "start"]).reset_index(drop=True)
    df.to_csv(os.path.join(out_dir, "genes.tsv"), sep="\t", index=False)

    with open(os.path.join(out_dir, "colors.txt"), "w") as fh:
        for i, fam in enumerate(sorted(modes), 1):
            fh.write("{},{}\n".format(fam, (i - 1) % 20 + 1))

    with open(os.path.join(out_dir, "README.txt"), "w") as fh:
        fh.write(blurb.rstrip() + "\n\n")
        fh.write("Files:\n"
                 "  genes.tsv    gene_id, gene_family, chromosome, start, end, strand "
                 "(1-based, inclusive)\n"
                 "  genome.txt   chromosome,length\n"
                 "  colors.txt   gene_family,palette-index (optional)\n\n")
        fh.write("Per-family spatial pattern used by the generator:\n")
        for fam in sorted(modes):
            fh.write("  {:<12s} {}\n".format(fam, modes[fam]))
    return len(df)


def build_synthetic():
    os.makedirs(SYNTH, exist_ok=True)
    n1 = build_synthetic_one(
        "arabidopsis_10", ARABIDOPSIS, 10, "AT", seed=20240829,
        blurb="Synthetic dataset: 10 randomly generated gene families placed on the "
              "five nuclear chromosomes of Arabidopsis thaliana (TAIR10 lengths). "
              "Each family follows one of four spatial patterns (tandem array, "
              "sub-telomeric, single-chromosome-enriched, dispersed) so the "
              "analytics have something to find. Not real genes.",
    )
    n2 = build_synthetic_one(
        "celegans_20", CELEGANS, 20, "CE", seed=13372025,
        blurb="Synthetic dataset: 20 randomly generated gene families placed on the "
              "six chromosomes of Caenorhabditis elegans (WBcel235 lengths). "
              "Each family follows one of four spatial patterns (tandem array, "
              "sub-telomeric, single-chromosome-enriched, dispersed). Not real genes.",
    )
    return n1, n2


def main():
    build_formats()
    n1, n2 = build_synthetic()
    print("wrote", os.path.relpath(FORMATS))
    print("wrote", os.path.relpath(os.path.join(SYNTH, "arabidopsis_10")),
          "({} genes)".format(n1))
    print("wrote", os.path.relpath(os.path.join(SYNTH, "celegans_20")),
          "({} genes)".format(n2))


if __name__ == "__main__":
    main()
