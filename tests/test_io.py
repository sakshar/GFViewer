import pandas as pd
import pytest

from gfviewer.errors import InputValidationError
from gfviewer.io import NORMALISED_COLUMNS, load_features


def _families(df):
    return set(df[df["kind"] == "gene"]["gene_family"])


def test_table_csv(genome, data_csv):
    df, warn = load_features(data_csv, genome)
    assert list(df.columns) == NORMALISED_COLUMNS
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}
    assert (df["start"] <= df["end"]).all()
    assert set(df["strand"]) <= {"+", "-", "."}


def test_table_xlsx(genome, data_xlsx):
    from gfviewer.genome import load_genome
    import os

    g = load_genome(os.path.join(os.path.dirname(data_xlsx), "chrs_test_1.txt"))
    df, _ = load_features(data_xlsx, g)
    assert len(_families(df)) == 18


def test_tsv_with_centromeres(genome_fasta):
    import os
    from gfviewer.genome import load_genome

    g = load_genome(genome_fasta)
    tsv = os.path.join(os.path.dirname(genome_fasta), "data_test_3.tsv")
    df, _ = load_features(tsv, g)
    assert (df["kind"] == "centromere").sum() == 5
    cen = df[df["kind"] == "centromere"]
    assert set(cen["strand"]) == {"."}
    assert (cen["gene_family"] == "centromere").all()


def test_bed_is_zero_based_converted(genome, bed_file, data_csv):
    bed_df, _ = load_features(bed_file, genome)
    csv_df, _ = load_features(data_csv, genome)
    csv_df = csv_df[csv_df["kind"] == "gene"].reset_index(drop=True)
    bed_df = bed_df[bed_df["kind"] == "gene"].reset_index(drop=True)
    # BED name column became the family; coordinates must line up 1:1 with the CSV
    merged = csv_df.merge(
        bed_df, on=["chromosome", "start", "end"], suffixes=("_csv", "_bed")
    )
    assert len(merged) == len(csv_df)
    assert (merged["gene_family_csv"] == merged["gene_family_bed"]).all()


def test_bed_single_file_name_is_family(genome, bed_file):
    df, _ = load_features(bed_file, genome)
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}


def test_per_family_beds(genome, per_family_beds):
    df, _ = load_features(per_family_beds, genome)
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}


def test_bed_plus_mapping(genome, named_bed, mapping_file):
    df, _ = load_features(named_bed, genome, mapping_file=mapping_file)
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}


def test_gff3(genome, gff3_file):
    df, _ = load_features(gff3_file, genome)
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}
    assert (df["start"] <= df["end"]).all()


def test_gtf(genome, gtf_file):
    df, _ = load_features(gtf_file, genome)
    assert _families(df) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}


def test_gff3_matches_csv_coordinates(genome, gff3_file, data_csv):
    gff, _ = load_features(gff3_file, genome)
    csv, _ = load_features(data_csv, genome)
    gcmp = gff[gff["kind"] == "gene"][["chromosome", "start", "end"]].sort_values(
        ["chromosome", "start"]
    ).reset_index(drop=True)
    ccmp = csv[csv["kind"] == "gene"][["chromosome", "start", "end"]].sort_values(
        ["chromosome", "start"]
    ).reset_index(drop=True)
    pd.testing.assert_frame_equal(gcmp, ccmp)


def test_unknown_chromosome_errors(genome, tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chrX,10,20,+\n"
    )
    with pytest.raises(InputValidationError):
        load_features(str(p), genome)


def test_unknown_chromosome_drop(genome, tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chr1,10,20,+\n"
        "g2,F1,chrX,10,20,+\n"
    )
    df, warn = load_features(str(p), genome, on_unknown_chrom="drop")
    assert len(df) == 1
    assert any("Dropped" in w for w in warn)


def test_missing_columns(genome, tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text("gene_id,chromosome,start,end\ng1,chr1,10,20\n")
    with pytest.raises(InputValidationError):
        load_features(str(p), genome)


def test_non_numeric_coords(genome, tmp_path):
    p = tmp_path / "bad.csv"
    p.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chr1,ten,20,+\n"
    )
    with pytest.raises(InputValidationError):
        load_features(str(p), genome)


def test_coordinate_beyond_chromosome_is_clipped(genome, tmp_path):
    p = tmp_path / "over.csv"
    p.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chr1,10,99999999999,+\n"
    )
    df, warn = load_features(str(p), genome, coord_bounds="clip")
    assert df.iloc[0]["end"] == genome["chr1"]
    assert any("past the chromosome" in w for w in warn)


def test_start_end_swap(genome, tmp_path):
    p = tmp_path / "swap.csv"
    p.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chr1,500,100,+\n"
    )
    df, warn = load_features(str(p), genome)
    row = df.iloc[0]
    assert row["start"] == 100 and row["end"] == 500
    assert any("swapped" in w for w in warn)


def test_header_synonyms(genome, tmp_path):
    p = tmp_path / "syn.tsv"
    p.write_text(
        "ID\tfamily\tchrom\tbegin\tstop\tstrand\n"
        "g1\tF1\tchr1\t10\t20\t+\n"
    )
    df, _ = load_features(str(p), genome)
    assert df.iloc[0]["gene_id"] == "g1"
    assert df.iloc[0]["gene_family"] == "F1"


def test_placed_column_and_unplaced_mode(tmp_path):
    from gfviewer.genome import load_genome

    gpath = tmp_path / "g.txt"
    gpath.write_text("chr1,1000000\nchr2,900000\nscaffold_7,40000\n")
    g = load_genome(str(gpath))
    csv = tmp_path / "d.csv"
    csv.write_text(
        "gene_id,gene_family,chromosome,start,end,strand\n"
        "g1,F1,chr1,1000,2000,+\n"
        "g2,F1,scaffold_7,100,200,+\n"       # in genome, but classified unplaced
        "g3,F1,ctg_stray,10,20,+\n"          # not in genome at all
    )
    # default: the truly-missing sequence errors
    with pytest.raises(InputValidationError):
        load_features(str(csv), g)

    df, warn = load_features(str(csv), g, on_unknown_chrom="unplaced")
    assert set(df.columns) >= {"placed"}
    placed = dict(zip(df["gene_id"], df["placed"]))
    assert placed["g1"] is True or placed["g1"] == True  # noqa: E712
    assert not placed["g2"]
    assert not placed["g3"]
    assert any("unplaced" in w for w in warn)


def test_multi_file_mixed_formats_use_basename(genome, tmp_path):
    a = tmp_path / "Alpha.bed"
    a.write_text("chr1\t10\t20\tg1\t0\t+\n")
    b = tmp_path / "Beta.gff3"
    b.write_text("chr2\tx\tgene\t30\t40\t.\t+\t.\tID=g2\n")
    df, _ = load_features([str(a), str(b)], genome)
    assert _families(df) == {"Alpha", "Beta"}
