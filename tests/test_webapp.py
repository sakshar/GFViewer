import io
import os
import time

import pytest

from gfviewer_web import create_app


@pytest.fixture
def app(tmp_path):
    app = create_app({"DATA_DIR": str(tmp_path / "jobs"), "WORKERS": 1,
                       "TESTING": True, "SECRET_KEY": "test"})
    return app


@pytest.fixture
def client(app):
    return app.test_client()


def _wait(client, job_id, timeout=60):
    end = time.time() + timeout
    while time.time() < end:
        r = client.get("/api/jobs/{}/status".format(job_id))
        data = r.get_json()
        if data["state"] in ("done", "error"):
            return data
        time.sleep(0.3)
    raise AssertionError("job did not finish")


def test_health(client):
    assert client.get("/api/health").get_json() == {"status": "ok"}


def test_home_and_help(client):
    assert client.get("/").status_code == 200
    assert client.get("/help").status_code == 200
    assert client.get("/about_us").status_code == 200


def test_install_page_and_guide(client):
    home = client.get("/")
    assert b"pip install gfviewer" in home.data
    assert b'href="/install"' in home.data

    page = client.get("/install")
    assert page.status_code == 200
    assert b"pip install gfviewer" in page.data
    assert b"gfviewer --help" in page.data
    assert b"pypi.org/project/gfviewer" in page.data

    md = client.get("/install/guide.md")
    assert md.status_code == 200
    assert md.mimetype == "text/markdown"
    assert md.data.startswith(b"# GFViewer")

    assert b'href="/install"' in client.get("/help").data


def _files(data_path, genome_path):
    return {
        "data": (open(data_path, "rb"), os.path.basename(data_path)),
        "genome": (open(genome_path, "rb"), os.path.basename(genome_path)),
    }


def test_job_lifecycle_and_downloads(client):
    root = os.path.join(os.path.dirname(__file__), "..", "static", "tests")
    r = client.post(
        "/api/jobs",
        data={
            **_files(os.path.join(root, "data_test_2.csv"),
                     os.path.join(root, "chrs_test_2.txt")),
            "with_analytics": "1",
            "permutations": "50",
            "formats": "png",
        },
        content_type="multipart/form-data",
    )
    assert r.status_code == 202
    job_id = r.get_json()["job_id"]

    data = _wait(client, job_id)
    assert data["state"] == "done", data.get("error")
    assert set(data["families"]) == {"MGF1", "MGF2", "MGF3", "MGF4", "MGF5", "MGF6"}
    assert data["chromosomes"] == ["chr1", "chr2", "chr3"]

    svg = client.get("/api/jobs/{}/svg".format(job_id))
    assert svg.status_code == 200 and b"gfv-legend" in svg.data

    for fmt in ("pdf", "png", "svg"):
        d = client.get("/api/jobs/{}/download/{}".format(job_id, fmt))
        assert d.status_code == 200 and len(d.data) > 0

    a = client.get("/api/jobs/{}/analytics/analytics_summary.json".format(job_id))
    assert a.status_code == 200 and b"n_families" in a.data

    gpf = client.get(
        "/api/jobs/{}/analytics/analytics_genes_per_family.csv".format(job_id)
    )
    assert gpf.status_code == 200 and gpf.data.startswith(b"gene_family,")

    for fmt in ("pdf", "svg", "png", "jpg", "eps"):
        c = client.get("/api/jobs/{}/chart/{}".format(job_id, fmt))
        assert c.status_code == 200 and len(c.data) > 0

    # new analytics tables + charts (these are produced for any >= 2-family run)
    for name in ("analytics_duplication_modes.csv", "analytics_chromosome_enrichment.csv",
                 "analytics_strand_bias.csv", "analytics_chromosome_richness.csv",
                 "analytics_positional_profile.csv", "analytics_family_proximity.csv"):
        r = client.get("/api/jobs/{}/analytics/{}".format(job_id, name))
        assert r.status_code == 200 and len(r.data) > 0
    for which in ("positional_profile", "family_proximity"):
        r = client.get("/api/jobs/{}/chart/png?name={}".format(job_id, which))
        assert r.status_code == 200 and len(r.data) > 0
    assert client.get("/api/jobs/{}/chart/png?name=bogus".format(job_id)).status_code == 400

    # bundle zips
    import zipfile

    b = client.get("/api/jobs/{}/bundle".format(job_id))
    assert b.status_code == 200
    zf = zipfile.ZipFile(io.BytesIO(b.data))
    names = zf.namelist()
    assert "analytics_summary.json" in names and any(n.startswith("gfviewer.") for n in names)

    ba = client.get("/api/jobs/{}/bundle?all=1".format(job_id))
    assert ba.status_code == 200
    zfa = zipfile.ZipFile(io.BytesIO(ba.data))
    na = zfa.namelist()
    assert "gfviewer.pdf" in na and "gfviewer.eps" in na
    assert "analytics_positional_profile.png" in na

    page = client.get("/jobs/{}".format(job_id))
    assert page.status_code == 200


def test_restyle(client):
    root = os.path.join(os.path.dirname(__file__), "..", "static", "tests")
    r = client.post(
        "/api/jobs",
        data=_files(os.path.join(root, "data_test_2.csv"),
                    os.path.join(root, "chrs_test_2.txt")),
        content_type="multipart/form-data",
    )
    job_id = r.get_json()["job_id"]
    assert _wait(client, job_id)["state"] == "done"

    resp = client.post(
        "/api/jobs/{}/restyle".format(job_id),
        json={"tick_style": "lollipop", "label_mode": "family",
              "orientation": "vertical", "only_chromosomes": ["chr1"],
              "only_families": ["MGF1", "MGF3"],
              "_colors": {"MGF1": [0.1, 0.2, 0.3]}},
    )
    assert resp.status_code == 200
    style = resp.get_json()["style"]
    assert style["tick_style"] == "lollipop"
    assert style["orientation"] == "vertical"
    assert style["only_chromosomes"] == ["chr1"]
    assert style["only_families"] == ["MGF1", "MGF3"]
    svg = client.get("/api/jobs/{}/svg".format(job_id))
    assert svg.status_code == 200 and b"gfv-legend" in svg.data
    assert b"gfv-chrom-chr1" in svg.data and b"gfv-chrom-chr2" not in svg.data
    assert b"gfv-fam-MGF1" in svg.data and b"gfv-fam-MGF2" not in svg.data

    # the restyle must persist: a fresh status read and an on-demand format
    # download both have to reflect the subset, not the original style
    st = client.get("/api/jobs/{}/status".format(job_id)).get_json()
    assert st["style"]["only_chromosomes"] == ["chr1"]
    assert st["style"]["only_families"] == ["MGF1", "MGF3"]
    assert st["style"]["tick_style"] == "lollipop"
    png = client.get("/api/jobs/{}/download/png".format(job_id))
    assert png.status_code == 200 and len(png.data) > 0

    # clearing the subsets restores every chromosome and family
    resp = client.post(
        "/api/jobs/{}/restyle".format(job_id),
        json={"only_chromosomes": [], "only_families": []},
    )
    assert resp.status_code == 200
    svg = client.get("/api/jobs/{}/svg".format(job_id))
    assert b"gfv-chrom-chr2" in svg.data and b"gfv-chrom-chr3" in svg.data
    assert b"gfv-fam-MGF2" in svg.data and b"gfv-fam-MGF5" in svg.data
    assert client.get("/api/jobs/{}/status".format(job_id)).get_json()["style"][
        "only_chromosomes"
    ] == []


def test_invalid_upload_returns_json_error(client):
    root = os.path.join(os.path.dirname(__file__), "..", "static", "tests")
    r = client.post(
        "/api/jobs",
        data={
            "data": (io.BytesIO(b"not,a,valid\nfile\n"), "junk.csv"),
            "genome": (open(os.path.join(root, "chrs_test_2.txt"), "rb"), "chrs_test_2.txt"),
        },
        content_type="multipart/form-data",
    )
    # submission accepted, job fails during processing with a clean message
    assert r.status_code == 202
    data = _wait(client, r.get_json()["job_id"])
    assert data["state"] == "error"
    assert "column" in data["error"].lower() or "family" in data["error"].lower()


def test_unknown_job_404(client):
    assert client.get("/api/jobs/nope/status").status_code == 404


def test_example_datasets(client):
    home = client.get("/")
    assert home.status_code == 200
    assert b"Example datasets" in home.data and b"celegans_20" in home.data

    import zipfile

    z = client.get("/datasets/arabidopsis_10/download")
    assert z.status_code == 200
    names = zipfile.ZipFile(io.BytesIO(z.data)).namelist()
    assert any(n.endswith("genes.tsv") for n in names)
    assert any(n.endswith("genome.txt") for n in names)

    allz = client.get("/datasets/download-all")
    assert allz.status_code == 200
    an = zipfile.ZipFile(io.BytesIO(allz.data)).namelist()
    assert any("synthetic/celegans_20/genes.tsv" in n for n in an)
    assert any(n.endswith("formats/genes.gff3") for n in an)

    assert client.get("/datasets/nope/download").status_code == 404

    # "Run" submits a real job and redirects to its results page
    r = client.get("/datasets/arabidopsis_10/run")
    assert r.status_code == 302
    job_id = r.headers["Location"].rstrip("/").split("/")[-1]
    data = _wait(client, job_id)
    assert data["state"] == "done", data.get("error")
    assert len(data["families"]) == 10


def test_show_titles_toggle_via_web(client):
    root = os.path.join(os.path.dirname(__file__), "..", "static", "tests")
    r = client.post(
        "/api/jobs",
        data={
            **_files(os.path.join(root, "data_test_2.csv"),
                     os.path.join(root, "chrs_test_2.txt")),
            "with_analytics": "1", "permutations": "30", "formats": "png",
            "title": "Babesia MO1",
            # show_titles checkbox left unchecked -> absent from the POST
        },
        content_type="multipart/form-data",
    )
    job_id = r.get_json()["job_id"]
    assert _wait(client, job_id)["state"] == "done"

    st = client.get("/api/jobs/{}/status".format(job_id)).get_json()
    assert st["style"]["show_titles"] is False

    fig = client.get("/api/jobs/{}/download/svg".format(job_id))
    assert b"Babesia MO1" not in fig.data and b"gfv-title" not in fig.data
    chart = client.get("/api/jobs/{}/chart/svg?name=positional_profile".format(job_id))
    assert b"Positional density profile" not in chart.data

    # flip it back on through the editor's restyle path
    client.post("/api/jobs/{}/restyle".format(job_id), json={"show_titles": True})
    fig = client.get("/api/jobs/{}/download/svg".format(job_id))
    assert b"Babesia MO1" in fig.data
    chart = client.get("/api/jobs/{}/chart/svg?name=positional_profile".format(job_id))
    assert b"Positional density profile" in chart.data


def test_usage_monitor(client):
    # a couple of page views + the stats endpoints
    client.get("/")
    client.get("/help")
    api = client.get("/api/stats")
    assert api.status_code == 200
    s = api.get_json()
    assert s["totals"]["page_views"] >= 2
    assert s["totals"]["visits"] >= 1
    assert any(p["path"] == "/help" for p in s["pages"])

    page = client.get("/stats")
    assert page.status_code == 200 and b"unique visitors" in page.data

    # a completed job bumps the job counters
    root = os.path.join(os.path.dirname(__file__), "..", "static", "tests")
    r = client.post(
        "/api/jobs",
        data=_files(os.path.join(root, "data_test_2.csv"),
                    os.path.join(root, "chrs_test_2.txt")),
        content_type="multipart/form-data",
    )
    assert _wait(client, r.get_json()["job_id"])["state"] == "done"
    t = client.get("/api/stats").get_json()["totals"]
    assert t["jobs_submitted"] >= 1 and t["jobs_completed"] >= 1

    # example-dataset download is counted
    before = client.get("/api/stats").get_json()["totals"]["downloads"]
    client.get("/datasets/arabidopsis_10/download")
    assert client.get("/api/stats").get_json()["totals"]["downloads"] == before + 1


def test_usage_stats_token_gate(tmp_path):
    app = create_app({"DATA_DIR": str(tmp_path / "jobs"), "WORKERS": 1,
                      "SECRET_KEY": "test", "STATS_TOKEN": "s3cret"})
    c = app.test_client()
    assert c.get("/api/stats").status_code == 403
    assert c.get("/stats").status_code == 403
    assert c.get("/api/stats?token=s3cret").status_code == 200
    assert c.get("/stats?token=s3cret").status_code == 200
