import gzip
import sqlite3
import zipfile

import numpy as np
import pytest

from dna_app.analysis import alt_count, analyse, categorize
from dna_app.ancestry import estimate_admixture
from dna_app.parsers import ParseStats, UnsupportedFormat, parse_file

from .conftest import simulate_genotypes, write_23andme


def _noop(*_):
    pass


def by_gene(result):
    return {f["gene"]: f for f in result["findings"]}


# ------------------------------------------------------------------ builders


def test_clinvar_db_contents(reference_dbs):
    con = sqlite3.connect(reference_dbs / "clinvar.sqlite")
    assert con.execute("SELECT COUNT(*) FROM clinvar").fetchone()[0] == 14  # 7 variants x 2 builds
    stars = dict(con.execute("SELECT gene, stars FROM clinvar WHERE assembly='GRCh37'"))
    assert stars["BRCA1"] == 3 and stars["CFTR"] == 4 and stars["HFE"] == 2


def test_1000g_keeps_clinvar_sites_and_splits_multiallelic(reference_dbs):
    con = sqlite3.connect(reference_dbs / "1000g.sqlite")
    assert con.execute("SELECT COUNT(*) FROM af WHERE rsid='rs80357906'").fetchone()[0] == 1
    assert con.execute("SELECT COUNT(*) FROM af WHERE rsid='rs6'").fetchone()[0] == 0  # rare, not ClinVar
    rows = con.execute("SELECT alt, sas FROM af WHERE rsid='rs5' ORDER BY alt").fetchall()
    assert rows == [("C", 0.2), ("G", 0.5)]


# ------------------------------------------------------------------ helpers


def test_categorize():
    assert categorize("Pathogenic/Likely pathogenic") == "pathogenic"
    assert categorize("Conflicting classifications of pathogenicity") == "conflicting"
    assert categorize("Benign/Likely benign") == "benign"
    assert categorize("drug response") == "drug"
    assert categorize("Uncertain significance") == "uncertain"


def test_alt_count_strand_handling():
    assert alt_count(("A", "G"), "G", "A", allow_flip=True) == 1
    assert alt_count(("C", "T"), "G", "A", allow_flip=True) == 1  # opposite strand
    assert alt_count(("C", "T"), "G", "A", allow_flip=False) is None
    assert alt_count(("T", "T"), "A", "T", allow_flip=True) == 2
    assert alt_count(("C", "C"), "A", "T", allow_flip=True) is None  # palindromic: no flip


def test_admixture_recovers_mixture(kg_sites):
    genos = simulate_genotypes(kg_sites, {"EUR": 0.5, "EAS": 0.5})
    g = np.array([x[-1] for x in genos])
    p = np.array([s[-1] for s in kg_sites])
    q, _ = estimate_admixture(g, p)
    q = dict(zip(("AFR", "AMR", "EAS", "EUR", "SAS"), q))
    assert abs(q["EUR"] - 0.5) < 0.1 and abs(q["EAS"] - 0.5) < 0.1


# ------------------------------------------------------------------ parsers


def test_parses_ancestry_and_myheritage(tmp_path):
    anc = tmp_path / "ancestry.txt"
    anc.write_text("#AncestryDNA raw data\nrsid\tchromosome\tposition\tallele1\tallele2\n"
                   "rs1\t1\t100\tA\tG\nrs2\t23\t200\tC\tC\nrs3\t26\t300\t0\t0\n")
    s = ParseStats()
    out = list(parse_file(str(anc), s))
    assert s.format == "ancestrydna" and s.no_calls == 1
    assert [(v.chrom, v.alleles) for v in out] == [("1", ("A", "G")), ("X", ("C", "C"))]

    mh = tmp_path / "mh.csv"
    mh.write_text('RSID,CHROMOSOME,POSITION,RESULT\n"rs1","1","100","AG"\n"rs4","MT","5","T"\n')
    s = ParseStats()
    out = list(parse_file(str(mh), s))
    assert s.format == "myheritage"
    assert out[1].chrom == "MT" and out[1].alleles == ("T", "T")


def test_parses_gzipped_vcf_with_build_detection(tmp_path):
    vcf = tmp_path / "x.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                 "chr1\t10\trs1\tA\tG,T\t.\tPASS\t.\tGT:DP\t1|2:30\n"
                 "chr1\t11\t.\tA\t.\t.\tPASS\t.\tGT\t0/0\n"
                 "chr1\t12\t.\tC\tG\t.\tPASS\t.\tGT\t./.\n")
    s = ParseStats()
    out = list(parse_file(str(vcf), s))
    assert s.assembly == "GRCh38" and s.sample == "S1" and s.no_calls == 1
    assert out == [("rs1", "1", 10, "A", ("G", "T"))]


def test_rejects_unknown_format(tmp_path):
    p = tmp_path / "junk.txt"
    p.write_text("hello world\n")
    with pytest.raises(UnsupportedFormat):
        list(parse_file(str(p), ParseStats()))


# ------------------------------------------------------------------ end to end


@pytest.mark.parametrize("wrap", ["plain", "zip"])
def test_23andme_end_to_end(configured, kg_sites, wrap):
    genos = simulate_genotypes(kg_sites, {"AFR": 0.8, "EUR": 0.2})
    path = configured / "genome.txt"
    write_23andme(path, kg_sites, genos)
    if wrap == "zip":
        zpath = configured / "genome.zip"
        with zipfile.ZipFile(zpath, "w") as z:
            z.write(path, "genome_v5_Full.txt")
        path = zpath
    r = analyse(str(path), configured / "work", _noop)

    assert r["file"]["format"] == "23andme" and r["file"]["assembly"] == "GRCh37"
    f = by_gene(r)
    assert f["BRCA1"]["zygosity"] == "heterozygous" and f["BRCA1"]["category"] == "pathogenic"
    assert f["HFE"]["zygosity"] == "homozygous"
    assert f["CYP2C19"]["category"] == "drug" and f["CYP2C19"]["opposite_strand"]
    assert not f["BRCA1"]["opposite_strand"]
    assert "APOE" not in f and "TTN" not in f and "CFTR" not in f
    assert r["clinvar"]["carried"]["benign"] == 1
    assert r["clinvar"]["checked"]["risk"] == 1
    assert f["HFE"]["frequencies"]["EUR"] == pytest.approx(0.05)

    anc = {p["code"]: p["proportion"] for p in r["ancestry"]["populations"]}
    assert r["ancestry"]["available"]
    assert anc["AFR"] == pytest.approx(0.8, abs=0.1)
    assert anc["EUR"] == pytest.approx(0.2, abs=0.1)


def test_vcf_grch38_matches_by_position_and_rsid(configured, kg_sites):
    genos = simulate_genotypes(kg_sites, {"SAS": 1.0})
    path = configured / "sample.vcf"
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n##reference=GRCh38\n"
                 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tME\n")
        fh.write("chr7\t117559590\t.\tATCT\tA\t.\tPASS\t.\tGT\t0/1\n")   # CFTR F508del, het
        fh.write("chr17\t43094464\t.\tC\tT,G\t.\tPASS\t.\tGT\t1/2\n")   # BRCA1, multi-allelic
        for chrom, pos, rsid, ref, alt, n in genos:  # GRCh38 VCF: ancestry joins by rsID
            gt = {0: "0/0", 1: "0/1", 2: "1/1"}[n]
            fh.write(f"{chrom}\t{pos + 5}\t{rsid}\t{ref}\t{alt}\t.\tPASS\t.\tGT\t{gt}\n")
    r = analyse(str(path), configured / "work", _noop)
    f = by_gene(r)
    assert r["file"]["assembly"] == "GRCh38"
    assert f["CFTR"]["zygosity"] == "heterozygous" and f["CFTR"]["stars"] == 4
    assert f["BRCA1"]["zygosity"] == "heterozygous"
    assert f["BRCA1"]["frequencies"]["EUR"] == pytest.approx(0.001)  # via rsID
    top = r["ancestry"]["populations"][0]
    assert top["code"] == "SAS" and top["proportion"] > 0.85


def test_without_reference_databases(configured, tmp_path, monkeypatch):
    from dna_app import config
    monkeypatch.setattr(config, "CLINVAR_DB", tmp_path / "missing1.sqlite")
    monkeypatch.setattr(config, "KG_DB", tmp_path / "missing2.sqlite")
    p = tmp_path / "g.txt"
    p.write_text("rs1\t1\t100\tAG\n")
    r = analyse(str(p), tmp_path / "work", _noop)
    assert not r["clinvar"]["available"] and not r["ancestry"]["available"]
    assert len(r["notes"]) == 2


def test_clinvar_keys(reference_dbs):
    from dna_app.build_db import ClinvarKeys
    keys = ClinvarKeys(reference_dbs / "clinvar.sqlite")
    assert keys.contains("rs80357906", "1", 1)
    assert keys.contains(None, "17", 41246481)
    assert not keys.contains("rs999999999", "17", 41246482)
    assert not keys.contains(None, "17", 43094464)  # GRCh38 position: 1000 Genomes is GRCh37
