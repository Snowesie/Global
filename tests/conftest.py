"""Synthetic ClinVar and 1000 Genomes files in the real upstream formats."""

from __future__ import annotations

import gzip

import numpy as np
import pytest

from dna_app import build_db, config

CLINVAR_HEADER = (
    "#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance\tClinSigSimple\t"
    "LastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDS\tPhenotypeList\tOrigin\t"
    "OriginSimple\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\t"
    "AlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\tGuidelines\tTestedInGTR\tOtherIDs\t"
    "SubmitterCategories\tVariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF"
)

# (variation_id, gene, significance, review, rsid, chrom, pos37, pos38, ref, alt, phenotype)
CLINVAR_ROWS = [
    (17661, "BRCA1", "Pathogenic", "reviewed by expert panel", 80357906, "17", 41246481, 43094464, "C", "T",
     "Hereditary breast ovarian cancer syndrome"),
    (10, "HFE", "Pathogenic/Likely pathogenic", "criteria provided, multiple submitters, no conflicts",
     1800562, "6", 26093141, 26092913, "G", "A", "Hereditary hemochromatosis"),
    (17000, "APOE", "risk factor", "criteria provided, single submitter", 429358, "19", 45411941, 44908684,
     "T", "C", "Alzheimer disease"),
    (37000, "CYP2C19", "drug response", "reviewed by expert panel", 4244285, "10", 96541616, 94781859,
     "G", "A", "Clopidogrel response"),
    (5000, "TTN", "Benign", "criteria provided, single submitter", 1001, "2", 179400000, 178535273,
     "A", "G", "not provided"),
    (6000, "MYH7", "Uncertain significance", "criteria provided, single submitter", 1002, "14", 23900000,
     23430791, "C", "G", "Hypertrophic cardiomyopathy"),
    (7000, "CFTR", "Pathogenic", "practice guideline", 113993960, "7", 117199644, 117559590, "ATCT", "A",
     "Cystic fibrosis"),
]


def write_clinvar(path):
    with gzip.open(path, "wt") as fh:
        fh.write(CLINVAR_HEADER + "\n")
        for vid, gene, sig, review, rs, chrom, p37, p38, ref, alt, pheno in CLINVAR_ROWS:
            for asm, pos in (("GRCh37", p37), ("GRCh38", p38)):
                row = [str(vid + 1), "single nucleotide variant", f"NM_000000.0:c.1{ref}>{alt}", "1", gene,
                       "HGNC:1", sig, "1", "Jan 01, 2024", str(rs), "-", "RCV0", "MedGen:C0", pheno,
                       "germline", "germline", asm, "NC_0", chrom, str(pos), str(pos), ref, alt, "1p1",
                       review, "3", "-", "N", "-", "3", str(vid), str(pos), ref, alt]
                fh.write("\t".join(row) + "\n")
        # A row ClinVar ships without VCF coordinates must be ignored.
        fh.write("\t".join(["1", "copy number loss", "x", "1", "X", "-", "Pathogenic", "1", "-", "-1", "-",
                            "-", "-", "-", "germline", "germline", "GRCh37", "-", "1", "1", "2", "na", "na",
                            "-", "no assertion criteria provided", "1", "-", "N", "-", "1", "999", "-1",
                            "na", "na"]) + "\n")


N_ANCESTRY_SITES = 4000
POPS = config.SUPERPOPS


def make_kg_sites(seed=1):
    """Random informative SNPs: chrom 1-22, positions far from the ClinVar sites."""
    rng = np.random.default_rng(seed)
    sites = []
    for i in range(N_ANCESTRY_SITES):
        chrom = str(1 + i % 22)
        pos = 1_000_000 + i * 997
        ref, alt = [("A", "G"), ("C", "T"), ("G", "T"), ("A", "C")][i % 4]
        f = rng.beta(0.6, 0.6, size=5)
        sites.append((chrom, pos, f"rs9{i:07d}", ref, alt, f))
    return sites


def write_kg(path, sites):
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for chrom, pos, rsid, ref, alt, f in sites:
            info = f"AF={f.mean():.4f};" + ";".join(f"{p}_AF={v:.4f}" for p, v in zip(POPS, f))
            fh.write(f"{chrom}\t{pos}\t{rsid}\t{ref}\t{alt}\t100\tPASS\t{info}\n")
        # ClinVar sites (rare, so only kept because they are in ClinVar) and a multi-allelic site.
        fh.write("17\t41246481\trs80357906\tC\tT\t100\tPASS\tAF=0.0002;AFR_AF=0;AMR_AF=0;EAS_AF=0;EUR_AF=0.001;SAS_AF=0\n")
        fh.write("6\t26093141\trs1800562\tG\tA\t100\tPASS\tAF=0.01;AFR_AF=0;AMR_AF=0.01;EAS_AF=0;EUR_AF=0.05;SAS_AF=0.01\n")
        fh.write("3\t500\trs5\tA\tC,G\t100\tPASS\tAF=0.2,0.3;AFR_AF=0.1,0.2;AMR_AF=0.2,0.3;EAS_AF=0.3,0.4;EUR_AF=0.2,0.1;SAS_AF=0.2,0.5\n")
        fh.write("3\t600\trs6\tA\tC\t100\tPASS\tAF=0.0001;AFR_AF=0;AMR_AF=0;EAS_AF=0;EUR_AF=0;SAS_AF=0.0005\n")


@pytest.fixture(scope="session")
def kg_sites():
    return make_kg_sites()


@pytest.fixture(scope="session")
def reference_dbs(tmp_path_factory, kg_sites):
    d = tmp_path_factory.mktemp("ref")
    write_clinvar(d / "variant_summary.txt.gz")
    write_kg(d / "sites.vcf.gz", kg_sites)
    build_db.build_clinvar(str(d / "variant_summary.txt.gz"), d / "clinvar.sqlite")
    build_db.build_1000g(str(d / "sites.vcf.gz"), d / "1000g.sqlite", 0.01, d / "clinvar.sqlite")
    return d


@pytest.fixture
def configured(reference_dbs, tmp_path, monkeypatch):
    monkeypatch.setattr(config, "CLINVAR_DB", reference_dbs / "clinvar.sqlite")
    monkeypatch.setattr(config, "KG_DB", reference_dbs / "1000g.sqlite")
    monkeypatch.setattr(config, "UPLOAD_DIR", tmp_path / "uploads")
    monkeypatch.setattr(config, "RESULTS_DIR", tmp_path / "results")
    return tmp_path


def simulate_genotypes(sites, mixture, seed=2):
    """Draw alt-allele counts for an individual with the given superpop mixture."""
    rng = np.random.default_rng(seed)
    q = np.array([mixture.get(p, 0.0) for p in POPS])
    out = []
    for chrom, pos, rsid, ref, alt, f in sites:
        p = float(f @ q)
        out.append((chrom, pos, rsid, ref, alt, int(rng.binomial(2, p))))
    return out


CLINVAR_GENOTYPES_23 = [
    ("rs80357906", "17", 41246481, "CT"),  # BRCA1 het pathogenic
    ("rs1800562", "6", 26093141, "AA"),    # HFE hom pathogenic
    ("rs429358", "19", 45411941, "TT"),    # APOE: does not carry risk allele
    ("rs4244285", "10", 96541616, "CT"),   # CYP2C19 reported on minus strand: G/A -> C/T
    ("rs1001", "2", 179400000, "AG"),      # benign, counted but not listed
    ("rs1002", "14", 23900000, "--"),      # no-call
]


def write_23andme(path, sites, genos):
    with open(path, "w") as fh:
        fh.write("# This data file generated by 23andMe\n# rsid\tchromosome\tposition\tgenotype\n")
        for rsid, chrom, pos, gt in CLINVAR_GENOTYPES_23:
            fh.write(f"{rsid}\t{chrom}\t{pos}\t{gt}\n")
        for chrom, pos, rsid, ref, alt, n in genos:
            gt = {0: ref + ref, 1: ref + alt, 2: alt + alt}[n]
            fh.write(f"{rsid}\t{chrom}\t{pos}\t{gt}\n")
