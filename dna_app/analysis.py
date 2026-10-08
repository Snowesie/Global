"""Annotate an uploaded genome with ClinVar and 1000 Genomes.

The upload is streamed once into a per-job SQLite table, then joined against
the reference databases (attached read-only) so even a 1 GB whole-genome VCF
never has to sit in memory.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Callable

import numpy as np

from . import config
from .ancestry import bootstrap_ci, estimate_admixture
from .parsers import ParseStats, parse_file

Progress = Callable[[str, float], None]

COMPLEMENT = str.maketrans("ACGT", "TGCA")
PALINDROMIC = {frozenset("AT"), frozenset("CG")}

# Order matters: the first matching rule names the finding's category.
CATEGORIES = [
    ("conflicting", "Conflicting classifications"),
    ("pathogenic", "Pathogenic / likely pathogenic"),
    ("risk", "Risk factor"),
    ("drug", "Drug response"),
    ("protective", "Protective"),
    ("association", "Association / trait"),
    ("uncertain", "Uncertain significance"),
    ("benign", "Benign / likely benign"),
    ("other", "Other / not provided"),
]
LISTED_CATEGORIES = {"pathogenic", "risk", "drug", "protective", "association", "conflicting", "uncertain"}


def categorize(significance: str) -> str:
    s = significance.lower()
    if "conflicting" in s:
        return "conflicting"
    if "pathogenic" in s and "benign" not in s:
        return "pathogenic"
    if "risk factor" in s:
        return "risk"
    if "drug response" in s:
        return "drug"
    if "protective" in s:
        return "protective"
    if "association" in s or "affects" in s:
        return "association"
    if "uncertain" in s:
        return "uncertain"
    if "benign" in s:
        return "benign"
    return "other"


def alt_count(alleles: tuple[str, str], ref: str, alt: str, allow_flip: bool) -> int | None:
    """How many copies of ``alt`` the genotype carries, or None if alleles don't fit.

    Array data carries no reference allele, so when the called bases are not
    the site's ref/alt we try the opposite strand (never for A/T or C/G SNPs,
    where a flip is undetectable).
    """
    pair = {ref, alt}
    if all(a in pair for a in alleles):
        return sum(a == alt for a in alleles)
    if allow_flip and len(ref) == 1 and len(alt) == 1 and frozenset(pair) not in PALINDROMIC:
        flipped = tuple(a.translate(COMPLEMENT) for a in alleles)
        if all(a in pair for a in flipped):
            return sum(a == alt for a in flipped)
    return None


def _db_meta(path: Path) -> dict | None:
    if not path.exists():
        return None
    con = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return dict(con.execute("SELECT key, value FROM meta"))
    except sqlite3.Error:
        return {}
    finally:
        con.close()


def database_status() -> dict:
    return {"clinvar": _db_meta(config.CLINVAR_DB), "kg": _db_meta(config.KG_DB)}


def _load_user_table(con: sqlite3.Connection, path: str, stats: ParseStats,
                     assembly: str | None, progress: Progress) -> None:
    con.execute(
        "CREATE TABLE user (rsid TEXT, chrom TEXT, pos INTEGER, ref TEXT, a1 TEXT, a2 TEXT)"
    )
    batch = []
    for v in parse_file(path, stats, assembly):
        batch.append((v.rsid, v.chrom, v.pos, v.ref, v.alleles[0], v.alleles[1]))
        if len(batch) >= 100_000:
            con.executemany("INSERT INTO user VALUES (?,?,?,?,?,?)", batch)
            batch.clear()
            progress(f"Reading file — {stats.variants:,} genotypes", 0.05)
    con.executemany("INSERT INTO user VALUES (?,?,?,?,?,?)", batch)
    progress("Indexing genotypes", 0.35)
    con.execute("CREATE INDEX user_pos ON user (chrom, pos)")
    con.execute("CREATE INDEX user_rsid ON user (rsid)")
    con.commit()


def _clinvar_findings(con: sqlite3.Connection, stats: ParseStats) -> tuple[list[dict], dict]:
    is_array = stats.format != "vcf"
    # Positional join on the file's assembly, plus rsID join for arrays where
    # the build may differ or coordinates have shifted between dbSNP releases.
    rows = con.execute(
        """
        SELECT u.rsid, u.chrom, u.pos, u.ref, u.a1, u.a2, c.*
        FROM user u JOIN cv.clinvar c
          ON c.assembly = ? AND c.chrom = u.chrom AND c.pos = u.pos
        """,
        (stats.assembly,),
    ).fetchall()
    if is_array:
        rows += con.execute(
            """
            SELECT u.rsid, u.chrom, u.pos, u.ref, u.a1, u.a2, c.*
            FROM user u JOIN cv.clinvar c
              ON c.rsid = u.rsid AND c.assembly = ?
            WHERE u.rsid IS NOT NULL AND NOT (c.chrom = u.chrom AND c.pos = u.pos)
            """,
            (stats.assembly,),
        ).fetchall()
    cols = ["u_rsid", "u_chrom", "u_pos", "u_ref", "a1", "a2"] + [
        d[0] for d in con.execute("SELECT * FROM cv.clinvar LIMIT 0").description
    ]
    seen: set[int] = set()
    findings: list[dict] = []
    counts = {key: 0 for key, _ in CATEGORIES}
    checked = {key: 0 for key, _ in CATEGORIES}
    for row in rows:
        r = dict(zip(cols, row))
        if r["variation_id"] in seen:
            continue
        if is_array and (len(r["ref"]) != 1 or len(r["alt"]) != 1):
            continue  # arrays can't genotype indels reliably
        if not is_array and r["u_ref"] != r["ref"]:
            continue
        if is_array:
            n_alt = alt_count((r["a1"], r["a2"]), r["ref"], r["alt"], allow_flip=True)
        else:  # VCF alleles are explicit, and may include a second, unrelated alt
            n_alt = (r["a1"] == r["alt"]) + (r["a2"] == r["alt"])
        if n_alt is None:
            continue
        seen.add(r["variation_id"])
        cat = categorize(r["significance"])
        checked[cat] += 1
        if n_alt == 0:
            continue
        counts[cat] += 1
        if cat not in LISTED_CATEGORIES:
            continue
        findings.append({
            "category": cat,
            "gene": r["gene"] or None,
            "name": r["name"],
            "type": r["type"],
            "rsid": r["rsid"] or r["u_rsid"],
            "chrom": r["chrom"],
            "pos": r["pos"],
            "ref": r["ref"],
            "alt": r["alt"],
            "genotype": f"{r['a1']}/{r['a2']}",
            "zygosity": "homozygous" if n_alt == 2 else "heterozygous",
            "opposite_strand": is_array and not {r["a1"], r["a2"]} <= {r["ref"], r["alt"]},
            "significance": r["significance"],
            "review_status": r["review_status"],
            "stars": r["stars"],
            "conditions": [p for p in (r["phenotypes"] or "").split("|") if p and p != "not provided"],
            "origin": r["origin"],
            "last_evaluated": r["last_evaluated"],
            "variation_id": r["variation_id"],
            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{r['variation_id']}/",
        })
    order = {key: i for i, (key, _) in enumerate(CATEGORIES)}
    findings.sort(key=lambda f: (order[f["category"]], -f["stars"], f["gene"] or "~"))
    return findings, {"carried": counts, "checked": checked}


def _attach_frequencies(con: sqlite3.Connection, findings: list[dict], assembly: str) -> None:
    for f in findings:
        if assembly == "GRCh37":
            row = con.execute(
                "SELECT af, afr, amr, eas, eur, sas FROM kg.af WHERE chrom=? AND pos=? AND ref=? AND alt=?",
                (f["chrom"], f["pos"], f["ref"], f["alt"]),
            ).fetchone()
        else:
            row = None
        if row is None and f["rsid"]:
            row = con.execute(
                "SELECT af, afr, amr, eas, eur, sas FROM kg.af WHERE rsid=? AND alt=?",
                (f["rsid"], f["alt"]),
            ).fetchone()
        f["frequencies"] = (
            dict(zip(("ALL",) + config.SUPERPOPS, row)) if row else None
        )


MAX_ANCESTRY_SITES = 200_000
MIN_ANCESTRY_SITES = 1_000


def _ancestry(con: sqlite3.Connection, stats: ParseStats, progress: Progress) -> dict:
    is_array = stats.format != "vcf"
    autosomes = tuple(str(i) for i in range(1, 23))
    ph = ",".join("?" * len(autosomes))
    informative = "(MAX(k.afr,k.amr,k.eas,k.eur,k.sas) - MIN(k.afr,k.amr,k.eas,k.eur,k.sas)) >= 0.1"
    snp = "length(k.ref)=1 AND length(k.alt)=1"
    wgs = (not is_array) and stats.variants >= 1_500_000
    note = None
    if stats.assembly == "GRCh37":
        join = "u.chrom = k.chrom AND u.pos = k.pos"
    else:
        join = "u.rsid = k.rsid"
    if wgs:
        # Whole-genome VCFs list only variant sites: a common site that is
        # absent from the file is homozygous reference, so sample reference
        # sites and left-join the genotypes.
        total = con.execute(f"SELECT COUNT(*) FROM kg.af k WHERE {snp} AND {informative}").fetchone()[0]
        step = max(1, total // MAX_ANCESTRY_SITES)
        sql = f"""
            SELECT u.ref, u.a1, u.a2, k.ref, k.alt, k.afr, k.amr, k.eas, k.eur, k.sas
            FROM kg.af k LEFT JOIN user u ON {join}
            WHERE k.chrom IN ({ph}) AND {snp} AND {informative} AND (k.rowid % ?) = 0
        """
        rows = con.execute(sql, (*autosomes, step)).fetchall()
        note = "Whole-genome VCF: sites absent from the file were treated as homozygous reference."
    else:
        sql = f"""
            SELECT u.ref, u.a1, u.a2, k.ref, k.alt, k.afr, k.amr, k.eas, k.eur, k.sas
            FROM user u JOIN kg.af k ON {join}
            WHERE u.chrom IN ({ph}) AND {snp} AND {informative}
        """
        rows = con.execute(sql, autosomes).fetchall()
        if not is_array:
            note = ("This VCF looks like an exome or panel; ancestry uses only the sites present "
                    "in the file and may be biased toward populations where those variants are common.")
    progress("Estimating ancestry", 0.8)
    genos, freqs = [], []
    for u_ref, a1, a2, ref, alt, *f in rows:
        if any(x is None for x in f):
            continue
        if a1 is None:
            n = 0
        else:
            if frozenset((ref, alt)) in PALINDROMIC and is_array:
                continue
            if is_array:
                n = alt_count((a1, a2), ref, alt, allow_flip=True)
            elif u_ref != ref:
                continue
            else:
                n = (a1 == alt) + (a2 == alt)
            if n is None:
                continue
        genos.append(n)
        freqs.append(f)
    if len(genos) < MIN_ANCESTRY_SITES:
        return {
            "available": False,
            "sites_used": len(genos),
            "reason": f"Only {len(genos):,} informative 1000 Genomes sites overlapped your data "
                      f"(need {MIN_ANCESTRY_SITES:,}).",
        }
    g = np.array(genos, dtype=np.int8)
    p = np.array(freqs, dtype=np.float64)
    if len(g) > MAX_ANCESTRY_SITES:
        idx = np.random.default_rng(0).choice(len(g), MAX_ANCESTRY_SITES, replace=False)
        g, p = g[idx], p[idx]
    q, ll = estimate_admixture(g, p)
    ci = bootstrap_ci(g, p)
    pops = [
        {
            "code": code,
            "name": config.SUPERPOP_NAMES[code],
            "proportion": float(q[i]),
            "ci_low": float(ci[0, i]),
            "ci_high": float(ci[1, i]),
        }
        for i, code in enumerate(config.SUPERPOPS)
    ]
    pops.sort(key=lambda x: -x["proportion"])
    return {
        "available": True,
        "sites_used": int(len(g)),
        "log_likelihood": ll,
        "populations": pops,
        "note": note,
    }


def analyse(upload_path: str, work_dir: Path, progress: Progress,
            assembly: str | None = None) -> dict:
    work_dir.mkdir(parents=True, exist_ok=True)
    tmp_db = work_dir / "genotypes.sqlite"
    if tmp_db.exists():
        tmp_db.unlink()
    con = sqlite3.connect(f"file:{tmp_db}", uri=True)  # uri=True so ATTACH can open read-only
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    stats = ParseStats()
    notes: list[str] = []
    try:
        progress("Reading file", 0.02)
        _load_user_table(con, upload_path, stats, assembly, progress)
        if stats.variants == 0:
            raise ValueError("No genotype calls were found in this file.")
        notes += stats.notes

        dbs = database_status()
        result: dict = {
            "file": {
                "format": stats.format,
                "assembly": stats.assembly,
                "sample": stats.sample,
                "lines": stats.lines,
                "genotypes": stats.variants,
                "no_calls": stats.no_calls,
                "skipped": stats.skipped,
            },
            "databases": dbs,
        }
        chroms = dict(con.execute("SELECT chrom, COUNT(*) FROM user GROUP BY chrom"))
        result["file"]["chromosomes"] = chroms
        if stats.format != "vcf" and stats.no_calls:
            result["file"]["call_rate"] = stats.variants / (stats.variants + stats.no_calls)

        findings: list[dict] = []
        if dbs["clinvar"] is not None:
            progress("Matching against ClinVar", 0.45)
            con.execute("ATTACH DATABASE ? AS cv", (f"file:{config.CLINVAR_DB}?mode=ro",))
            findings, counts = _clinvar_findings(con, stats)
            result["clinvar"] = {"available": True, **counts}
        else:
            result["clinvar"] = {"available": False}
            notes.append("ClinVar database not installed — run `python -m dna_app.build_db clinvar`.")

        if dbs["kg"] is not None:
            con.execute("ATTACH DATABASE ? AS kg", (f"file:{config.KG_DB}?mode=ro",))
            progress("Looking up 1000 Genomes allele frequencies", 0.65)
            _attach_frequencies(con, findings, stats.assembly)
            result["ancestry"] = _ancestry(con, stats, progress)
        else:
            result["ancestry"] = {"available": False,
                                  "reason": "1000 Genomes database not installed."}
            notes.append("1000 Genomes database not installed — run `python -m dna_app.build_db 1000g`.")

        result["findings"] = findings
        result["categories"] = [{"key": k, "label": label} for k, label in CATEGORIES]
        result["notes"] = notes
        progress("Done", 1.0)
        return result
    finally:
        con.close()
        tmp_db.unlink(missing_ok=True)
