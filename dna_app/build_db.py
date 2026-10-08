"""Build the local reference databases the analyser reads.

    python -m dna_app.build_db clinvar            # downloads variant_summary.txt.gz
    python -m dna_app.build_db 1000g              # downloads the phase 3 sites VCF (~1.4 GB)
    python -m dna_app.build_db clinvar --source /path/to/variant_summary.txt.gz

Both commands accept ``--source`` (local path or URL) and ``--out``.
Build ClinVar first: the 1000 Genomes build keeps every ClinVar site even
when it is rarer than ``--min-af``.
"""

from __future__ import annotations

import argparse
import gzip
import io
import os
import sqlite3
import sys
import time
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from . import config
from .parsers import normalize_chrom

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
}


def _open_source(source: str, cache_dir: Path) -> io.TextIOBase:
    if source.startswith(("http://", "https://", "ftp://")):
        cache_dir.mkdir(parents=True, exist_ok=True)
        local = cache_dir / source.rsplit("/", 1)[-1]
        if not local.exists():
            print(f"Downloading {source} -> {local}", file=sys.stderr)
            tmp = local.with_suffix(local.suffix + ".part")
            with urllib.request.urlopen(source) as resp, open(tmp, "wb") as out:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                last = 0.0
                while chunk := resp.read(1 << 20):
                    out.write(chunk)
                    done += len(chunk)
                    if time.time() - last > 2:
                        last = time.time()
                        pct = f" ({done * 100 // total}%)" if total else ""
                        print(f"  {done / 1e6:,.0f} MB{pct}", file=sys.stderr)
            tmp.rename(local)
        source = str(local)
    if source.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(source, "rb"), encoding="utf-8", errors="replace")
    return open(source, "r", encoding="utf-8", errors="replace")


def _write_meta(con: sqlite3.Connection, **items: str) -> None:
    con.execute("CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT)")
    items.setdefault("built_at", datetime.now(timezone.utc).isoformat(timespec="seconds"))
    con.executemany("INSERT OR REPLACE INTO meta VALUES (?, ?)", items.items())


def _fresh_db(out: Path) -> sqlite3.Connection:
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(out.suffix + ".building")
    if tmp.exists():
        tmp.unlink()
    con = sqlite3.connect(tmp)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    return con


def _finish_db(con: sqlite3.Connection, out: Path) -> None:
    con.commit()
    con.execute("ANALYZE")
    con.close()
    os.replace(out.with_suffix(out.suffix + ".building"), out)


# --------------------------------------------------------------------------- ClinVar


def build_clinvar(source: str, out: Path) -> int:
    con = _fresh_db(out)
    con.execute(
        """CREATE TABLE clinvar (
            variation_id INTEGER, allele_id INTEGER, assembly TEXT,
            chrom TEXT, pos INTEGER, ref TEXT, alt TEXT, rsid TEXT,
            gene TEXT, name TEXT, type TEXT, significance TEXT,
            review_status TEXT, stars INTEGER, phenotypes TEXT,
            origin TEXT, last_evaluated TEXT
        )"""
    )
    n = 0
    with _open_source(source, out.parent / "downloads") as fh:
        header = fh.readline().lstrip("#").rstrip("\r\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        required = ["AlleleID", "Assembly", "Chromosome", "PositionVCF",
                    "ReferenceAlleleVCF", "AlternateAlleleVCF", "ClinicalSignificance"]
        missing = [c for c in required if c not in col]
        if missing:
            raise SystemExit(f"Not a ClinVar variant_summary file (missing columns: {missing})")

        def get(row: list[str], name: str) -> str:
            i = col.get(name)
            return row[i] if i is not None and i < len(row) else ""

        batch = []
        for line in fh:
            row = line.rstrip("\r\n").split("\t")
            assembly = get(row, "Assembly")
            if assembly not in ("GRCh37", "GRCh38"):
                continue
            chrom = normalize_chrom(get(row, "Chromosome"))
            ref = get(row, "ReferenceAlleleVCF").upper()
            alt = get(row, "AlternateAlleleVCF").upper()
            try:
                pos = int(get(row, "PositionVCF"))
            except ValueError:
                pos = -1
            if chrom is None or pos <= 0 or ref in ("", "NA") or alt in ("", "NA"):
                continue
            rs = get(row, "RS# (dbSNP)")
            rsid = f"rs{rs}" if rs and rs != "-1" else None
            review = get(row, "ReviewStatus")
            batch.append((
                int(get(row, "VariationID") or 0), int(get(row, "AlleleID") or 0), assembly,
                chrom, pos, ref, alt, rsid,
                get(row, "GeneSymbol"), get(row, "Name"), get(row, "Type"),
                get(row, "ClinicalSignificance"), review, STARS.get(review, 0),
                get(row, "PhenotypeList"), get(row, "OriginSimple"), get(row, "LastEvaluated"),
            ))
            if len(batch) >= 50_000:
                con.executemany(f"INSERT INTO clinvar VALUES ({','.join('?' * 17)})", batch)
                n += len(batch)
                batch.clear()
                print(f"  {n:,} ClinVar records", file=sys.stderr)
        con.executemany(f"INSERT INTO clinvar VALUES ({','.join('?' * 17)})", batch)
        n += len(batch)
    print("Indexing…", file=sys.stderr)
    con.execute("CREATE INDEX clinvar_pos ON clinvar (assembly, chrom, pos)")
    con.execute("CREATE INDEX clinvar_rsid ON clinvar (rsid)")
    _write_meta(con, source=source, records=str(n))
    _finish_db(con, out)
    print(f"Wrote {n:,} ClinVar records to {out}", file=sys.stderr)
    return n


# --------------------------------------------------------------------------- 1000 Genomes


def _clinvar_keys(clinvar_db: Path) -> tuple[set[str], set[tuple[str, int]]]:
    if not clinvar_db.exists():
        return set(), set()
    con = sqlite3.connect(clinvar_db)
    rsids = {r for (r,) in con.execute("SELECT DISTINCT rsid FROM clinvar WHERE rsid IS NOT NULL")}
    positions = set(con.execute("SELECT DISTINCT chrom, pos FROM clinvar WHERE assembly='GRCh37'"))
    con.close()
    return rsids, positions


def build_1000g(source: str, out: Path, min_af: float, clinvar_db: Path) -> int:
    """Load per-superpopulation allele frequencies from the phase 3 sites VCF.

    Sites are kept if any superpopulation frequency is at least ``min_af`` or
    the site is in ClinVar (so rare pathogenic variants still get frequencies).
    """
    keep_rsids, keep_pos = _clinvar_keys(clinvar_db)
    if not keep_rsids:
        print("Note: ClinVar database not found; rare ClinVar sites may be filtered out.", file=sys.stderr)
    con = _fresh_db(out)
    con.execute(
        """CREATE TABLE af (
            chrom TEXT, pos INTEGER, rsid TEXT, ref TEXT, alt TEXT,
            af REAL, afr REAL, amr REAL, eas REAL, eur REAL, sas REAL
        )"""
    )
    keys = ("AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF")
    n = 0
    batch = []
    with _open_source(source, out.parent / "downloads") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t", 8)
            chrom = normalize_chrom(cols[0])
            if chrom is None:
                continue
            pos = int(cols[1])
            rsid = next((i.lower() for i in cols[2].split(";") if i.lower().startswith("rs")), None)
            ref = cols[3].upper()
            alts = cols[4].upper().split(",")
            info: dict[str, list[str]] = {}
            for kv in cols[7].split(";"):
                k, _, v = kv.partition("=")
                if k in keys:
                    info[k] = v.split(",")
            if "AF" not in info:
                continue
            in_clinvar = (rsid in keep_rsids) or ((chrom, pos) in keep_pos)
            for i, alt in enumerate(alts):
                if alt.startswith("<"):
                    continue
                freqs = []
                for k in keys:
                    vals = info.get(k)
                    try:
                        freqs.append(float(vals[i]) if vals else None)
                    except (ValueError, IndexError):
                        freqs.append(None)
                pop_max = max((f for f in freqs[1:] if f is not None), default=0.0)
                pop_min = min((f for f in freqs[1:] if f is not None), default=0.0)
                common = pop_max >= min_af and pop_min <= 1 - min_af
                if not (common or in_clinvar):
                    continue
                batch.append((chrom, pos, rsid, ref, alt, *freqs))
            if len(batch) >= 100_000:
                con.executemany("INSERT INTO af VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
                n += len(batch)
                batch.clear()
                print(f"  {n:,} sites (at {chrom}:{pos})", file=sys.stderr)
    con.executemany("INSERT INTO af VALUES (?,?,?,?,?,?,?,?,?,?,?)", batch)
    n += len(batch)
    print("Indexing…", file=sys.stderr)
    con.execute("CREATE INDEX af_rsid ON af (rsid)")
    con.execute("CREATE INDEX af_pos ON af (chrom, pos)")
    _write_meta(con, source=source, records=str(n), assembly="GRCh37",
                release="1000 Genomes Phase 3 (20130502)", min_af=str(min_af))
    _finish_db(con, out)
    print(f"Wrote {n:,} 1000 Genomes sites to {out}", file=sys.stderr)
    return n


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("clinvar", help="Build the ClinVar database")
    c.add_argument("--source", default=config.CLINVAR_URL)
    c.add_argument("--out", type=Path, default=config.CLINVAR_DB)
    g = sub.add_parser("1000g", help="Build the 1000 Genomes allele-frequency database")
    g.add_argument("--source", default=config.KG_SITES_URL)
    g.add_argument("--out", type=Path, default=config.KG_DB)
    g.add_argument("--min-af", type=float, default=0.01,
                   help="Drop sites rarer than this in every superpopulation (ClinVar sites are always kept)")
    g.add_argument("--clinvar-db", type=Path, default=config.CLINVAR_DB)
    a = p.parse_args(argv)
    if a.cmd == "clinvar":
        build_clinvar(a.source, a.out)
    else:
        build_1000g(a.source, a.out, a.min_af, a.clinvar_db)


if __name__ == "__main__":
    main()
