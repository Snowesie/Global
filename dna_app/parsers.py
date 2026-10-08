"""Readers for consumer raw-data exports and VCF files.

Every reader yields ``Variant`` tuples so the rest of the pipeline does not
care where the data came from. Supported inputs:

* 23andMe / LivingDNA  (rsid, chromosome, position, genotype — tab separated)
* AncestryDNA          (rsid, chromosome, position, allele1, allele2)
* MyHeritage / FTDNA   (CSV: RSID, CHROMOSOME, POSITION, RESULT)
* VCF (.vcf / .vcf.gz) — first sample column is used
* Any of the above wrapped in .gz or .zip
"""

from __future__ import annotations

import gzip
import io
import re
import zipfile
from dataclasses import dataclass, field
from typing import IO, Iterator, NamedTuple

VALID_CHROMS = {str(i) for i in range(1, 23)} | {"X", "Y", "MT"}

# AncestryDNA encodes the sex chromosomes and mitochondria numerically.
_NUMERIC_CHROM = {"23": "X", "24": "Y", "25": "X", "26": "MT"}
_BASES = set("ACGT")


class Variant(NamedTuple):
    rsid: str | None
    chrom: str
    pos: int
    ref: str | None  # None for genotype-array data (reference allele unknown)
    alleles: tuple[str, ...]  # the two called alleles, e.g. ("A", "G")


@dataclass
class ParseStats:
    format: str = "unknown"
    assembly: str | None = None
    sample: str | None = None
    lines: int = 0
    variants: int = 0
    no_calls: int = 0
    skipped: int = 0
    notes: list[str] = field(default_factory=list)


class UnsupportedFormat(ValueError):
    pass


def normalize_chrom(raw: str) -> str | None:
    c = raw.strip().strip('"').upper()
    if c.startswith("CHR"):
        c = c[3:]
    if c in ("M", "MITO"):
        c = "MT"
    if c == "XY":
        c = "X"
    c = _NUMERIC_CHROM.get(c, c)
    return c if c in VALID_CHROMS else None


def open_text(path: str) -> IO[str]:
    """Open a possibly gzip- or zip-compressed file as text."""
    with open(path, "rb") as fh:
        magic = fh.read(4)
    if magic[:2] == b"\x1f\x8b":
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    if magic == b"PK\x03\x04":
        zf = zipfile.ZipFile(path)
        members = [
            m for m in zf.infolist()
            if not m.is_dir() and not m.filename.startswith("__MACOSX") and m.file_size > 0
        ]
        if not members:
            raise UnsupportedFormat("The zip archive is empty.")
        member = max(members, key=lambda m: m.file_size)
        raw = zf.open(member)
        if member.filename.endswith(".gz"):
            raw = gzip.GzipFile(fileobj=raw)
        return io.TextIOWrapper(raw, encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def _sniff(path: str) -> tuple[str, list[str]]:
    """Return (format, header_lines) by peeking at the start of the file."""
    header: list[str] = []
    with open_text(path) as fh:
        for i, line in enumerate(fh):
            if i > 5000:
                break
            s = line.strip()
            if not s:
                continue
            if s.startswith("##fileformat=VCF"):
                return "vcf", header
            if s.startswith("#"):
                header.append(s)
                if s.startswith("#CHROM"):
                    return "vcf", header
                continue
            low = s.lower().replace('"', "")
            if low.startswith("rsid\tchromosome\tposition\tallele1"):
                return "ancestrydna", header
            if low.startswith("rsid,chromosome,position,result"):
                return "myheritage", header
            if low.startswith("rsid\tchromosome\tposition\tgenotype") or low.startswith("# rsid"):
                return "23andme", header
            # Headerless data line: decide on shape.
            if "," in s and s.count(",") >= 3:
                return "myheritage", header
            parts = s.split("\t")
            if len(parts) == 4:
                return "23andme", header
            if len(parts) == 5:
                return "ancestrydna", header
            if len(parts) >= 8:
                return "vcf", header
            break
    raise UnsupportedFormat(
        "Could not recognise this file. Supported formats: 23andMe, AncestryDNA, "
        "MyHeritage, FamilyTreeDNA, LivingDNA raw data, and VCF."
    )


def _guess_array_assembly(header: list[str]) -> str:
    text = " ".join(header).lower()
    if "grch38" in text or "build 38" in text or "hg38" in text:
        return "GRCh38"
    # 23andMe, AncestryDNA, MyHeritage and FTDNA all export build 37.
    return "GRCh37"


_VCF_CONTIG_LEN = re.compile(r"##contig=<ID=(?:chr)?1,.*?length=(\d+)", re.I)


def _guess_vcf_assembly(meta: list[str]) -> str | None:
    text = "\n".join(meta)
    m = _VCF_CONTIG_LEN.search(text)
    if m:
        length = int(m.group(1))
        if length == 249250621:
            return "GRCh37"
        if length == 248956422:
            return "GRCh38"
    low = text.lower()
    if any(k in low for k in ("grch38", "hg38", "b38")):
        return "GRCh38"
    if any(k in low for k in ("grch37", "hg19", "b37", "hs37d5", "human_g1k_v37")):
        return "GRCh37"
    return None


def parse_file(path: str, stats: ParseStats, assembly_override: str | None = None) -> Iterator[Variant]:
    fmt, header = _sniff(path)
    stats.format = fmt
    if fmt == "vcf":
        yield from _parse_vcf(path, stats, assembly_override)
    else:
        stats.assembly = assembly_override or _guess_array_assembly(header)
        yield from _parse_array(path, fmt, stats)


def _split_genotype(gt: str) -> tuple[str, ...] | None:
    gt = gt.strip().strip('"').upper()
    if len(gt) == 2 and gt[0] in _BASES and gt[1] in _BASES:
        return (gt[0], gt[1])
    if len(gt) == 1 and gt in _BASES:  # hemizygous X/Y/MT calls
        return (gt, gt)
    return None


def _parse_array(path: str, fmt: str, stats: ParseStats) -> Iterator[Variant]:
    sep = "," if fmt == "myheritage" else "\t"
    with open_text(path) as fh:
        for line in fh:
            stats.lines += 1
            if not line.strip() or line.startswith("#"):
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\r\n").split(sep)]
            if parts[0].lower() == "rsid":
                continue
            if fmt == "ancestrydna":
                if len(parts) < 5:
                    stats.skipped += 1
                    continue
                gt = parts[3] + parts[4]
            else:
                if len(parts) < 4:
                    stats.skipped += 1
                    continue
                gt = parts[3]
            chrom = normalize_chrom(parts[1])
            try:
                pos = int(parts[2])
            except ValueError:
                stats.skipped += 1
                continue
            if chrom is None:
                stats.skipped += 1
                continue
            alleles = _split_genotype(gt)
            if alleles is None:
                stats.no_calls += 1
                continue
            rsid = parts[0] if parts[0].lower().startswith("rs") else None
            stats.variants += 1
            yield Variant(rsid.lower() if rsid else None, chrom, pos, None, alleles)


def _parse_vcf(path: str, stats: ParseStats, assembly_override: str | None) -> Iterator[Variant]:
    meta: list[str] = []
    gt_idx: int | None = None
    sample_col = 9
    with open_text(path) as fh:
        for line in fh:
            stats.lines += 1
            if line.startswith("##"):
                meta.append(line.rstrip())
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\r\n").split("\t")
                if len(cols) > 9:
                    stats.sample = cols[9]
                    if len(cols) > 10:
                        stats.notes.append(
                            f"VCF has {len(cols) - 9} samples; analysing the first ({cols[9]})."
                        )
                else:
                    stats.notes.append("VCF has no sample columns; treating every site as heterozygous.")
                stats.assembly = assembly_override or _guess_vcf_assembly(meta) or "GRCh37"
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 5:
                stats.skipped += 1
                continue
            chrom = normalize_chrom(cols[0])
            if chrom is None:
                stats.skipped += 1
                continue
            try:
                pos = int(cols[1])
            except ValueError:
                stats.skipped += 1
                continue
            ref = cols[3].upper()
            alts = cols[4].upper().split(",")
            all_alleles = [ref] + alts
            if len(cols) > sample_col:
                fmt_keys = cols[8].split(":")
                if gt_idx is None or (gt_idx < len(fmt_keys) and fmt_keys[gt_idx] != "GT"):
                    gt_idx = fmt_keys.index("GT") if "GT" in fmt_keys else None
                if gt_idx is None:
                    stats.skipped += 1
                    continue
                fields = cols[sample_col].split(":")
                gt = fields[gt_idx] if gt_idx < len(fields) else "."
                idxs = re.split(r"[/|]", gt)
                if any(i in (".", "") for i in idxs):
                    stats.no_calls += 1
                    continue
                try:
                    called = tuple(all_alleles[int(i)] for i in idxs)
                except (ValueError, IndexError):
                    stats.skipped += 1
                    continue
                if len(called) == 1:
                    called = (called[0], called[0])
            else:
                called = (ref, alts[0])
            if cols[4] in (".", "<NON_REF>", "<*>"):
                # Reference-only / gVCF block: no alternate allele to annotate.
                stats.skipped += 1
                continue
            rsid = None
            for ident in cols[2].split(";"):
                if ident.lower().startswith("rs"):
                    rsid = ident.lower()
                    break
            stats.variants += 1
            yield Variant(rsid, chrom, pos, ref, called)
    if stats.assembly is None:
        stats.assembly = assembly_override or "GRCh37"
