# Genome Reader

A self-hosted web app that reads your raw DNA file (up to 1 GB) and interprets it with two public
resources:

- **[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/)** (NCBI): which of your variants have a
  clinical classification (pathogenic, risk factor, drug response, and so on), your zygosity,
  the review-status stars and the linked conditions.
- **[1000 Genomes Project Phase 3](https://www.internationalgenome.org/)**: how common each finding
  is in the five superpopulations (African, Admixed American, East Asian, European, South Asian),
  plus a continental ancestry estimate with bootstrap confidence intervals.

Everything runs on your own machine. Uploads are deleted as soon as the analysis finishes.

> **Not medical advice.** Consumer genotyping chips have high false-positive rates for rare
> pathogenic variants. Confirm anything important with clinical testing and a genetic counsellor.

## Supported files

| Source | Format |
| --- | --- |
| 23andMe, LivingDNA | `rsid  chromosome  position  genotype` (tab) |
| AncestryDNA | `rsid  chromosome  position  allele1  allele2` (tab) |
| MyHeritage, FamilyTreeDNA | CSV `RSID,CHROMOSOME,POSITION,RESULT` |
| Sequencing (WGS/WES, Nebula, Dante, …) | VCF / `.vcf.gz` (first sample is used) |

Any of these may be uploaded as-is or inside a `.gz` or `.zip`. The genome build (GRCh37/GRCh38)
is detected from the file header, or you can set it in the UI.

## Quick start

```bash
pip install -r requirements.txt

# 1. Build the reference databases (needs internet; one-off)
python -m dna_app.build_db clinvar   # ~400 MB download, ~1 GB SQLite, a few minutes
python -m dna_app.build_db 1000g     # ~1.4 GB download, ~1-2 GB SQLite, 20-40 minutes

# 2. Run the app
uvicorn dna_app.server:app --host 0.0.0.0 --port 8000
# open http://localhost:8000
```

Build ClinVar before 1000 Genomes: the 1000 Genomes build drops sites rarer than 1% in every
superpopulation (`--min-af`), but always keeps sites that appear in ClinVar so rare pathogenic
variants still get frequencies. Both builders accept `--source` to use a file you've already
downloaded:

```bash
python -m dna_app.build_db clinvar --source ~/Downloads/variant_summary.txt.gz
python -m dna_app.build_db 1000g --source ~/Downloads/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz
```

Re-run `build_db clinvar` periodically; ClinVar is updated every month.

### Docker

```bash
docker build -t genome-reader .
docker run --rm -v dna-data:/data genome-reader python -m dna_app.build_db clinvar
docker run --rm -v dna-data:/data genome-reader python -m dna_app.build_db 1000g
docker run -p 8000:8000 -v dna-data:/data genome-reader
```

If you put it behind a reverse proxy, raise the body-size limit (for example nginx
`client_max_body_size 1g;`).

## How it works

1. **Upload.** The browser sends the file as a raw request body. The server streams it to disk and
   stops at 1 GiB, whether or not the client sent a `Content-Length`.
2. **Parse.** The format is sniffed from the first lines. Genotypes go into a per-job SQLite table,
   so a whole-genome VCF never has to fit in memory. A 900 MB, 6-million-variant VCF takes under a
   minute.
3. **ClinVar.** Genotypes are joined to ClinVar's `variant_summary` on chromosome/position for the
   file's build, and on rsID too for chip data. A finding is reported only when you carry the
   classified allele:
   - VCF alleles are compared exactly.
   - Chip genotypes have no reference allele, so they are checked on both strands, except A/T and
     C/G SNPs, where a strand flip can't be detected.
   - Indels are skipped for chip data because chips can't call them reliably.
4. **1000 Genomes.** Each finding gets its allele frequency per superpopulation. Ancestry comes
   from a supervised admixture model on up to 200,000 ancestry-informative autosomal SNPs. The
   superpopulation frequencies are held fixed, and your mixture proportions are fitted by EM (the
   same update ADMIXTURE uses in supervised mode). For whole-genome VCFs, which list only
   variant sites, a common site missing from the file is counted as homozygous reference.

The 1000 Genomes Phase 3 sites file is on GRCh37. GRCh37 files are matched on position. GRCh38
files are matched on rsID, so a GRCh38 VCF needs rsIDs in its `ID` column to get frequencies and
ancestry.

## Configuration

| Variable | Default | |
| --- | --- | --- |
| `DNA_DATA_DIR` | `./data` | Where uploads, results and databases live |
| `DNA_CLINVAR_DB` | `$DNA_DATA_DIR/clinvar.sqlite` | |
| `DNA_1000G_DB` | `$DNA_DATA_DIR/1000g.sqlite` | |
| `DNA_MAX_UPLOAD_BYTES` | `1073741824` | Upload limit (1 GiB) |
| `DNA_KEEP_UPLOADS` | `0` | Set to `1` to keep raw uploads after analysis |

## API

| Method | Path | |
| --- | --- | --- |
| `GET` | `/api/status` | Upload limit and database versions |
| `POST` | `/api/upload?filename=…&assembly=auto\|GRCh37\|GRCh38` | Raw file body → `{job_id}` (202) |
| `GET` | `/api/jobs/{id}` | `queued` / `running` / `done` / `error`, with progress |
| `GET` | `/api/jobs/{id}/result` | Full JSON report |
| `DELETE` | `/api/jobs/{id}` | Delete the report |

Job IDs are random 144-bit tokens, and anyone who has one can read that report. There is no user
login, so don't expose the server publicly without putting authentication in front of it.

## Development

```bash
pip install -r requirements-dev.txt
pytest
```

The tests build small synthetic ClinVar and 1000 Genomes files in the real upstream formats,
simulate genomes with known ancestry mixtures, and check the full pipeline: parsing, ClinVar
matching, strand flips, multi-allelic VCF sites, frequencies, ancestry recovery and the upload
size limit.
