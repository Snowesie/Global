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
docker run -p 8000:8000 -v dna-data:/data -e DNA_AUTO_BUILD=1 -e DNA_APP_PASSWORD=choose-one genome-reader
```

With `DNA_AUTO_BUILD=1`, the container downloads and builds any missing database in the
background on first boot. The page shows progress and holds uploads until the build finishes.

## Deploy to Render (public link)

[![Deploy to Render](https://render.com/images/deploy-to-render-button.svg)](https://render.com/deploy?repo=https://github.com/Snowesie/Global)

The repo includes a [`render.yaml`](render.yaml) Blueprint: a Docker web service with a 10 GB
persistent disk for the databases and uploads.

1. Merge this branch into the default branch. The deploy button reads `render.yaml` from there.
   You can also choose the branch yourself under **New → Blueprint** in the Render dashboard.
2. Click the button (or **New → Blueprint**) and connect the GitHub repo.
3. When Render asks for `DNA_APP_PASSWORD`, pick a password. Anyone opening the site gets a
   browser login prompt; any username works with that password.
4. Create the service. You get a link like `https://genome-reader-xxxx.onrender.com`.
5. On first boot the app downloads ClinVar and 1000 Genomes. That takes roughly 30–60 minutes,
   and the page shows progress while it runs. The databases are stored on the disk, so later
   restarts and deploys skip this step.

**Cost:** persistent disks need a paid instance, so this is roughly the Starter instance plus
10 GB of disk per month (see Render's pricing). Starter (512 MB RAM) handles chip files and
typical VCFs. Switch `plan` to `standard` if whole-genome VCFs run out of memory.

**Updating ClinVar:** delete `/data/clinvar.sqlite` and `/data/1000g.sqlite` from the Render
shell and restart the service. The 1000 Genomes database is rebuilt too, because it keeps the
rare sites that appear in ClinVar.

The browser uploads in 32 MB chunks and retries failed chunks, so 1 GB files work behind
proxies that cap request size.

## How it works

1. **Upload.** The browser sends the file in 32 MB chunks, each streamed straight to disk, and
   retries any chunk that fails. Files over 1 GiB are refused.
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
| `DNA_APP_PASSWORD` | unset | Require this password for every page and API call (HTTP Basic auth; any username) |
| `DNA_AUTO_BUILD` | unset | With `start.sh`/Docker: set to `1` to build missing databases on boot |

## API

| Method | Path | |
| --- | --- | --- |
| `GET` | `/api/health` | Health check (never needs the password) |
| `GET` | `/api/status` | Upload limit, database versions, build progress |
| `POST` | `/api/upload?filename=…&assembly=auto\|GRCh37\|GRCh38` | Whole file as the raw body → `{job_id}` (202). Handy for `curl --data-binary @file` |
| `POST` | `/api/uploads?size=…&filename=…&assembly=…` | Start a chunked upload → `{upload_id, chunk_bytes}` |
| `PUT` | `/api/uploads/{id}?offset=…` | Send the next chunk (raw body) → `{received}` |
| `GET` | `/api/uploads/{id}` | Bytes received so far, for resuming |
| `POST` | `/api/uploads/{id}/complete` | Queue the analysis → `{job_id}` (202) |
| `GET` | `/api/jobs/{id}` | `queued` / `running` / `done` / `error`, with progress |
| `GET` | `/api/jobs/{id}/result` | Full JSON report |
| `DELETE` | `/api/jobs/{id}` | Delete the report |

Job IDs are random 144-bit tokens. With `DNA_APP_PASSWORD` set, everyone shares the one
password, and anyone who has both the password and a job ID can read that report. Always set a
password on a public deployment.

## Development

```bash
pip install -r requirements-dev.txt
pytest
```

The tests build small synthetic ClinVar and 1000 Genomes files in the real upstream formats,
simulate genomes with known ancestry mixtures, and check the full pipeline: parsing, ClinVar
matching, strand flips, multi-allelic VCF sites, frequencies, ancestry recovery and the upload
size limit.
