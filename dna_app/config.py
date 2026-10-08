import os
from pathlib import Path

MAX_UPLOAD_BYTES = int(os.environ.get("DNA_MAX_UPLOAD_BYTES", 1024 ** 3))  # 1 GiB

DATA_DIR = Path(os.environ.get("DNA_DATA_DIR", Path(__file__).resolve().parent.parent / "data"))
UPLOAD_DIR = DATA_DIR / "uploads"
RESULTS_DIR = DATA_DIR / "results"
CLINVAR_DB = Path(os.environ.get("DNA_CLINVAR_DB", DATA_DIR / "clinvar.sqlite"))
KG_DB = Path(os.environ.get("DNA_1000G_DB", DATA_DIR / "1000g.sqlite"))

# When set, every page and API call requires this password (HTTP Basic auth).
APP_PASSWORD = os.environ.get("DNA_APP_PASSWORD") or None

# Delete the raw upload once it has been analysed (results are kept).
DELETE_UPLOADS = os.environ.get("DNA_KEEP_UPLOADS", "0") != "1"

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
KG_SITES_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/"
    "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5c.20130502.sites.vcf.gz"
)

SUPERPOPS = ("AFR", "AMR", "EAS", "EUR", "SAS")
SUPERPOP_NAMES = {
    "AFR": "African",
    "AMR": "Admixed American",
    "EAS": "East Asian",
    "EUR": "European",
    "SAS": "South Asian",
}
