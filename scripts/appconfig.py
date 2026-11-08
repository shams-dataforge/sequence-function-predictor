"""
appconfig.py  -  settings that differ between your laptop and the institute server.

Change a value here, or (on a server) set the environment variable named in brackets,
so no code has to be edited.
"""
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def _env(name: str, default: str) -> str:
    return os.environ.get(name, default)


# ---- name shown at the top of the page
APP_NAME = _env("SEQ2FUNC_APP_NAME", "Sequence Function Predictor")   # [SEQ2FUNC_APP_NAME]
APP_VERSION = "1.0.0"
INSTITUTE = _env("SEQ2FUNC_INSTITUTE", "")        # e.g. "National Institute of ..."  [SEQ2FUNC_INSTITUTE]
CONTACT_EMAIL = _env("SEQ2FUNC_CONTACT", "")      # shown in the footer / Help      [SEQ2FUNC_CONTACT]

# ---- folders
DB_DIR = Path(_env("SEQ2FUNC_DB_DIR", str(Path.home() / "eggnog_data")))             # eggNOG database
RUNS_DIR = Path(_env("SEQ2FUNC_RUNS_DIR", str(PROJECT_ROOT / "data" / "runs")))      # saved results
UPLOAD_DIR = Path(_env("SEQ2FUNC_UPLOAD_DIR", str(PROJECT_ROOT / "data" / "uploads")))
CACHE_DIR = Path(_env("SEQ2FUNC_CACHE_DIR", str(PROJECT_ROOT / "data" / "cache")))   # KEGG name lists

# ---- protection of a shared server
MAX_SEQUENCES = int(_env("SEQ2FUNC_MAX_SEQUENCES", "20000"))   # most sequences allowed in one run
# "1" = everybody sees all saved runs (fine on your own laptop).
# "0" = each visitor only sees runs started in their own browser session, and opens older
#       ones with the private Run ID.  Use "0" when the app is hosted for several users.
SHARED_RUNS = _env("SEQ2FUNC_SHARED_RUNS", "1") == "1"

# ---- automatic deletion of saved runs and uploaded files
# Runs and uploads are deleted from the server this many days after they were made. 0 = never delete.
RETENTION_DAYS = int(_env("SEQ2FUNC_RETENTION_DAYS", "2"))   # [SEQ2FUNC_RETENTION_DAYS]
