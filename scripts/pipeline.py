"""
pipeline.py  -  the analysis ENGINE (no web-page code in here).

What it does, in order
  1. checks and cleans your FASTA file
  2. finds genes / ORFs when the input is DNA (Prodigal, or a longest-ORF finder)
  3. runs eggNOG-mapper on the proteins
  4. saves everything inside ONE folder per run:  data/runs/<date>_<name>/

Turning those raw files into readable tables happens in results.py, so changing a
threshold in the app never means re-running the (slow) search.
The same functions will later be reused by the command-line tool.
"""
from __future__ import annotations

import json
import os
import re
import secrets
import shutil
import signal
import statistics
import subprocess
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Optional

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

Log = Callable[[str], None]

import appconfig  # noqa: E402

PROJECT_ROOT = appconfig.PROJECT_ROOT
RUNS_DIR = appconfig.RUNS_DIR
DEFAULT_DB_DIR = appconfig.DB_DIR

# --------------------------------------------------------------------------- choices
INPUT_TRANSCRIPT = "transcript"
INPUT_PROKARYOTE = "prokaryote"
INPUT_PROTEIN = "protein"
INPUT_TYPE_LABELS = {
    INPUT_TRANSCRIPT: "DNA - transcripts, mRNA or coding sequences (any organism)",
    INPUT_PROKARYOTE: "DNA - bacterial, archaeal or viral genomes / contigs / genes",
    INPUT_PROTEIN: "Protein sequences",
}

# label shown in the app -> value given to eggNOG-mapper's --tax_scope (None = automatic)
TAX_SCOPES = {
    "Auto - let eggNOG choose per sequence (recommended)": None,
    "Bacteria": "bacteria",
    "Archaea": "archaea",
    "Eukaryota (all)": "eukaryota",
    "Plants (Viridiplantae)": "33090,2759,1",
    "Animals (Metazoa)": "33208,2759,1",
    "Fungi": "4751,2759,1",
}
GO_EVIDENCE = {
    "Non-electronic (eggNOG default)": "non-electronic",
    "Experimental only (strictest)": "experimental",
    "All, including electronic IEA (largest coverage)": "all",
}
SENSMODES = ["fast", "default", "mid-sensitive", "sensitive",
             "more-sensitive", "very-sensitive", "ultra-sensitive"]
TARGET_ORTHOLOGS = ["all", "one2one", "many2one", "one2many", "many2many"]

MIN_DNA_NT = 30
MIN_PROTEIN_AA = 10
_IUPAC_DNA = set("ACGTURYKMSWBDHVN")
_IUPAC_TO_N = str.maketrans("RYKMSWBDHV", "N" * 10)
_AMINO = set("ACDEFGHIKLMNPQRSTVWYX")
_PRODIGAL_HDR = re.compile(r"^(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#")


class PipelineError(RuntimeError):
    """A problem the user can fix (missing database, wrong file type, ...)."""


@dataclass
class RunSettings:
    input_type: str = INPUT_PROKARYOTE
    cpu: int = 4
    sensmode: str = "fast"
    evalue: float = 0.001            # eggNOG-mapper default
    min_score: float = 60.0          # eggNOG-mapper default
    min_pident: float = 0.0          # 0 = no filter (eggNOG-mapper default)
    min_query_cover: float = 0.0     # 0 = no filter
    min_subject_cover: float = 0.0   # 0 = no filter
    tax_scope: Optional[str] = None  # None = automatic (eggNOG-mapper default)
    target_orthologs: str = "all"    # eggNOG-mapper default
    go_evidence: str = "non-electronic"  # eggNOG-mapper default
    min_orf_aa: int = 100            # transcripts only
    genetic_code: int = 1            # transcripts only (1 = standard)
    orf_both_strands: bool = True    # transcripts only; False = forward strand only

    def to_dict(self) -> dict:
        return asdict(self)


# --------------------------------------------------------------------------- helpers
def slugify(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(name)).strip("_")[:50] or "run"


def run_command(cmd: list, log: Log, cwd: Optional[Path] = None) -> None:
    """Run an external program and stream its output to `log` line by line."""
    log("$ " + " ".join(str(c) for c in cmd))
    try:
        proc = subprocess.Popen([str(c) for c in cmd], stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1,
                                encoding="utf-8", errors="replace", cwd=cwd,
                                start_new_session=True)   # own process group, so ALL its children can be stopped
    except FileNotFoundError:
        raise PipelineError(f"Program not found: {cmd[0]}. Did you run "
                            "`conda activate bioannotate` first?")
    tail: list[str] = []
    assert proc.stdout is not None
    try:
        for line in proc.stdout:
            line = line.rstrip()
            if line:
                log(line)
                tail.append(line)
                tail = tail[-15:]
        code = proc.wait()
    except BaseException:          # e.g. the user closed the page: stop the program AND its children
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except Exception:
            proc.kill()
        proc.wait()
        raise
    if code != 0:
        raise PipelineError(f"{cmd[0]} stopped with error code {code}. "
                            "Last output:\n" + "\n".join(tail))


def check_environment(input_type: str, data_dir) -> list[str]:
    """Return a list of problems that would stop a run (empty list = all good)."""
    problems = []
    if shutil.which("emapper.py") is None:
        problems.append("eggNOG-mapper (`emapper.py`) was not found. Start the app from "
                        "the conda environment:  `conda activate bioannotate`.")
    if input_type == INPUT_PROKARYOTE and shutil.which("prodigal") is None:
        problems.append("`prodigal` was not found. Start the app from the conda "
                        "environment:  `conda activate bioannotate`.")
    d = Path(data_dir).expanduser()
    for f in ("eggnog.db", "eggnog_proteins.dmnd"):
        if not (d / f).exists():
            problems.append(f"eggNOG database file `{f}` was not found in `{d}`.")
    return problems



def _parse_fasta(path):
    """Read a FASTA file; a clear message instead of a cryptic error for non-FASTA files."""
    try:
        yield from SeqIO.parse(str(path), "fasta")
    except ValueError as exc:
        raise PipelineError("This file does not look like a FASTA file (the first line must "
                            "start with '>').") from exc

def _iter_for_inspection(path):
    yield from _parse_fasta(path)


# --------------------------------------------------------------------------- file check
def inspect_fasta(path) -> dict:
    """Quick look at a FASTA file: how many sequences, DNA or protein, lengths, warnings."""
    n, lengths, ids, dups = 0, [], set(), 0
    dna_like = sampled = 0
    try:
        records = list(_iter_for_inspection(path))
    except PipelineError as exc:
        return {"n_sequences": 0, "kind": "empty", "min_len": 0, "median_len": 0, "max_len": 0,
                "duplicate_ids": 0, "warnings": [str(exc)]}
    for rec in records:
        n += 1
        lengths.append(len(rec.seq))
        if rec.id in ids:
            dups += 1
        ids.add(rec.id)
        if sampled < 2_000_000:                       # sample, so huge files stay fast
            s = str(rec.seq).upper().replace("-", "").replace("*", "")
            sampled += len(s)
            dna_like += sum(s.count(c) for c in "ACGTUN")
    warnings = []
    if n == 0:
        kind = "empty"
        warnings.append("No sequences found. A FASTA file has header lines starting with '>'.")
    else:
        kind = "dna" if sampled and dna_like / sampled >= 0.9 else "protein"
    if dups:
        warnings.append(f"{dups} duplicate sequence IDs - they will be renamed automatically.")
    return {
        "n_sequences": n, "kind": kind,
        "min_len": min(lengths) if lengths else 0,
        "median_len": int(statistics.median(lengths)) if lengths else 0,
        "max_len": max(lengths) if lengths else 0,
        "duplicate_ids": dups, "warnings": warnings,
    }


# --------------------------------------------------------------------------- preparing sequences
def _unique_id(rid: str, seen: set, stats: dict) -> str:
    if rid in seen:
        k = 2
        while f"{rid}_dup{k}" in seen:
            k += 1
        rid = f"{rid}_dup{k}"
        stats["renamed_duplicates"] += 1
    seen.add(rid)
    return rid


def prepare_dna(input_fasta, out_fasta, log: Log) -> dict:
    """Clean a DNA FASTA: RNA->DNA, ambiguity codes -> N, drop junk, unique IDs."""
    stats = {"input": 0, "kept": 0, "too_short": 0, "invalid": 0,
             "renamed_duplicates": 0, "no_orf": 0}
    kept, seen = [], set()
    for rec in _parse_fasta(input_fasta):
        stats["input"] += 1
        seq = str(rec.seq).upper().replace("U", "T").replace("-", "").replace(" ", "")
        if len(seq) < MIN_DNA_NT:
            stats["too_short"] += 1
            continue
        if set(seq) - _IUPAC_DNA:
            stats["invalid"] += 1
            continue
        seq = seq.translate(_IUPAC_TO_N)
        kept.append(SeqRecord(Seq(seq), id=_unique_id(rec.id, seen, stats), description=""))
    stats["kept"] = len(kept)
    if not kept:
        raise PipelineError("No usable DNA sequences were found. Are these really DNA "
                            "sequences (letters A, C, G, T/U, N)? If they are protein "
                            "sequences, choose 'Protein sequences' as the sequence type.")
    SeqIO.write(kept, str(out_fasta), "fasta")
    log(f"      {stats['kept']} of {stats['input']} sequences kept "
        f"({stats['too_short']} too short, {stats['invalid']} with invalid characters)")
    return stats


def prepare_proteins(input_fasta, out_faa, log: Log):
    """Clean a protein FASTA. Returns (gene_map table, stats)."""
    stats = {"input": 0, "kept": 0, "too_short": 0, "invalid": 0,
             "renamed_duplicates": 0, "no_orf": 0}
    kept, rows, seen = [], [], set()
    for rec in _parse_fasta(input_fasta):
        stats["input"] += 1
        seq = str(rec.seq).upper().replace("*", "").replace("-", "").replace(" ", "")
        seq = re.sub(r"[UOBZJ]", "X", seq)
        if len(seq) < MIN_PROTEIN_AA:
            stats["too_short"] += 1
            continue
        if set(seq) - _AMINO:
            stats["invalid"] += 1
            continue
        rid = _unique_id(rec.id, seen, stats)
        kept.append(SeqRecord(Seq(seq), id=rid, description=""))
        rows.append({"protein_id": rid, "source_sequence": rid, "start": "", "end": "",
                     "strand": "", "protein_length_aa": len(seq), "gene_call": "Provided protein"})
    stats["kept"] = len(kept)
    if not kept:
        raise PipelineError("No usable protein sequences were found. If your file contains "
                            "DNA, choose one of the DNA sequence types instead.")
    SeqIO.write(kept, str(out_faa), "fasta")
    log(f"      {stats['kept']} of {stats['input']} proteins kept "
        f"({stats['too_short']} too short, {stats['invalid']} with invalid characters)")
    return pd.DataFrame(rows), stats


def run_prodigal(dna_fasta, run_dir: Path, out_faa, log: Log, stats: dict) -> pd.DataFrame:
    """Gene prediction for bacteria / archaea / viruses. Returns the gene_map table."""
    raw = run_dir / "prodigal_proteins.faa"
    run_command(["prodigal", "-i", dna_fasta, "-a", raw, "-d", run_dir / "genes.fna",
                 "-o", run_dir / "prodigal_coords.gbk", "-p", "meta", "-q"], log)
    kept, rows = [], []
    for rec in SeqIO.parse(str(raw), "fasta"):
        m = _PRODIGAL_HDR.match(rec.description)
        pid = rec.id
        source = pid.rsplit("_", 1)[0]
        start, end, strand = (int(m.group(2)), int(m.group(3)), int(m.group(4))) if m else ("", "", "")
        seq = str(rec.seq).rstrip("*").replace("*", "X")
        kept.append(SeqRecord(Seq(seq), id=pid, description=""))
        rows.append({"protein_id": pid, "source_sequence": source, "start": start, "end": end,
                     "strand": strand, "protein_length_aa": len(seq), "gene_call": "Prodigal"})
    if not kept:
        raise PipelineError("Prodigal did not find any genes. If these are transcripts / mRNA "
                            "(or from a plant, animal or fungus), choose the transcript option.")
    SeqIO.write(kept, str(out_faa), "fasta")
    log(f"      Prodigal predicted {len(kept)} proteins")
    return pd.DataFrame(rows)


def longest_orf(seq: str, table: int = 1, min_aa: int = 100, both_strands: bool = True):
    """Longest ORF that starts with Met, searched on both strands and in all 6 frames.
    Returns (protein, start, end, strand, has_stop_codon) or None."""
    L = len(seq)
    best = None
    strands = [(1, seq)]
    if both_strands:
        strands.append((-1, str(Seq(seq).reverse_complement())))
    for strand, s in strands:
        for frame in range(3):
            usable = (L - frame) // 3 * 3
            if usable < 3:
                continue
            aa = str(Seq(s[frame:frame + usable]).translate(table=table))
            parts = aa.split("*")
            offset = 0
            for i, part in enumerate(parts):
                m = part.find("M")
                if m != -1:
                    orf = part[m:]
                    if len(orf) >= min_aa and (best is None or len(orf) > len(best[0])):
                        has_stop = i < len(parts) - 1
                        aa_start = offset + m
                        nt_start = frame + 3 * aa_start + 1
                        nt_end = frame + 3 * (aa_start + len(orf) + (1 if has_stop else 0))
                        if strand == -1:
                            nt_start, nt_end = L - nt_end + 1, L - nt_start + 1
                        best = (orf, nt_start, nt_end, strand, has_stop)
                offset += len(part) + 1
    return best


def prepare_transcripts(clean_fasta, out_faa, settings: RunSettings, log: Log,
                        stats: dict) -> pd.DataFrame:
    """Translate the longest ORF of every transcript / CDS. Returns the gene_map table."""
    kept, rows = [], []
    for rec in SeqIO.parse(str(clean_fasta), "fasta"):
        res = longest_orf(str(rec.seq), settings.genetic_code, settings.min_orf_aa,
                          settings.orf_both_strands)
        if res is None:
            stats["no_orf"] += 1
            continue
        prot, start, end, strand, complete = res
        pid = f"{rec.id}_ORF1"
        kept.append(SeqRecord(Seq(prot), id=pid, description=""))
        rows.append({"protein_id": pid, "source_sequence": rec.id, "start": start, "end": end,
                     "strand": strand, "protein_length_aa": len(prot),
                     "gene_call": "Longest ORF" if complete else "Longest ORF (no stop codon)"})
    if not kept:
        raise PipelineError(f"No open reading frame of at least {settings.min_orf_aa} amino acids "
                            "was found. Try a smaller minimum ORF length in the settings, or check "
                            "the sequence type.")
    SeqIO.write(kept, str(out_faa), "fasta")
    log(f"      {len(kept)} ORFs found ({stats['no_orf']} sequences had none long enough)")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- eggNOG-mapper
def build_emapper_command(proteins_faa, out_dir, tmp_dir, data_dir, s: RunSettings) -> list:
    cmd = ["emapper.py", "-i", proteins_faa, "-o", "eggnog", "--output_dir", out_dir,
           "--temp_dir", tmp_dir, "--data_dir", data_dir, "-m", "diamond",
           "--itype", "proteins", "--cpu", s.cpu, "--sensmode", s.sensmode,
           "--evalue", s.evalue, "--score", s.min_score,
           "--seed_ortholog_evalue", s.evalue, "--seed_ortholog_score", s.min_score,
           "--target_orthologs", s.target_orthologs, "--go_evidence", s.go_evidence,
           "--override"]
    if s.min_pident > 0:
        cmd += ["--pident", s.min_pident]
    if s.min_query_cover > 0:
        cmd += ["--query_cover", s.min_query_cover]
    if s.min_subject_cover > 0:
        cmd += ["--subject_cover", s.min_subject_cover]
    if s.tax_scope:
        cmd += ["--tax_scope", s.tax_scope]
    return [str(c) for c in cmd]


def emapper_version_line() -> str:
    try:
        out = subprocess.run(["emapper.py", "--version"], capture_output=True, text=True,
                             timeout=120)
        for line in (out.stdout + "\n" + out.stderr).splitlines():
            if line.strip().startswith("emapper-"):
                return line.strip()
    except Exception:
        pass
    return ""


# --------------------------------------------------------------------------- protecting a shared server
def count_sequences(path) -> int:
    """Number of FASTA records (fast: counts header lines)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return sum(1 for line in fh if line.startswith(">"))


def _lock_is_stale(lock: Path) -> bool:
    try:
        pid, started = lock.read_text().split()[:2]
        if time.time() - float(started) > 24 * 3600:
            return True
        os.kill(int(pid), 0)             # raises if that process no longer exists
        return False
    except ProcessLookupError:
        return True
    except PermissionError:              # process exists but belongs to someone else
        return False
    except Exception:
        return True


@contextmanager
def analysis_lock(runs_dir):
    """Only ONE analysis may run at a time, so several users cannot exhaust the server's memory."""
    runs_dir = Path(runs_dir)
    runs_dir.mkdir(parents=True, exist_ok=True)
    lock = runs_dir / ".analysis.lock"
    for _attempt in range(2):
        try:
            fd = os.open(str(lock), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            os.write(fd, f"{os.getpid()} {time.time()}".encode())
            os.close(fd)
            break
        except FileExistsError:
            if _lock_is_stale(lock):
                lock.unlink(missing_ok=True)
                continue
            raise PipelineError("Another analysis is running on this server right now. "
                                "Please wait a few minutes and press Run again.")
    else:
        raise PipelineError("Could not start the analysis (the server is busy). Please try again.")
    try:
        yield
    finally:
        lock.unlink(missing_ok=True)


def prodigal_version_line() -> str:
    try:
        out = subprocess.run(["prodigal", "-v"], capture_output=True, text=True, timeout=30)
        for line in (out.stdout + "\n" + out.stderr).splitlines():
            if line.strip().lower().startswith("prodigal"):
                return line.strip()
    except Exception:
        pass
    return ""


def find_run(run_id: str, runs_dir=RUNS_DIR):
    """Folder of a finished run from its Run ID (or None). Refuses anything that is not a plain name."""
    run_id = (run_id or "").strip()
    if not run_id or Path(run_id).name != run_id or run_id.startswith("."):
        return None
    d = Path(runs_dir) / run_id
    f = d / "run_info.json"
    try:
        if f.exists() and json.loads(f.read_text()).get("status") == "finished":
            return d
    except Exception:
        pass
    return None


# --------------------------------------------------------------------------- the whole run
def _write_info(run_dir: Path, info: dict) -> None:
    (run_dir / "run_info.json").write_text(json.dumps(info, indent=2, default=str))


def run_full_pipeline(input_fasta, run_name: str, settings: RunSettings,
                      data_dir=DEFAULT_DB_DIR, runs_dir=RUNS_DIR, log: Log = print,
                      max_sequences: Optional[int] = None) -> Path:
    """Run everything and return the folder that holds all results of this run.
    Refuses files above the sequence limit and runs only one analysis at a time."""
    limit = appconfig.MAX_SEQUENCES if max_sequences is None else max_sequences
    try:
        n = count_sequences(input_fasta)
    except OSError as exc:
        raise PipelineError(f"Could not read the input file: {exc}")
    if limit and n > limit:
        raise PipelineError(f"This file has {n:,} sequences, but the limit per run on this server is "
                            f"{limit:,}. Please split the file into smaller parts.")
    with analysis_lock(runs_dir):
        return _run_locked(input_fasta, run_name, settings, data_dir, runs_dir, log)


def _run_locked(input_fasta, run_name, settings, data_dir, runs_dir, log) -> Path:
    input_fasta, data_dir, runs_dir = Path(input_fasta), Path(data_dir).expanduser(), Path(runs_dir)
    problems = check_environment(settings.input_type, data_dir)
    if problems:
        raise PipelineError("\n".join(problems))

    run_dir = runs_dir / f"{datetime.now():%Y%m%d_%H%M%S}_{slugify(run_name)}_{secrets.token_hex(6)}"
    run_dir.mkdir(parents=True, exist_ok=False)
    info = {"run_name": run_name, "run_id": run_dir.name,
            "created": f"{datetime.now():%Y-%m-%d %H:%M}", "status": "running",
            "input_file": input_fasta.name, "input_type": settings.input_type,
            "input_type_label": INPUT_TYPE_LABELS.get(settings.input_type, settings.input_type),
            "settings": settings.to_dict(), "database_dir": str(data_dir),
            "app_version": appconfig.APP_VERSION}
    _write_info(run_dir, info)
    t0 = time.time()
    total = 2 if settings.input_type == INPUT_PROTEIN else 3
    try:
        shutil.copy2(input_fasta, run_dir / "input_original.fasta")
        proteins = run_dir / "proteins.faa"

        log(f"[1/{total}] Checking the input sequences...")
        if settings.input_type == INPUT_PROTEIN:
            gene_map, stats = prepare_proteins(input_fasta, proteins, log)
        else:
            clean = run_dir / "clean_dna.fasta"
            stats = prepare_dna(input_fasta, clean, log)
            if settings.input_type == INPUT_PROKARYOTE:
                log(f"[2/{total}] Predicting genes with Prodigal...")
                gene_map = run_prodigal(clean, run_dir, proteins, log, stats)
            else:
                log(f"[2/{total}] Translating the longest ORF of each sequence...")
                gene_map = prepare_transcripts(clean, proteins, settings, log, stats)
        gene_map.to_csv(run_dir / "gene_map.tsv", sep="\t", index=False)

        log(f"[{total}/{total}] Running eggNOG-mapper on {len(gene_map)} proteins "
            "(the first run after starting the computer is slowest: the database must be loaded)...")
        out_dir, tmp_dir = run_dir / "eggnog", run_dir / "tmp"
        out_dir.mkdir(), tmp_dir.mkdir()
        run_command(build_emapper_command(proteins, out_dir, tmp_dir, data_dir, settings), log)
        shutil.rmtree(tmp_dir, ignore_errors=True)

        info.update(status="finished", emapper_version=emapper_version_line(),
                    prodigal_version=(prodigal_version_line() if settings.input_type == INPUT_PROKARYOTE else ""),
                    counts={"input_sequences": stats["input"], "sequences_used": stats["kept"],
                            "proteins": len(gene_map), **{k: v for k, v in stats.items()
                                                         if k not in ("input", "kept")}},
                    seconds=round(time.time() - t0, 1))
        _write_info(run_dir, info)
        log(f"Finished in {info['seconds']} s. Results saved in {run_dir}")
        return run_dir
    except Exception as exc:
        info.update(status="failed", error=str(exc))
        _write_info(run_dir, info)
        raise


# --------------------------------------------------------------------------- automatic deletion
def _newest_mtime(path: Path) -> float:
    """Most recent modification time of a folder and everything inside it."""
    newest = path.stat().st_mtime
    for p in path.rglob("*"):
        try:
            newest = max(newest, p.stat().st_mtime)
        except OSError:
            pass
    return newest


def cleanup_old_data(runs_dir=RUNS_DIR, upload_dir=None, days: Optional[float] = None) -> dict:
    """Delete saved runs and uploaded files that are older than `days` days.

    Runs: only folders that contain a run_info.json are touched; the age is counted from the time that file was
    last written (when the run finished). Uploads: the per-visitor folders inside `upload_dir`; the age is the
    newest file inside. `days` defaults to appconfig.RETENTION_DAYS; 0 or less means nothing is deleted."""
    days = appconfig.RETENTION_DAYS if days is None else days
    upload_dir = appconfig.UPLOAD_DIR if upload_dir is None else upload_dir
    removed = {"runs": 0, "uploads": 0}
    if not days or days <= 0:
        return removed
    cutoff = time.time() - days * 86400
    runs_dir, upload_dir = Path(runs_dir), Path(upload_dir)
    if runs_dir.is_dir():
        for d in runs_dir.iterdir():
            f = d / "run_info.json"
            try:
                if d.name.startswith(".") or not d.is_dir() or not f.exists():
                    continue
                if f.stat().st_mtime < cutoff:
                    shutil.rmtree(d, ignore_errors=True)
                    removed["runs"] += 1
            except OSError:
                continue
    if upload_dir.is_dir():
        for d in upload_dir.iterdir():
            try:
                if d.name.startswith(".") or not d.is_dir():
                    continue
                if _newest_mtime(d) < cutoff:
                    shutil.rmtree(d, ignore_errors=True)
                    removed["uploads"] += 1
            except OSError:
                continue
    return removed


def list_runs(runs_dir=RUNS_DIR) -> list[dict]:
    """Finished runs, newest first."""
    out = []
    runs_dir = Path(runs_dir)
    if not runs_dir.exists():
        return out
    for d in sorted(runs_dir.iterdir(), reverse=True):
        f = d / "run_info.json"
        if f.exists():
            try:
                info = json.loads(f.read_text())
            except Exception:
                continue
            if info.get("status") == "finished":
                out.append({"path": str(d), "run_id": d.name, "name": info.get("run_name", d.name),
                            "created": info.get("created", ""),
                            "proteins": info.get("counts", {}).get("proteins")})
    return out
