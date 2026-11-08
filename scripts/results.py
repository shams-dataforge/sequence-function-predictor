"""
results.py  -  turns the raw files of ONE run into readable tables.

  * reads eggNOG-mapper's output (annotations + seed-ortholog statistics)
  * adds descriptions: GO term names/definitions, KEGG names, COG categories
  * gives every sequence a High / Medium / Low confidence level (thresholds adjustable)
  * GO / KEGG enrichment (over-representation) tests
  * builds the Excel workbook

Nothing here runs a search, so it is fast and can be repeated with new thresholds.
"""
from __future__ import annotations

import io
import json
import re
import urllib.request
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

ASPECT_LABEL = {"molecular_function": "Molecular Function",
                "biological_process": "Biological Process",
                "cellular_component": "Cellular Component"}
ASPECT_ORDER = ["Molecular Function", "Biological Process", "Cellular Component", "Unknown"]
ASPECT_COLUMN = {"molecular_function": "GO_molecular_function",
                 "biological_process": "GO_biological_process",
                 "cellular_component": "GO_cellular_component"}

COG_CATEGORIES = {
    "J": "Translation, ribosomal structure and biogenesis", "A": "RNA processing and modification",
    "K": "Transcription", "L": "Replication, recombination and repair",
    "B": "Chromatin structure and dynamics", "D": "Cell cycle control, cell division, chromosome partitioning",
    "Y": "Nuclear structure", "V": "Defense mechanisms", "T": "Signal transduction mechanisms",
    "M": "Cell wall/membrane/envelope biogenesis", "N": "Cell motility", "Z": "Cytoskeleton",
    "W": "Extracellular structures", "U": "Intracellular trafficking, secretion, and vesicular transport",
    "O": "Posttranslational modification, protein turnover, chaperones",
    "X": "Mobilome: prophages, transposons", "C": "Energy production and conversion",
    "G": "Carbohydrate transport and metabolism", "E": "Amino acid transport and metabolism",
    "F": "Nucleotide transport and metabolism", "H": "Coenzyme transport and metabolism",
    "I": "Lipid transport and metabolism", "P": "Inorganic ion transport and metabolism",
    "Q": "Secondary metabolites biosynthesis, transport and catabolism",
    "R": "General function prediction only", "S": "Function unknown",
}

GO_OBO_URL = "https://purl.obolibrary.org/obo/go/go-basic.obo"
KEGG_LISTS = {"ko": "https://rest.kegg.jp/list/ko",
              "pathway": "https://rest.kegg.jp/list/pathway",
              "module": "https://rest.kegg.jp/list/module"}

STAT_COLS = ["evalue", "bitscore", "identity_pct", "query_cov_pct", "subject_cov_pct", "coverage_pct"]
NAME_NA = "(name not available)"

# The columns shown in the "Compact" view (also used for the Excel sheets). Edit these lists to
# change what most people see; the "Full" view always shows everything.
COMPACT_SUMMARY = ["sequence_id", "source_sequence", "confidence", "evalue", "bitscore",
                   "identity_pct", "coverage_pct", "gene_name", "description",
                   "GO_molecular_function", "GO_biological_process", "GO_cellular_component",
                   "KEGG_KO", "KEGG_pathway"]
COMPACT_GO = ["sequence_id", "aspect", "GO_ID", "GO_term", "confidence", "evalue", "bitscore",
              "identity_pct", "coverage_pct"]
COMPACT_KO = ["sequence_id", "KO", "KO_name", "confidence", "evalue", "bitscore",
              "identity_pct", "coverage_pct"]
COMPACT_PATHWAY = ["sequence_id", "pathway_ID", "pathway_name", "confidence", "evalue",
                   "bitscore", "identity_pct", "coverage_pct"]


def compact(df: pd.DataFrame, cols: list) -> pd.DataFrame:
    """Keep only the listed columns (those that exist). 'source_sequence' is dropped when it
    only repeats 'sequence_id' (protein input)."""
    out = df[[c for c in cols if c in df.columns]]
    if "source_sequence" in out.columns and "sequence_id" in out.columns \
            and len(out) and (out["source_sequence"] == out["sequence_id"]).all():
        out = out.drop(columns="source_sequence")
    return out


# =========================================================================== thresholds
@dataclass(frozen=True)
class Thresholds:
    """Rules for the High / Medium / Low confidence levels. Coverage = the smaller of
    query coverage and subject coverage. Starting values are common practice, not a
    universal standard - adjust them to your own needs and report them in your methods."""
    high_evalue: float = 1e-30
    high_identity: float = 60.0
    high_coverage: float = 70.0
    med_evalue: float = 1e-10
    med_identity: float = 40.0
    med_coverage: float = 50.0

    def as_tuple(self) -> tuple:
        return (self.high_evalue, self.high_identity, self.high_coverage,
                self.med_evalue, self.med_identity, self.med_coverage)


def classify_confidence(df: pd.DataFrame, th: Thresholds, has_hit: pd.Series) -> pd.Series:
    e = pd.to_numeric(df["evalue"], errors="coerce")
    pid = pd.to_numeric(df["identity_pct"], errors="coerce")
    cov = pd.concat([pd.to_numeric(df["query_cov_pct"], errors="coerce"),
                     pd.to_numeric(df["subject_cov_pct"], errors="coerce")],
                    axis=1).min(axis=1, skipna=False)
    avail = e.notna() | pid.notna() | cov.notna()

    def meets(e_th, id_th, cov_th):
        # a missing statistic is ignored; at least one statistic must exist
        return avail & (e.isna() | (e <= e_th)) & (pid.isna() | (pid >= id_th)) & \
               (cov.isna() | (cov >= cov_th))

    hi = meets(th.high_evalue, th.high_identity, th.high_coverage)
    med = meets(th.med_evalue, th.med_identity, th.med_coverage)
    out = np.where(~has_hit, "No hit",
          np.where(~avail, "Not scored",
          np.where(hi, "High", np.where(med, "Medium", "Low"))))
    return pd.Series(out, index=df.index)


# =========================================================================== GO / KEGG lookups
def ensure_go_obo(path) -> Path:
    """Make sure go-basic.obo exists (downloads it once if it is missing)."""
    path = Path(path)
    if path.exists() and path.stat().st_size > 1_000_000:
        return path
    tmp = path.with_suffix(".tmp")
    req = urllib.request.Request(GO_OBO_URL, headers={"User-Agent": "seq2func/1.0"})
    with urllib.request.urlopen(req, timeout=120) as r, open(tmp, "wb") as f:
        f.write(r.read())
    if tmp.stat().st_size < 1_000_000:
        tmp.unlink(missing_ok=True)
        raise RuntimeError("Downloaded GO file looks incomplete.")
    tmp.replace(path)
    return path


def load_go_dag(path):
    from goatools.obo_parser import GODag
    return GODag(str(path), optional_attrs=["defn"], load_obsolete=True, prt=None)


def _clean_definition(text: str) -> str:
    m = re.match(r'^"(.*?)"(?:\s*\[.*\])?\s*$', text or "", flags=re.S)
    return (m.group(1) if m else (text or "")).strip()


def go_info(go_dag, go_id: str) -> dict:
    """name / aspect / definition for one GO ID (never raises)."""
    term = go_dag.get(go_id) if go_dag is not None else None
    if term is None:
        return {"name": "(not found in the loaded ontology)", "aspect": "Unknown",
                "aspect_key": "", "definition": "", "obsolete": False}
    name = term.name or ""
    obsolete = bool(getattr(term, "is_obsolete", False))
    if obsolete and not name.lower().startswith("obsolete"):
        name = "obsolete " + name
    ns = getattr(term, "namespace", "") or ""
    return {"name": name, "aspect": ASPECT_LABEL.get(ns, "Unknown"), "aspect_key": ns,
            "definition": _clean_definition(getattr(term, "defn", "")), "obsolete": obsolete}


def go_parents(go_dag, go_id: str) -> list[str]:
    term = go_dag.get(go_id) if go_dag is not None else None
    if term is None:
        return []
    return sorted(f"{p.item_id} ({p.name})" for p in term.parents)


def _read_kegg_file(path: Path) -> dict:
    out = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            out[k.split(":", 1)[-1].strip()] = v.strip()
    return out


def load_kegg_names(cache_dir, refresh: bool = False, timeout: int = 90) -> dict:
    """KEGG ID -> readable name (KO, pathway, module). Downloaded once, then cached.
    KEGG's web service is free for academic use - check its terms before public hosting."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    result = {"ko": {}, "pathway": {}, "module": {}, "messages": []}
    for key, url in KEGG_LISTS.items():
        f = cache_dir / f"kegg_{key}.tsv"
        if refresh or not f.exists() or f.stat().st_size == 0:
            try:
                req = urllib.request.Request(url, headers={"User-Agent": "seq2func/1.0"})
                with urllib.request.urlopen(req, timeout=timeout) as r:
                    f.write_bytes(r.read())
            except Exception as exc:
                result["messages"].append(f"Could not download KEGG {key} names ({exc}).")
                continue
        try:
            result[key] = _read_kegg_file(f)
        except Exception as exc:
            result["messages"].append(f"Could not read KEGG {key} names ({exc}).")
    return result


# =========================================================================== reading eggNOG output
def _read_hash_table(path: Path, first_header_col: set):
    header, rows, comments = None, [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip("\r\n")
            if not line.strip():
                continue
            if line.startswith("##"):
                comments.append(line)
            elif line.startswith("#"):
                cols = line[1:].split("\t")
                if header is None and cols[0] in first_header_col:
                    header = cols
                else:
                    comments.append(line)
            else:
                rows.append(line.split("\t"))
    return header, rows, comments


def _frame(header, rows) -> pd.DataFrame:
    n = len(header)
    fixed = [(r + [""] * (n - len(r)))[:n] for r in rows]
    return pd.DataFrame(fixed, columns=header)


ANN_EXPECTED = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
                "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
                "KEGG_Pathway", "KEGG_Module", "PFAMs"]


def read_annotations(path) -> tuple[pd.DataFrame, dict]:
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=ANN_EXPECTED), {}
    header, rows, comments = _read_hash_table(path, {"query"})
    if header is None:
        return pd.DataFrame(columns=ANN_EXPECTED), {}
    df = _frame(header, rows)
    for c in ANN_EXPECTED:
        if c not in df.columns:
            df[c] = ""
    for c in df.columns:
        if c not in ("evalue", "score"):
            df[c] = df[c].replace("-", "")
    df["evalue"] = pd.to_numeric(df["evalue"], errors="coerce")
    df["score"] = pd.to_numeric(df["score"], errors="coerce")
    meta = {"emapper_comment": next((c for c in comments if "emapper-" in c), "")}
    return df.drop_duplicates("query"), meta


def read_seed_orthologs(path) -> pd.DataFrame:
    cols = ["qseqid", "sseqid", "evalue", "bitscore", "qstart", "qend", "sstart", "send",
            "pident", "qcov", "scov"]
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=cols)
    header, rows, _ = _read_hash_table(path, {"qseqid"})
    if header is None:
        header = cols
    df = _frame(header, rows)
    for c in cols:
        if c not in df.columns:
            df[c] = np.nan
    for c in ("evalue", "bitscore", "pident", "qcov", "scov"):
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df.drop_duplicates("qseqid")


# =========================================================================== parsing helpers
def _split(s) -> list[str]:
    return [x.strip() for x in str(s).split(",") if x.strip() and x.strip() != "-"]


def parse_go(s) -> list[str]:
    return list(dict.fromkeys(x for x in _split(s) if x.startswith("GO:")))


def parse_ko(s) -> list[str]:
    return list(dict.fromkeys(x.replace("ko:", "") for x in _split(s) if re.search(r"K\d{5}", x)))


def parse_pathways(s) -> list[str]:
    out = []
    for x in _split(s):
        m = re.search(r"(\d{5})$", x)
        if m:
            out.append("map" + m.group(1))
    return list(dict.fromkeys(out))


def parse_modules(s) -> list[str]:
    return list(dict.fromkeys(x for x in _split(s) if re.fullmatch(r"M\d{5}", x)))


def _with_name(i: str, names: dict) -> str:
    n = names.get(i, "")
    return f"{i} ({n})" if n else i


def _cog_text(letters: str) -> str:
    letters = "".join(ch for ch in str(letters) if ch.isalpha())
    return "; ".join(f"{ch}: {COG_CATEGORIES.get(ch, '?')}" for ch in letters)


# =========================================================================== building the tables
@dataclass
class ResultTables:
    summary: pd.DataFrame
    go_long: pd.DataFrame
    go_dict: pd.DataFrame
    ko_long: pd.DataFrame
    ko_dict: pd.DataFrame
    pathway_long: pd.DataFrame
    pathway_dict: pd.DataFrame
    info: dict
    counts: dict
    notes: list = field(default_factory=list)
    thresholds: Thresholds = field(default_factory=Thresholds)


def build_tables(run_dir, thresholds: Thresholds, go_dag, kegg: Optional[dict],
                 show_obsolete: bool = False) -> ResultTables:
    run_dir = Path(run_dir)
    info = json.loads((run_dir / "run_info.json").read_text())
    kegg = kegg or {"ko": {}, "pathway": {}, "module": {}, "messages": []}
    notes = list(kegg.get("messages", []))
    if go_dag is None:
        notes.append("The Gene Ontology file could not be loaded - GO term names are not available.")

    gm = pd.read_csv(run_dir / "gene_map.tsv", sep="\t", dtype=str, keep_default_na=False)
    for c in ("start", "end", "protein_length_aa"):
        gm[c] = pd.to_numeric(gm[c], errors="coerce")
    gm["strand"] = pd.to_numeric(gm["strand"], errors="coerce")

    ann, _meta = read_annotations(run_dir / "eggnog" / "eggnog.emapper.annotations")
    seed = read_seed_orthologs(run_dir / "eggnog" / "eggnog.emapper.seed_orthologs")
    if ann.empty:
        notes.append("eggNOG-mapper found no ortholog matches with the chosen cut-offs.")
    if not ann.empty and seed["pident"].isna().all():
        notes.append("Identity / coverage values were not found, so confidence levels use the "
                     "E-value only.")

    a = ann.rename(columns={"query": "protein_id", "Preferred_name": "gene_name",
                            "Description": "description", "max_annot_lvl": "taxonomic_level",
                            "score": "bitscore", "COG_category": "cog_category",
                            "PFAMs": "PFAM_domains", "GOs": "_GOs", "KEGG_ko": "_KO",
                            "KEGG_Pathway": "_PATH", "KEGG_Module": "_MOD"})
    a["_hit"] = True
    keep = ["protein_id", "_hit", "evalue", "bitscore", "seed_ortholog", "taxonomic_level",
            "gene_name", "description", "cog_category", "EC", "PFAM_domains",
            "_GOs", "_KO", "_PATH", "_MOD"]
    s = seed.rename(columns={"qseqid": "protein_id", "pident": "identity_pct",
                             "qcov": "query_cov_pct", "scov": "subject_cov_pct"})
    s = s[["protein_id", "identity_pct", "query_cov_pct", "subject_cov_pct"]]

    sm = gm.merge(a[keep], on="protein_id", how="left").merge(s, on="protein_id", how="left")
    sm["_hit"] = sm["_hit"].eq(True)
    for c in ("seed_ortholog", "taxonomic_level", "gene_name", "description", "cog_category",
              "EC", "PFAM_domains", "_GOs", "_KO", "_PATH", "_MOD"):
        sm[c] = sm[c].fillna("")
    sm = sm.rename(columns={"protein_id": "sequence_id"})
    sm["coverage_pct"] = sm[["query_cov_pct", "subject_cov_pct"]].min(axis=1, skipna=False)

    go_lists = sm["_GOs"].map(parse_go)
    ko_lists = sm["_KO"].map(parse_ko)
    path_lists = sm["_PATH"].map(parse_pathways)
    mod_lists = sm["_MOD"].map(parse_modules)

    memo: dict = {}

    def ginfo(g):
        if g not in memo:
            memo[g] = go_info(go_dag, g)
        return memo[g]

    hidden_obsolete = 0
    if not show_obsolete:          # obsolete terms are not valid annotations: hide them by default
        def _drop_obsolete(ids):
            nonlocal hidden_obsolete
            keep_ids = [g for g in ids if not ginfo(g)["obsolete"]]
            hidden_obsolete += len(ids) - len(keep_ids)
            return keep_ids
        go_lists = go_lists.map(_drop_obsolete)

    # wide GO columns:  "GO:0005737 (cytoplasm); GO:..."
    for key, col in ASPECT_COLUMN.items():
        sm[col] = [
            "; ".join(f"{g} ({ginfo(g)['name']})" for g in ids if ginfo(g)["aspect_key"] == key)
            for ids in go_lists]
    sm["GO_unclassified"] = ["; ".join(f"{g} ({ginfo(g)['name']})" for g in ids
                                        if ginfo(g)["aspect_key"] == "") for ids in go_lists]
    sm["n_GO_terms"] = go_lists.map(len)
    sm["KEGG_KO"] = [("; ".join(_with_name(k, kegg["ko"]) for k in ks)) for ks in ko_lists]
    sm["KEGG_pathway"] = [("; ".join(_with_name(p, kegg["pathway"]) for p in ps)) for ps in path_lists]
    sm["KEGG_module"] = [("; ".join(_with_name(m, kegg["module"]) for m in ms)) for ms in mod_lists]
    sm["COG_description"] = sm["cog_category"].map(_cog_text)

    has_ann = (sm["n_GO_terms"] > 0) | ko_lists.map(bool) | path_lists.map(bool)
    sm["annotation_status"] = np.where(~sm["_hit"], "No ortholog hit",
                              np.where(has_ann, "Annotated", "Ortholog hit, no GO/KEGG"))
    sm["confidence"] = classify_confidence(sm, thresholds, sm["_hit"])

    summary_cols = ["sequence_id", "source_sequence", "start", "end", "strand",
                    "protein_length_aa", "gene_call", "annotation_status", "confidence",
                    "evalue", "bitscore", "identity_pct", "query_cov_pct", "subject_cov_pct",
                    "coverage_pct", "seed_ortholog", "taxonomic_level", "gene_name", "description",
                    "cog_category", "COG_description", "EC", "PFAM_domains",
                    "GO_molecular_function", "GO_biological_process", "GO_cellular_component",
                    "GO_unclassified", "n_GO_terms", "KEGG_KO", "KEGG_pathway", "KEGG_module"]
    summary = sm[summary_cols].copy()
    for c in ("start", "end", "strand", "protein_length_aa"):
        summary[c] = pd.to_numeric(summary[c], errors="coerce").astype("Int64")
    stats = summary[["sequence_id"] + STAT_COLS[:0] + ["confidence"] + STAT_COLS +
                    ["seed_ortholog"]].copy()

    # ---------------------------------------------------------------- GO long table
    rows = []
    for pid, src, ids in zip(sm["sequence_id"], sm["source_sequence"], go_lists):
        for g in ids:
            i = ginfo(g)
            rows.append((pid, src, i["aspect"], g, i["name"], i["definition"],
                         "yes" if i["obsolete"] else ""))
    go_long = pd.DataFrame(rows, columns=["sequence_id", "source_sequence", "aspect", "GO_ID",
                                          "GO_term", "definition", "obsolete"])
    go_long = go_long.drop_duplicates(["sequence_id", "GO_ID"])
    go_long = go_long.merge(stats, on="sequence_id", how="left")
    go_long["aspect"] = pd.Categorical(go_long["aspect"], categories=ASPECT_ORDER, ordered=True)
    go_long = go_long.sort_values(["aspect", "sequence_id", "GO_ID"]).reset_index(drop=True)
    go_long["aspect"] = go_long["aspect"].astype(str)

    n_obs = int((go_long["obsolete"] == "yes").sum()) if len(go_long) else 0
    n_unknown = int((go_long["aspect"] == "Unknown").sum()) if len(go_long) else 0
    if hidden_obsolete:
        notes.append(f"{hidden_obsolete} GO annotations use terms that are obsolete in the current "
                     "ontology and are hidden. Tick 'Show obsolete GO terms' in the sidebar to see them.")
    if n_obs:
        notes.append(f"{n_obs} GO annotations use terms that are obsolete in the current ontology; "
                     "they are kept and flagged in the 'obsolete' column.")
    if n_unknown:
        notes.append(f"{n_unknown} GO annotations were not found in the loaded ontology file "
                     "(it may be older than eggNOG's data).")

    def dictionary(df, key, cols_first):
        if df.empty:
            return pd.DataFrame(columns=cols_first + ["n_sequences", "n_high_confidence"])
        g = df.groupby(key, sort=False)
        d = g.agg(n_sequences=("sequence_id", "nunique"),
                  n_high_confidence=("confidence", lambda x: int((x == "High").sum())))
        first = g[[c for c in cols_first if c != key]].first()
        return first.join(d).reset_index()[cols_first + ["n_sequences", "n_high_confidence"]] \
            .sort_values("n_sequences", ascending=False).reset_index(drop=True)

    go_dict = dictionary(go_long, "GO_ID", ["GO_ID", "GO_term", "aspect", "definition"])

    # ---------------------------------------------------------------- KEGG long tables
    rows = [(pid, src, k, kegg["ko"].get(k, "")) for pid, src, ks in
            zip(sm["sequence_id"], sm["source_sequence"], ko_lists) for k in ks]
    ko_long = pd.DataFrame(rows, columns=["sequence_id", "source_sequence", "KO", "KO_name"])
    ko_long = ko_long.merge(stats, on="sequence_id", how="left")
    ko_long["KO_name"] = ko_long["KO_name"].replace("", NAME_NA)
    ko_dict = dictionary(ko_long, "KO", ["KO", "KO_name"])

    rows = [(pid, src, p, kegg["pathway"].get(p, "")) for pid, src, ps in
            zip(sm["sequence_id"], sm["source_sequence"], path_lists) for p in ps]
    pathway_long = pd.DataFrame(rows, columns=["sequence_id", "source_sequence", "pathway_ID",
                                               "pathway_name"])
    pathway_long = pathway_long.merge(stats, on="sequence_id", how="left")
    pathway_long["pathway_name"] = pathway_long["pathway_name"].replace("", NAME_NA)
    pathway_dict = dictionary(pathway_long, "pathway_ID", ["pathway_ID", "pathway_name"])

    if not kegg["pathway"] and len(pathway_long):
        notes.append("KEGG pathway names are not available offline - showing IDs only.")

    conf = summary["confidence"].value_counts()
    counts = {
        "input_sequences": info.get("counts", {}).get("input_sequences"),
        "proteins_analysed": int(len(summary)),
        "with_ortholog_hit": int(sm["_hit"].sum()),
        "annotated": int((summary["annotation_status"] == "Annotated").sum()),
        "with_GO": int((summary["n_GO_terms"] > 0).sum()),
        "with_GO_MF": int((summary["GO_molecular_function"] != "").sum()),
        "with_GO_BP": int((summary["GO_biological_process"] != "").sum()),
        "with_GO_CC": int((summary["GO_cellular_component"] != "").sum()),
        "with_KO": int(ko_lists.map(bool).sum()),
        "with_pathway": int(path_lists.map(bool).sum()),
        "high": int(conf.get("High", 0)), "medium": int(conf.get("Medium", 0)),
        "low": int(conf.get("Low", 0)), "no_hit": int(conf.get("No hit", 0)),
    }
    return ResultTables(summary, go_long, go_dict, ko_long, ko_dict, pathway_long, pathway_dict,
                        info, counts, notes, thresholds)


# =========================================================================== enrichment
def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg false-discovery-rate correction."""
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.minimum(ranked, 1.0)
    return q


def _enrich(assoc: dict, study: set, min_count: int):
    """One-sided hypergeometric test (same as Fisher's exact) for every term.
    assoc: sequence -> set of terms.  Population = sequences that have >= 1 term."""
    from scipy.stats import hypergeom
    pop = set(assoc)
    study = study & pop
    M, N = len(pop), len(study)
    if M == 0 or N == 0:
        return pd.DataFrame(), M, N
    pop_counts = Counter(t for terms in assoc.values() for t in terms)
    study_counts = Counter(t for g in study for t in assoc[g])
    rows = []
    for t, k in study_counts.items():
        if k < min_count:
            continue
        n = pop_counts[t]
        p = float(hypergeom.sf(k - 1, M, n, N))
        rows.append({"term": t, "study_count": k, "study_size": N, "population_count": n,
                     "population_size": M, "fold_enrichment": (k / N) / (n / M), "p_value": p})
    return pd.DataFrame(rows), M, N


@dataclass
class EnrichmentResult:
    go: pd.DataFrame
    kegg: pd.DataFrame
    study_requested: int
    study_matched: int
    not_found: list
    fdr_cutoff: float
    min_count: int
    notes: list = field(default_factory=list)


def resolve_study_ids(tables: ResultTables, ids: list[str]):
    """Match pasted IDs to sequence IDs (or to source sequence IDs)."""
    sm = tables.summary
    valid = set(sm["sequence_id"])
    by_source = sm.groupby("source_sequence")["sequence_id"].apply(list).to_dict()
    matched, missing = set(), []
    for i in ids:
        if i in valid:
            matched.add(i)
        elif i in by_source:
            matched.update(by_source[i])
        else:
            missing.append(i)
    return matched, missing


def run_enrichment(tables: ResultTables, go_dag, ids: list[str], min_count: int = 2,
                   fdr_cutoff: float = 0.05) -> EnrichmentResult:
    ids = list(dict.fromkeys(i for i in ids if i))
    matched, missing = resolve_study_ids(tables, ids)
    notes = []

    # ---- GO (annotations are propagated to parent terms: the "true path rule")
    anc_cache: dict = {}

    def with_ancestors(g):
        if g not in anc_cache:
            term = go_dag.get(g) if go_dag is not None else None
            anc = set(term.get_all_parents()) if term is not None else set()
            anc_cache[g] = {g} | anc
        return anc_cache[g]

    gl = tables.go_long
    gl = gl[(gl["obsolete"] != "yes") & (gl["aspect"] != "Unknown")] if len(gl) else gl
    frames = []
    for aspect in ASPECT_ORDER[:3]:
        sub = gl[gl["aspect"] == aspect] if len(gl) else gl
        assoc: dict = {}
        for pid, g in zip(sub["sequence_id"], sub["GO_ID"]):
            assoc.setdefault(pid, set()).update(with_ancestors(g))
        df, _M, _N = _enrich(assoc, matched, min_count)
        if df.empty:
            continue
        df["FDR"] = bh_fdr(df["p_value"].values)
        df.insert(0, "aspect", aspect)
        frames.append(df)
    if frames:
        go_res = pd.concat(frames, ignore_index=True)
        go_res["GO_ID"] = go_res["term"]
        go_res["GO_term"] = [go_info(go_dag, t)["name"] for t in go_res["term"]]
        go_res["significant"] = go_res["FDR"] <= fdr_cutoff
        go_res = go_res[["aspect", "GO_ID", "GO_term", "study_count", "study_size",
                         "population_count", "population_size", "fold_enrichment",
                         "p_value", "FDR", "significant"]].sort_values(["FDR", "p_value"])
    else:
        go_res = pd.DataFrame()

    # ---- KEGG pathways
    pl = tables.pathway_long
    assoc = {}
    for pid, p in zip(pl["sequence_id"], pl["pathway_ID"]):
        assoc.setdefault(pid, set()).add(p)
    kdf, _M, _N = _enrich(assoc, matched, min_count)
    if not kdf.empty:
        kdf["FDR"] = bh_fdr(kdf["p_value"].values)
        names = dict(zip(pl["pathway_ID"], pl["pathway_name"]))
        kdf.insert(0, "pathway_ID", kdf["term"])
        kdf.insert(1, "pathway_name", [names.get(t, "") for t in kdf["term"]])
        kdf["significant"] = kdf["FDR"] <= fdr_cutoff
        kdf = kdf.drop(columns=["term"]).sort_values(["FDR", "p_value"])

    n_annot = len(set(tables.go_long["sequence_id"])) if len(tables.go_long) else 0
    if matched and n_annot and len(matched) >= 0.9 * n_annot:
        notes.append("Your group of interest is almost the whole background - enrichment "
                     "results will not be informative.")
    if not matched:
        notes.append("None of the IDs you gave matched sequences in this run.")
    return EnrichmentResult(go_res, kdf, len(ids), len(matched), missing, fdr_cutoff,
                            min_count, notes)


# =========================================================================== methods text
_SCOPE_NAMES = {"bacteria": "Bacteria", "archaea": "Archaea", "eukaryota": "Eukaryota",
                "33090,2759,1": "Viridiplantae (plants)", "33208,2759,1": "Metazoa (animals)",
                "4751,2759,1": "Fungi"}


def methods_text(t: "ResultTables", go_version: str = "", enr: "Optional[EnrichmentResult]" = None,
                 app_name: str = "Sequence Function Predictor") -> str:
    """A ready-to-paste Methods paragraph built from the settings that were really used."""
    info, th = t.info, t.thresholds
    s = info.get("settings", {}) or {}

    def g(key, default):
        v = s.get(key)
        return default if v is None else v

    ver = info.get("emapper_version", "") or ""
    m = re.search(r"emapper-([\d.]+)", ver)
    d = re.search(r"Installed eggNOG DB version:\s*([\d.]+)", ver)
    emapper = f"eggNOG-mapper v{m.group(1)}" if m else "eggNOG-mapper"
    db = f"the eggNOG {d.group(1)} database" if d else "the eggNOG database"

    typ = info.get("input_type", "")
    if typ == "protein":
        first = "Protein sequences were analysed as provided."
    elif typ == "prokaryote":
        pv = (info.get("prodigal_version", "") or "").split(":")[0].replace(" V", " v") or "Prodigal"
        first = f"Genes were predicted with {pv} in anonymous (metagenomic) mode and translated into proteins."
    else:
        strands = "both strands" if g("orf_both_strands", True) else "the forward strand"
        first = (f"The longest open reading frame (at least {g('min_orf_aa', 100)} amino acids; genetic code "
                 f"table {g('genetic_code', 1)}; {strands} searched) of each sequence was translated into a protein.")

    filt = [f"E-value <= {float(g('evalue', 0.001)):g}", f"bit score >= {float(g('min_score', 60)):g}"]
    if float(g("min_pident", 0)) > 0:
        filt.append(f"identity >= {float(g('min_pident', 0)):g}%")
    if float(g("min_query_cover", 0)) > 0:
        filt.append(f"query coverage >= {float(g('min_query_cover', 0)):g}%")
    if float(g("min_subject_cover", 0)) > 0:
        filt.append(f"subject coverage >= {float(g('min_subject_cover', 0)):g}%")
    scope = g("tax_scope", None)
    scope_txt = ("chosen automatically for each sequence" if not scope
                 else _SCOPE_NAMES.get(str(scope), str(scope)))

    text = (f"{first} Proteins were annotated by orthology with {emapper} using {db}. Sequence searches were "
            f"performed with DIAMOND in '{g('sensmode', 'fast')}' mode ({', '.join(filt)}). The taxonomic scope was "
            f"{scope_txt}; orthologs of type '{g('target_orthologs', 'all')}' were used to transfer annotations, and "
            f"GO terms were transferred using '{g('go_evidence', 'non-electronic')}' evidence. "
            f"Each protein was given a confidence level from its best ortholog match: High (E-value <= "
            f"{th.high_evalue:g}, identity >= {th.high_identity:g}%, coverage >= {th.high_coverage:g}%), Medium "
            f"(E-value <= {th.med_evalue:g}, identity >= {th.med_identity:g}%, coverage >= {th.med_coverage:g}%) or "
            f"Low (passing the search thresholds but not the Medium rule), where coverage is the smaller of the "
            f"query and subject coverage. GO term names and definitions come from the Gene Ontology"
            f"{' (' + go_version + ')' if go_version else ''} and KEGG names from the KEGG REST API.")
    if enr is not None and enr.study_matched:
        text += (f" Enrichment of GO terms and KEGG pathways in a group of {enr.study_matched} sequences was tested "
                 f"against all annotated sequences of the run with a one-sided hypergeometric test (annotations "
                 f"propagated to parent GO terms) and Benjamini-Hochberg correction (FDR <= {enr.fdr_cutoff:g}).")
    text += (f" The results are computational predictions and should be validated experimentally. "
             f"Analysis was carried out with {app_name} v{info.get('app_version', '')}.").replace(" v.", ".")
    return text


# =========================================================================== exports
def to_csv_bytes(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8-sig")


def _readme_rows(t: ResultTables, enr: Optional[EnrichmentResult]) -> list:
    info, th, st = t.info, t.thresholds, t.info.get("settings", {})
    rows = [("Sequence -> Function Predictor", "Results workbook"), ("", ""),
            ("Run name", info.get("run_name")), ("Date", info.get("created")),
            ("Run ID", info.get("run_id")), ("Software version", info.get("app_version", "")),
            ("Input file", info.get("input_file")), ("Sequence type", info.get("input_type_label")),
            ("eggNOG-mapper", info.get("emapper_version", "")),
            ("Prodigal", info.get("prodigal_version", "")), ("", ""),
            ("SEARCH SETTINGS", "")]
    for k in ("sensmode", "evalue", "min_score", "min_pident", "min_query_cover",
              "min_subject_cover", "tax_scope", "target_orthologs", "go_evidence",
              "min_orf_aa", "genetic_code"):
        rows.append((k, st.get(k)))
    rows += [("", ""), ("CONFIDENCE RULES (coverage = smaller of query / subject coverage)", ""),
             ("High", f"E-value <= {th.high_evalue:g}, identity >= {th.high_identity:g}%, "
                      f"coverage >= {th.high_coverage:g}%"),
             ("Medium", f"E-value <= {th.med_evalue:g}, identity >= {th.med_identity:g}%, "
                        f"coverage >= {th.med_coverage:g}%"),
             ("Low", "Passed the search cut-offs but not the Medium rule"),
             ("No hit", "No eggNOG ortholog found with the chosen cut-offs"), ("", ""),
             ("COLUMN GUIDE", ""),
             ("evalue", "Expected number of matches this good by chance. Smaller = more reliable."),
             ("bitscore", "Alignment score that does not depend on database size. Larger = better."),
             ("identity_pct", "Percent of identical residues in the alignment to the matched ortholog."),
             ("query_cov_pct", "Percent of YOUR sequence covered by the alignment."),
             ("subject_cov_pct", "Percent of the matched eggNOG ortholog covered by the alignment."),
             ("taxonomic_level", "Evolutionary level of the ortholog group used for annotation."),
             ("", ""), ("COUNTS", "")]
    rows += [(k, v) for k, v in t.counts.items()]
    if t.notes:
        rows += [("", ""), ("NOTES", "")] + [("", n) for n in t.notes]
    if enr is not None:
        rows += [("", ""), ("ENRICHMENT", ""),
                 ("IDs given / matched", f"{enr.study_requested} / {enr.study_matched}"),
                 ("Test", "One-sided hypergeometric test, Benjamini-Hochberg FDR per GO aspect"),
                 ("FDR cut-off", enr.fdr_cutoff), ("Minimum count", enr.min_count)]
    return rows


def export_excel(t: ResultTables, enr: Optional[EnrichmentResult] = None) -> bytes:
    """Multi-sheet Excel workbook as bytes."""
    sheets = [("Summary", compact(t.summary, COMPACT_SUMMARY))]
    gl = t.go_long
    glc = compact(gl, COMPACT_GO)
    sheets.append(("GO_all", glc))
    for sheet, aspect in (("GO_MF", "Molecular Function"), ("GO_BP", "Biological Process"),
                          ("GO_CC", "Cellular Component")):
        sheets.append((sheet, glc[glc["aspect"] == aspect] if len(glc) else glc))
    sheets += [("GO_dictionary", t.go_dict),
               ("KEGG_pathways", compact(t.pathway_long, COMPACT_PATHWAY)),
               ("KEGG_pathway_dictionary", t.pathway_dict),
               ("KEGG_KO", compact(t.ko_long, COMPACT_KO)), ("KEGG_KO_dictionary", t.ko_dict),
               ("Summary_full", t.summary)]
    if enr is not None:
        sheets += [("Enrichment_GO", enr.go if len(enr.go) else pd.DataFrame({"note": ["no results"]})),
                   ("Enrichment_KEGG", enr.kegg if len(enr.kegg) else pd.DataFrame({"note": ["no results"]}))]

    buf = io.BytesIO()
    try:
        import xlsxwriter  # noqa: F401
        engine = "xlsxwriter"
    except ImportError:
        engine = "openpyxl"
    with pd.ExcelWriter(buf, engine=engine) as xw:
        readme = pd.DataFrame(_readme_rows(t, enr), columns=["item", "value"])
        readme.to_excel(xw, sheet_name="README", index=False)
        for name, df in sheets:
            df.to_excel(xw, sheet_name=name[:31], index=False)
        if engine == "xlsxwriter":
            wb = xw.book
            head = wb.add_format({"bold": True, "font_color": "white", "bg_color": "#1B6E5A",
                                  "border": 1, "valign": "vcenter"})
            green = wb.add_format({"bg_color": "#C6EFCE", "font_color": "#006100"})
            amber = wb.add_format({"bg_color": "#FFEB9C", "font_color": "#7F6000"})
            red = wb.add_format({"bg_color": "#FFC7CE", "font_color": "#9C0006"})
            sci = wb.add_format({"num_format": "0.00E+00"})
            for name, df in [("README", readme)] + sheets:
                ws = xw.sheets[name[:31]]
                for j, col in enumerate(df.columns):
                    ws.write(0, j, col, head)
                    sample = df[col].head(300).map(lambda v: "" if pd.isna(v) else str(v)) if len(df) else []
                    width = max([len(str(col)) + 2] + [len(x) for x in sample]) if len(df) else len(str(col)) + 2
                    ws.set_column(j, j, min(max(width, 8), 60), sci if col in ("evalue", "p_value", "FDR") else None)
                    if col == "confidence" and len(df):
                        for val, fmt in (("High", green), ("Medium", amber), ("Low", red)):
                            ws.conditional_format(1, j, len(df), j, {"type": "cell", "criteria": "==",
                                                                    "value": f'"{val}"', "format": fmt})
                ws.freeze_panes(1, 0)
                if len(df.columns):
                    ws.autofilter(0, 0, max(len(df), 1), len(df.columns) - 1)
    return buf.getvalue()
