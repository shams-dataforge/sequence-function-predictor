"""
app.py  -  the web interface (Streamlit).   Start it with:   streamlit run app.py
The analysis lives in scripts/pipeline.py and scripts/results.py; wording is in scripts/helptext.py;
names, folders and server limits are in scripts/appconfig.py.
"""
import html
import os
import shutil
import sys
import uuid
from pathlib import Path

import numpy as np
import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))

import appconfig as conf   # noqa: E402
import helptext as ht      # noqa: E402
import pipeline as pl      # noqa: E402
import results as rs       # noqa: E402
from ncbi_fetch import download_sequences_batched  # noqa: E402

st.set_page_config(page_title=conf.APP_NAME, page_icon="🧬", layout="wide",
                   initial_sidebar_state="expanded")

# ============================================================================ look & feel
_inst = f'<span class="chip">{html.escape(conf.INSTITUTE)}</span>' if conf.INSTITUTE else ""
st.markdown("""
<style>
:root {--brand:#1b6e5a; --brand-dark:#0f4c5c; --soft:rgba(27,110,90,.08); --line:rgba(27,110,90,.28);}
.block-container {padding-top: 4.6rem; padding-bottom: 2.5rem; max-width: 1480px;}
.stAppDeployButton, [data-testid="stMainMenu"], [data-testid="stDecoration"], footer {display: none !important;}
.hero {display:flex; flex-wrap:wrap; justify-content:space-between; align-items:center; gap:1rem;
       background: linear-gradient(135deg,#0f4c5c 0%,#1b6e5a 55%,#2a9d8f 100%); border-radius:16px;
       padding:1.5rem 2rem; margin-bottom:1rem; box-shadow:0 6px 18px rgba(15,76,92,.20);}
.hero-title {color:#ffffff; font-size:2rem; font-weight:700; line-height:1.15; letter-spacing:-.01em;}
.hero-sub {color:#e3f4f0; font-size:1.02rem; margin-top:.4rem; max-width:860px; line-height:1.45;}
.hero-side {display:flex; flex-wrap:wrap; gap:.4rem; align-items:center;}
.chip {display:inline-block; background:rgba(255,255,255,.17); color:#ffffff; border:1px solid rgba(255,255,255,.38);
       border-radius:999px; padding:.2rem .8rem; font-size:.8rem; white-space:nowrap;}
.steps {display:grid; grid-template-columns:repeat(auto-fit,minmax(200px,1fr)); gap:.7rem; margin-bottom:1.1rem;}
.step {border:1px solid var(--line); background:var(--soft); border-radius:12px; padding:.65rem .95rem;}
.step b {display:block; font-size:.92rem; margin-bottom:.1rem;}
.step span {font-size:.8rem; opacity:.82;}
[data-testid="stMetric"] {background:var(--soft); border:1px solid var(--line); border-radius:12px;
       padding:.7rem 1rem; box-shadow:0 1px 3px rgba(0,0,0,.05);}
[data-testid="stMetricLabel"] p {font-size:.82rem;}
.stTabs [data-baseweb="tab-list"] {gap:.3rem;}
.stTabs [data-baseweb="tab"] {font-weight:600; padding:.55rem 1.05rem; border-radius:8px 8px 0 0;}
.stButton > button, .stDownloadButton > button {border-radius:10px; font-weight:600;}
section[data-testid="stSidebar"] {border-right:1px solid var(--line);}
section[data-testid="stSidebar"] .block-container {padding-top:2rem;}
.side-brand {background:linear-gradient(135deg,#0f4c5c 0%,#1b6e5a 55%,#2a9d8f 100%); border-radius:14px; padding:.9rem 1rem;
       margin-bottom:1.1rem; box-shadow:0 4px 12px rgba(15,76,92,.22);}
.side-brand b {display:block; color:#fff; font-size:1.05rem; line-height:1.2;}
.side-brand span {color:#e3f4f0; font-size:.74rem;}
.side-title {font-size:.72rem; font-weight:700; letter-spacing:.12em; text-transform:uppercase; opacity:.62; margin:0 0 .5rem .2rem;}
.st-key-menu div[role="radiogroup"] {gap:.5rem; width:100%;}
.st-key-menu div[role="radiogroup"] > label {width:100%; margin:0; padding:.7rem .95rem; border-radius:12px;
       border:1px solid var(--line); background:var(--soft); cursor:pointer; transition:all .15s ease;}
.st-key-menu div[role="radiogroup"] > label:hover {border-color:var(--brand); transform:translateX(3px);}
.st-key-menu div[role="radiogroup"] > label > div:first-child {display:none;}
.st-key-menu div[role="radiogroup"] > label p {font-weight:600; font-size:.95rem;}
.st-key-menu div[role="radiogroup"] > label:has(input:checked) {border-color:transparent;
       background:linear-gradient(135deg,#0f4c5c 0%,#1b6e5a 60%,#2a9d8f 100%); box-shadow:0 4px 12px rgba(15,76,92,.25);}
.st-key-menu div[role="radiogroup"] > label:has(input:checked) p {color:#ffffff;}
.param-card {border:1px solid var(--line); background:var(--soft); border-radius:12px; padding:.7rem .9rem; margin-top:1rem;}
.param-card .ph {font-size:.72rem; font-weight:700; letter-spacing:.12em; text-transform:uppercase; opacity:.62; margin-bottom:.4rem;}
.param-row {padding:.3rem 0; border-top:1px solid var(--line);}
.param-row:first-of-type {border-top:none;}
.param-row span {display:block; font-size:.72rem; opacity:.65;}
.param-row b {font-size:.84rem; font-weight:600; line-height:1.3;}
.param-note {font-size:.76rem; opacity:.8; margin-top:.6rem; line-height:1.4;}
.app-footer {text-align:center; font-size:.78rem; opacity:.72; margin-top:2.6rem; padding-top:1rem;
       border-top:1px solid var(--line); line-height:1.6;}
</style>
""", unsafe_allow_html=True)
st.markdown(f"""<div class="hero"><div><div class="hero-title">🧬 {html.escape(conf.APP_NAME)}</div>
<div class="hero-sub">Predict Molecular Function, Biological Process, Cellular Component and KEGG pathways from DNA or protein sequences, with statistics, confidence levels and full descriptions.</div></div>
<div class="hero-side">{_inst}<span class="chip">version {conf.APP_VERSION}</span></div></div>
<div class="steps">
<div class="step"><b>① Get sequences</b><span>Download from NCBI or use your own FASTA file</span></div>
<div class="step"><b>② Run analysis</b><span>Choose the sequence type, then press Run</span></div>
<div class="step"><b>③ Results</b><span>Statistics, confidence, GO and KEGG with descriptions</span></div>
<div class="step"><b>④ Enrichment</b><span>Optional: find over-represented functions</span></div>
</div>""", unsafe_allow_html=True)

BADGE = {"High": "🟢 High", "Medium": "🟡 Medium", "Low": "🔴 Low",
         "No hit": "⚪ No hit", "Not scored": "⚪ Not scored"}
COLOR_MAP = {"High": "#2e9e6b", "Medium": "#e0b030", "Low": "#d9534f", "No hit": "#b5bdbb"}
STAT_CFG = {
    "evalue": st.column_config.NumberColumn("E-value", format="%.2e", help="Smaller = more reliable match"),
    "bitscore": st.column_config.NumberColumn("Bit score", format="%.1f", help="Larger = better alignment"),
    "identity_pct": st.column_config.NumberColumn("Identity %", format="%.1f"),
    "coverage_pct": st.column_config.NumberColumn(
        "Coverage %", format="%.1f",
        help="Smaller of query and subject coverage - the value used for the confidence level"),
    "query_cov_pct": st.column_config.NumberColumn("Query cov %", format="%.1f"),
    "subject_cov_pct": st.column_config.NumberColumn("Subject cov %", format="%.1f"),
    "protein_length_aa": st.column_config.NumberColumn("Length (aa)"),
}
TH_DEFAULTS = dict(th_he=-30, th_hi=60, th_hc=70, th_me=-10, th_mi=40, th_mc=50)
H = ht.HELP
MH = ht.MENU_HELP
M_SIMPLE, M_ADV, M_NCBI = "🧪  1 · Simple prediction", "🔬  2 · Advanced analysis", "🌐  3 · NCBI fetch"


def stat_cfg(df):
    return {k: v for k, v in STAT_CFG.items() if k in df.columns}


def show_df(df, **kw):
    """st.dataframe that works on old and new Streamlit versions."""
    try:
        st.dataframe(df, width="stretch", hide_index=True, **kw)
    except Exception:
        st.dataframe(df, use_container_width=True, hide_index=True, **kw)


def bar(df, x, y, color=None):
    try:
        import altair as alt
        enc = dict(x=alt.X(x, sort=None, title=None), y=alt.Y(y, title=None))
        if color:
            enc["color"] = alt.Color(color, legend=None, scale=alt.Scale(
                domain=list(COLOR_MAP), range=list(COLOR_MAP.values())))
        chart = alt.Chart(df).mark_bar(cornerRadiusTopLeft=4, cornerRadiusTopRight=4).encode(**enc)
        try:
            st.altair_chart(chart, width="stretch")
        except Exception:
            st.altair_chart(chart, use_container_width=True)
    except Exception:
        st.bar_chart(df.set_index(x)[y])


def safe_upload_name(name: str) -> str:
    suffix = Path(name).suffix.lower()
    return pl.slugify(Path(name).stem) + (suffix if suffix in {".fasta", ".fa", ".fna", ".faa", ".txt"} else ".fasta")


def render_footer():
    _contact = f" · Contact: {html.escape(conf.CONTACT_EMAIL)}" if conf.CONTACT_EMAIL else ""
    _inst_f = f"{html.escape(conf.INSTITUTE)} · " if conf.INSTITUTE else ""
    st.markdown(f"""<div class="app-footer">{_inst_f}{html.escape(conf.APP_NAME)} v{conf.APP_VERSION}{_contact}<br>
Built on eggNOG-mapper, DIAMOND, Prodigal, the Gene Ontology and KEGG · Predictions are computational and should be validated experimentally</div>""",
                unsafe_allow_html=True)


# ============================================================================ cached loaders
@st.cache_resource(show_spinner="Loading the Gene Ontology (first time only)...")
def _go_dag():
    path = ROOT / "go-basic.obo"
    rs.ensure_go_obo(path)
    return rs.load_go_dag(path)


@st.cache_resource(show_spinner="Loading KEGG names (first time only)...")
def _kegg():
    k = rs.load_kegg_names(conf.CACHE_DIR)
    if not (k["ko"] or k["pathway"] or k["module"]):
        raise RuntimeError("; ".join(k["messages"]) or "KEGG names unavailable")
    return k


def go_dag():
    try:
        return _go_dag(), None
    except Exception as exc:
        return None, f"Gene Ontology file unavailable: {exc}"


def kegg():
    try:
        return _kegg()
    except Exception as exc:
        return {"ko": {}, "pathway": {}, "module": {}, "messages": [f"KEGG names unavailable: {exc}"]}


@st.cache_data(show_spinner="Building result tables...")
def _tables(run_dir: str, th: tuple, stamp: float, show_obsolete: bool):
    return rs.build_tables(run_dir, rs.Thresholds(*th), go_dag()[0], kegg(), show_obsolete)


@st.cache_data(show_spinner=False)
def _inspect(path: str, size: int, mtime: float):
    return pl.inspect_fasta(path)


def load_active_tables(th: rs.Thresholds, show_obsolete: bool = False):
    run = st.session_state.get("active_run")
    if not run or not (Path(run) / "run_info.json").exists():
        return None
    stamp = (Path(run) / "run_info.json").stat().st_mtime
    return _tables(run, th.as_tuple(), stamp, show_obsolete)


# ============================================================================ sidebar
_max_cpu = os.cpu_count() or 4
SETTING_DEFAULTS = dict(
    cpu_n=min(6, _max_cpu), sensmode=pl.SENSMODES[0], evalue=1e-3, min_score=60.0, min_pident=0.0,
    min_qcov=0.0, min_scov=0.0, tax_label=list(pl.TAX_SCOPES)[0], orth=pl.TARGET_ORTHOLOGS[0],
    go_label=list(pl.GO_EVIDENCE)[0], min_orf=100, gcode=1, both=True, show_obsolete=False,
    ncbi_term="Escherichia coli[Organism] AND 500:3000[Sequence Length]", ncbi_db="nucleotide",
    ncbi_max=300, ncbi_name="ncbi_fetched.fasta", ncbi_email="your_email@example.com")
for k, v in {**TH_DEFAULTS, **SETTING_DEFAULTS}.items():
    st.session_state.setdefault(k, v)
st.session_state.setdefault("sid", uuid.uuid4().hex[:10])   # private folder for this visitor's uploads
st.session_state.setdefault("my_runs", [])
# A control that is not drawn in the chosen menu option would lose its value; re-assigning keeps it.
for k in [*TH_DEFAULTS, *SETTING_DEFAULTS, "seq_type", "ncbi_example", "detail", "enr_ids"]:
    if st.session_state.get(k) is not None:
        st.session_state[k] = st.session_state[k]


@st.cache_resource(ttl=3600, show_spinner=False)
def _cleanup(_hour_key: int = 0):
    """Really delete saved runs and uploads older than conf.RETENTION_DAYS (checked at most once an hour)."""
    return pl.cleanup_old_data(conf.RUNS_DIR, conf.UPLOAD_DIR, conf.RETENTION_DAYS)


_cleanup()
up_dir = conf.UPLOAD_DIR / st.session_state["sid"]


def reset_thresholds():
    st.session_state.update(TH_DEFAULTS)


def open_in_advanced(run_dir):
    st.session_state["active_run"] = run_dir
    st.session_state.pop("enrichment", None)
    st.session_state["menu"] = M_ADV


def send_fetched(target):
    """Use the file fetched in option 3 in option 1 or 2 without uploading it again."""
    st.session_state["menu"] = target
    if target == M_SIMPLE:
        st.session_state["simple_input"] = st.session_state.get("fetched_path")
    else:
        st.session_state["source"] = "Use the file fetched from NCBI"
        if st.session_state.get("fetched_db") == "protein":
            st.session_state["seq_type"] = pl.INPUT_TYPE_LABELS[pl.INPUT_PROTEIN]


TH_SIMPLE = rs.Thresholds(float(f"1e{TH_DEFAULTS['th_he']}"), float(TH_DEFAULTS["th_hi"]),
                          float(TH_DEFAULTS["th_hc"]), float(f"1e{TH_DEFAULTS['th_me']}"),
                          float(TH_DEFAULTS["th_mi"]), float(TH_DEFAULTS["th_mc"]))

with st.sidebar:
    st.markdown(f'<div class="side-brand"><b>🧬 {html.escape(conf.APP_NAME)}</b><span>version {conf.APP_VERSION}</span></div>'
                '<div class="side-title">Menu</div>', unsafe_allow_html=True)
    menu = st.radio("Menu", [M_SIMPLE, M_ADV, M_NCBI], key="menu", label_visibility="collapsed", help=MH["menu"])

# ---------------------------------------------------------------------------- OPTION 1: SIMPLE PREDICTION
if menu == M_SIMPLE:
    ss = pl.RunSettings(input_type=pl.INPUT_PROTEIN, cpu=min(6, os.cpu_count() or 4))
    _tax = next((k for k, v in pl.TAX_SCOPES.items() if v == ss.tax_scope), "Auto")
    _go = next((k for k, v in pl.GO_EVIDENCE.items() if v == ss.go_evidence), ss.go_evidence)
    _off = lambda x: "off" if not x else f"{x:g}"
    _params = [
        ("Sequence type", pl.INPUT_TYPE_LABELS[pl.INPUT_PROTEIN]),
        ("Search sensitivity", ss.sensmode),
        ("E-value cut-off", f"{ss.evalue:g}"),
        ("Minimum bit score", f"{ss.min_score:g}"),
        ("Min. identity % / query cov. % / subject cov. %",
         f"{_off(ss.min_pident)} / {_off(ss.min_query_cover)} / {_off(ss.min_subject_cover)}"),
        ("Taxonomic scope", _tax),
        ("Ortholog types used", ss.target_orthologs),
        ("GO evidence", _go),
        ("High confidence", f"E-value ≤ {TH_SIMPLE.high_evalue:g}, identity ≥ {TH_SIMPLE.high_identity:g}%, "
                            f"coverage ≥ {TH_SIMPLE.high_coverage:g}%"),
        ("Medium confidence", f"E-value ≤ {TH_SIMPLE.med_evalue:g}, identity ≥ {TH_SIMPLE.med_identity:g}%, "
                              f"coverage ≥ {TH_SIMPLE.med_coverage:g}%"),
    ]
    if conf.RETENTION_DAYS > 0:
        _policy = (f"🗑 Results are deleted from the server automatically {conf.RETENTION_DAYS} "
                   f"day{'s' if conf.RETENTION_DAYS != 1 else ''} after the run. Download your files before then.")
    else:
        _policy = "Results stay on the server until the administrator removes them."
    with st.sidebar:
        st.markdown(
            '<div class="param-card"><div class="ph">Default parameters</div>'
            + "".join(f'<div class="param-row"><span>{html.escape(k)}</span><b>{html.escape(str(v))}</b></div>'
                      for k, v in _params)
            + f'<div class="param-note">{html.escape(_policy)}<br>⏱ The first run after the server starts is the '
              'slowest (about 20 minutes). Keep the page open while it runs.</div></div>',
            unsafe_allow_html=True)
        st.divider()
        st.markdown("### 📁 Saved runs")
        _runs = pl.list_runs(conf.RUNS_DIR)
        if not conf.SHARED_RUNS:
            _runs = [r for r in _runs if r["path"] in st.session_state["my_runs"]]
        if _runs:
            _labels = {f"{r['name']}  ·  {r['created']}": r["path"] for r in _runs}
            _pick = st.selectbox("Open an earlier run", list(_labels), index=None, placeholder="Choose a run",
                                 help="Runs you have made. Results are rebuilt from the saved files, so nothing is re-run.")
            if st.button("Open run", disabled=_pick is None):
                st.session_state["simple_run"] = st.session_state["active_run"] = _labels[_pick]
                st.session_state.pop("enrichment", None)
        else:
            st.caption("No saved runs yet." if conf.SHARED_RUNS else "Runs you start appear here.")
        _rid = st.text_input("Open a run by its ID", key="simple_open_id", placeholder="paste a Run ID",
                             help="Every run has a private Run ID (shown with the results). "
                                  "Keep it to reopen your results later.")
        if st.button("Open by ID", disabled=not _rid.strip()):
            _found = pl.find_run(_rid, conf.RUNS_DIR)
            if _found is None:
                st.error("No finished run with that ID.")
            else:
                st.session_state["simple_run"] = st.session_state["active_run"] = str(_found)
                st.session_state.pop("enrichment", None)

    input_path, default_name, ready = None, "simple_run", False
    with st.container(border=True):
        up = st.file_uploader("Protein FASTA file", type=["fasta", "fa", "faa", "txt"], key="simple_up",
                              help=MH["simple_upload"])
        if up is not None:
            up_dir.mkdir(parents=True, exist_ok=True)
            input_path = up_dir / safe_upload_name(up.name)
            if not input_path.exists() or input_path.stat().st_size != up.size:
                input_path.write_bytes(up.getbuffer())
            default_name = Path(up.name).stem
        else:
            sent = st.session_state.get("simple_input")
            if sent and Path(sent).exists():
                input_path = Path(sent)
                default_name = input_path.stem
                st.info(f"Using the protein file fetched from NCBI: **{input_path.name}**. "
                        "Upload a file here to use that one instead.")
        if input_path is not None:
            _st = Path(input_path).stat()
            fi = _inspect(str(input_path), _st.st_size, _st.st_mtime)
            for w in fi["warnings"]:
                st.warning(w)
            if fi["kind"] == "empty":
                pass
            elif fi["n_sequences"] > conf.MAX_SEQUENCES:
                st.error(f"This file has {fi['n_sequences']:,} sequences. The limit per run on this server is "
                         f"{conf.MAX_SEQUENCES:,}. Please split the file into smaller parts.")
            elif fi["kind"] == "dna":
                st.error("This file looks like **DNA**. Simple prediction accepts protein sequences only - "
                         "use option 2, Advanced analysis, for DNA.")
            else:
                ready = True
                st.caption(f"{fi['n_sequences']:,} protein sequences found.")
        env_problems = pl.check_environment(pl.INPUT_PROTEIN, conf.DB_DIR)
        for p_ in env_problems:
            st.error(p_)
        run_clicked = st.button("🚀 Run", type="primary", disabled=not ready or bool(env_problems),
                                help=MH["simple_run"])

    if run_clicked:
        lines: list = []
        with st.status("Running analysis...", expanded=True) as status:
            box = st.empty()

            def log(msg):
                lines.append(msg)
                box.code("\n".join(lines[-14:]), language="text")
            try:
                run_dir = pl.run_full_pipeline(input_path, default_name.strip() or "simple_run", ss,
                                               data_dir=conf.DB_DIR, runs_dir=conf.RUNS_DIR, log=log)
                status.update(label="Finished", state="complete")
                st.session_state["simple_run"] = str(run_dir)
                st.session_state["active_run"] = str(run_dir)
                st.session_state["my_runs"].append(str(run_dir))
                st.session_state.pop("enrichment", None)
            except pl.PipelineError as exc:
                status.update(label="Could not finish", state="error")
                st.error(str(exc))
            except Exception as exc:  # unexpected
                status.update(label="Could not finish", state="error")
                st.error(f"Unexpected error: {exc}")

    sr = st.session_state.get("simple_run")
    if sr and (Path(sr) / "run_info.json").exists():
        tb = _tables(sr, TH_SIMPLE.as_tuple(), (Path(sr) / "run_info.json").stat().st_mtime, False)
        nm = tb.info.get("run_name", "run")
        _df = tb.summary
        _txt = lambda col: (_df[col].fillna("").astype(str).str.strip() if col in _df.columns
                            else pd.Series("", index=_df.index))
        _go_parts = [_txt(c) for c in ("GO_molecular_function", "GO_biological_process", "GO_cellular_component")]
        _simple = pd.DataFrame({
            "sequence_id": _df["sequence_id"],
            "GO terms": pd.Series(["; ".join(x for x in vals if x) for vals in zip(*_go_parts)], index=_df.index),
            "KEGG pathways": _txt("KEGG_pathway"),
            "confidence": _df["confidence"].map(lambda v: BADGE.get(v, v)),
            "evalue": pd.to_numeric(_df["evalue"], errors="coerce") if "evalue" in _df.columns else None,
            "identity_pct": pd.to_numeric(_df["identity_pct"], errors="coerce") if "identity_pct" in _df.columns else None,
            "coverage_pct": pd.to_numeric(_df["coverage_pct"], errors="coerce") if "coverage_pct" in _df.columns else None,
        })
        with st.container(border=True):
            st.markdown(f"#### Results: {nm}")
            st.caption(f"Run ID: `{Path(sr).name}` · {len(_simple):,} sequences")
            st.button("🔍 Show detailed analysis", type="primary", on_click=open_in_advanced, args=(sr,),
                      help=MH["simple_detail"])
            show_df(_simple, height=480, column_config=stat_cfg(_simple))
            d1, d2, d3 = st.columns(3)
            d1.download_button("⬇️ GO annotations (CSV)", rs.to_csv_bytes(rs.compact(tb.go_long, rs.COMPACT_GO)),
                               f"{nm}_GO.csv", "text/csv", help=MH["simple_dl_go"])
            d2.download_button("⬇️ KEGG pathways (CSV)",
                               rs.to_csv_bytes(rs.compact(tb.pathway_long, rs.COMPACT_PATHWAY)),
                               f"{nm}_KEGG.csv", "text/csv", help=MH["simple_dl_kegg"])
            d3.download_button("⬇️ Full summary (CSV)", rs.to_csv_bytes(tb.summary),
                               f"{nm}_summary_full.csv", "text/csv", help=MH["simple_dl_summary"])
    render_footer()
    st.stop()

# ---------------------------------------------------------------------------- OPTION 3: NCBI FETCH
if menu == M_NCBI:
    with st.sidebar:
        st.divider()
        st.markdown("### ⚙️ NCBI settings")
        db = st.radio("Database", ["nucleotide", "protein"], horizontal=True, key="ncbi_db", help=H["ncbi_db"])
        retmax = st.number_input("Max sequences to fetch", 1, 5000, key="ncbi_max", help=H["ncbi_max"])
        fname = st.text_input("Save as", key="ncbi_name", help=H["ncbi_name"])
        email = st.text_input("Your e-mail (NCBI asks for it)", key="ncbi_email", help=H["ncbi_email"])

    st.markdown("#### Download sequences from NCBI")
    st.caption("Fetch sequences, download them as a FASTA file, or send them straight to Simple prediction or "
               "Advanced analysis.")

    def load_example():
        ex = st.session_state.get("ncbi_example")
        if ex in ht.NCBI_EXAMPLES:
            st.session_state["ncbi_term"], st.session_state["ncbi_db"], _ = ht.NCBI_EXAMPLES[ex]

    with st.container(border=True):
        st.selectbox("Start from an example search", list(ht.NCBI_EXAMPLES), index=None, key="ncbi_example",
                     placeholder="Optional - choose an example", on_change=load_example,
                     help="Fills in the search term below. Remember to choose the matching sequence type in Advanced analysis.")
        term = st.text_input("NCBI search term", key="ncbi_term", help=H["ncbi_term"])
        ex = st.session_state.get("ncbi_example")
        if ex in ht.NCBI_EXAMPLES:
            st.info(f"For this example choose the sequence type: **{ht.NCBI_EXAMPLES[ex][2]}**")
        if st.button("Fetch sequences", type="primary", help="Downloads the sequences and keeps them for the next step."):
            up_dir.mkdir(parents=True, exist_ok=True)
            out = up_dir / safe_upload_name(fname)
            with st.spinner("Fetching from NCBI..."):
                try:
                    n = download_sequences_batched(term, str(out), database=db, retmax=int(retmax), email=email)
                except Exception as exc:
                    n = -1
                    st.error(f"NCBI request failed: {exc}")
            if n == 0:
                st.warning("No sequences found - try a simpler search term.")
            elif n > 0:
                st.session_state["fetched_path"] = str(out)
                st.session_state["fetched_db"] = db
                st.session_state["fetched_n"] = n

    fp = st.session_state.get("fetched_path")
    if fp and Path(fp).exists():
        fdb = st.session_state.get("fetched_db")
        with st.container(border=True):
            st.success(f"Fetched {st.session_state.get('fetched_n', '')} {fdb} sequences and saved them as "
                       f"**{Path(fp).name}**.")
            b1, b2, b3 = st.columns(3)
            b1.download_button("⬇️ Download FASTA", Path(fp).read_bytes(), Path(fp).name, "text/plain",
                               help=MH["ncbi_download"])
            b2.button("➡️ Send to Simple prediction", on_click=send_fetched, args=(M_SIMPLE,),
                      disabled=fdb != "protein", help=MH["ncbi_send_simple"])
            b3.button("➡️ Send to Advanced analysis", on_click=send_fetched, args=(M_ADV,),
                      help=MH["ncbi_send_adv"])
            if fdb != "protein":
                st.caption("Simple prediction accepts protein sequences only. Fetch from the 'protein' database "
                           "to use it, or send these DNA sequences to Advanced analysis.")
    render_footer()
    st.stop()

# ---------------------------------------------------------------------------- OPTION 2: ADVANCED ANALYSIS (sidebar settings below)


with st.sidebar:
    st.markdown("### ⚙️ Settings")
    type_label = st.selectbox("Sequence type", list(pl.INPUT_TYPE_LABELS.values()), index=None, key="seq_type",
                              placeholder="Choose the type of your sequences", help=H["sequence_type"])
    input_type = next((k for k, v in pl.INPUT_TYPE_LABELS.items() if v == type_label), None)

    max_cpu = os.cpu_count() or 4
    cpu = st.slider("CPU threads", 1, max_cpu, key="cpu_n", help=H["cpu"]) if max_cpu > 1 else 1

    with st.expander("Search settings (advanced)"):
        st.caption("Defaults are eggNOG-mapper's own standard values. Hover over ? for an explanation.")
        sensmode = st.selectbox("Search sensitivity", pl.SENSMODES, key="sensmode", help=H["sensmode"])
        evalue = st.select_slider("E-value cut-off", [1e-1, 1e-2, 1e-3, 1e-5, 1e-10, 1e-20, 1e-30],
                                  key="evalue", format_func=lambda v: f"{v:g}", help=H["evalue"])
        min_score = st.number_input("Minimum bit score", 0.0, 2000.0, step=5.0, key="min_score", help=H["min_score"])
        min_pident = st.number_input("Minimum identity % (0 = off)", 0.0, 100.0, step=5.0, key="min_pident", help=H["min_pident"])
        min_qcov = st.number_input("Minimum query coverage % (0 = off)", 0.0, 100.0, step=5.0, key="min_qcov", help=H["min_qcov"])
        min_scov = st.number_input("Minimum subject coverage % (0 = off)", 0.0, 100.0, step=5.0, key="min_scov", help=H["min_scov"])
        tax_label = st.selectbox("Taxonomic scope", list(pl.TAX_SCOPES), key="tax_label", help=H["tax_scope"])
        orth = st.selectbox("Ortholog types used", pl.TARGET_ORTHOLOGS, key="orth", help=H["target_orthologs"])
        go_label = st.selectbox("GO evidence", list(pl.GO_EVIDENCE), key="go_label", help=H["go_evidence"])
    with st.expander("Gene finding (DNA input)"):
        min_orf = st.number_input("Minimum ORF length, amino acids (transcripts)", 10, 2000, step=10, key="min_orf", help=H["min_orf"])
        gcode = st.number_input("Genetic code / translation table (transcripts)", 1, 33, step=1, key="gcode", help=H["genetic_code"])
        both = st.checkbox("Search both strands (transcripts)", key="both", help=H["both_strands"])

    with st.expander("Confidence thresholds"):
        st.caption("High = passes all three High rules; otherwise Medium = passes all Medium rules; otherwise Low. "
                   "Starting values are common practice, not a universal standard.")
        st.markdown("**High**")
        st.slider("E-value ≤ 10^", -100, -1, key="th_he", help=H["th_he"])
        st.slider("Identity ≥ %", 0, 100, key="th_hi", help=H["th_hi"])
        st.slider("Coverage ≥ %", 0, 100, key="th_hc", help=H["th_hc"])
        st.markdown("**Medium**")
        st.slider("E-value ≤ 10^ ", -100, -1, key="th_me", help=H["th_me"])
        st.slider("Identity ≥ % ", 0, 100, key="th_mi", help=H["th_mi"])
        st.slider("Coverage ≥ % ", 0, 100, key="th_mc", help=H["th_mc"])
        st.button("Reset to defaults", on_click=reset_thresholds, help="Back to the starting values.")
        if (st.session_state.th_he > st.session_state.th_me or st.session_state.th_hi < st.session_state.th_mi
                or st.session_state.th_hc < st.session_state.th_mc):
            st.warning("High rules should be at least as strict as Medium rules.")

    with st.expander("Display"):
        st.checkbox("Show obsolete GO terms", key="show_obsolete", help=H["show_obsolete"])

    st.divider()
    st.markdown("### 📁 Saved runs")
    runs = pl.list_runs(conf.RUNS_DIR)
    if not conf.SHARED_RUNS:
        runs = [r for r in runs if r["path"] in st.session_state["my_runs"]]
    if runs:
        labels = {f"{r['name']}  ·  {r['created']}": r["path"] for r in runs}
        pick_run = st.selectbox("Open an earlier run", list(labels), index=None, placeholder="Choose a run",
                                help="Runs you have made. Results are rebuilt from the saved files, so nothing is re-run.")
        if st.button("Open run", disabled=pick_run is None):
            st.session_state["active_run"] = labels[pick_run]
            st.session_state.pop("enrichment", None)
    else:
        st.caption("No saved runs yet." if conf.SHARED_RUNS else "Runs you start appear here.")
    rid = st.text_input("Open a run by its ID", key="open_id", placeholder="paste a Run ID",
                        help="Every run has a private Run ID (shown after the run and under Run info). "
                             "Keep it to reopen your results later.")
    if st.button("Open by ID", disabled=not rid.strip()):
        found = pl.find_run(rid, conf.RUNS_DIR)
        if found is None:
            st.error("No finished run with that ID.")
        else:
            st.session_state["active_run"] = str(found)
            st.session_state.pop("enrichment", None)
    if st.session_state.get("active_run"):
        st.caption(f"Showing: **{Path(st.session_state['active_run']).name}**")

    with st.expander("System status"):
        for lbl, ok in [("eggNOG-mapper", bool(shutil.which("emapper.py"))),
                        ("Prodigal", bool(shutil.which("prodigal"))),
                        ("eggNOG database", (conf.DB_DIR / "eggnog.db").exists()
                         and (conf.DB_DIR / "eggnog_proteins.dmnd").exists())]:
            st.markdown(f"{'✅' if ok else '❌'} {lbl}")
        st.caption(f"Version {conf.APP_VERSION} · {os.cpu_count()} CPU cores · "
                   f"limit {conf.MAX_SEQUENCES:,} sequences per run")

th = rs.Thresholds(float(f"1e{st.session_state.th_he}"), float(st.session_state.th_hi),
                   float(st.session_state.th_hc), float(f"1e{st.session_state.th_me}"),
                   float(st.session_state.th_mi), float(st.session_state.th_mc))
settings = pl.RunSettings(
    input_type=input_type or pl.INPUT_PROKARYOTE, cpu=cpu, sensmode=sensmode, evalue=float(evalue),
    min_score=float(min_score), min_pident=float(min_pident), min_query_cover=float(min_qcov),
    min_subject_cover=float(min_scov), tax_scope=pl.TAX_SCOPES[tax_label], target_orthologs=orth,
    go_evidence=pl.GO_EVIDENCE[go_label], min_orf_aa=int(min_orf), genetic_code=int(gcode),
    orf_both_strands=bool(both))

# ============================================================================ tabs
tab_run, tab_res, tab_enr, tab_help = st.tabs(
    ["② Run analysis", "③ Results", "④ Enrichment", "📘 Help"])
show_obs = bool(st.session_state.get("show_obsolete", False))
tables = load_active_tables(th, show_obs)

# ---------------------------------------------------------------------------- ② RUN
with tab_run:
    input_path, default_name = None, "my_run"
    with st.container(border=True):
        st.markdown("##### 1 · Choose your sequences")
        options = ["Upload a FASTA file", "Paste sequences"]
        if st.session_state.get("fetched_path") and Path(st.session_state["fetched_path"]).exists():
            options.append("Use the file fetched from NCBI")
        if st.session_state.get("source") not in options:
            st.session_state.pop("source", None)
        source = st.radio("Input", options, horizontal=True, label_visibility="collapsed", key="source",
                          help=H["source"])
        if source == "Upload a FASTA file":
            up = st.file_uploader("FASTA file", type=["fasta", "fa", "fna", "faa", "txt"], help=H["upload"])
            if up is not None:
                up_dir.mkdir(parents=True, exist_ok=True)
                input_path = up_dir / safe_upload_name(up.name)
                if not input_path.exists() or input_path.stat().st_size != up.size:
                    input_path.write_bytes(up.getbuffer())
                default_name = Path(up.name).stem
        elif source == "Paste sequences":
            pasted = st.text_area("FASTA text", height=170, placeholder=">seq1\nATGGCC...\n>seq2\nATGAAA...",
                                  help=H["paste"])
            if pasted.strip():
                up_dir.mkdir(parents=True, exist_ok=True)
                input_path = up_dir / "pasted_sequences.fasta"
                input_path.write_text(pasted if pasted.lstrip().startswith(">") else ">sequence1\n" + pasted)
                default_name = "pasted"
        else:
            input_path = Path(st.session_state["fetched_path"])
            default_name = input_path.stem

    if input_path is None:
        st.info("Upload or paste sequences to begin (or fetch them in option 3). New here? Open the 📘 Help tab for a 1-minute guide.")
    else:
        ready = False
        with st.container(border=True):
            st.markdown("##### 2 · Check")
            _st = Path(input_path).stat()
            fi = _inspect(str(input_path), _st.st_size, _st.st_mtime)
            m = st.columns(5)
            m[0].metric("Sequences", f"{fi['n_sequences']:,}")
            m[1].metric("Looks like", {"dna": "DNA", "protein": "Protein"}.get(fi["kind"], "-"))
            m[2].metric("Shortest", f"{fi['min_len']:,}")
            m[3].metric("Median length", f"{fi['median_len']:,}")
            m[4].metric("Longest", f"{fi['max_len']:,}")
            for w in fi["warnings"]:
                st.warning(w)
            if input_type is None:
                st.info("👈 Choose the **sequence type** in the sidebar to continue.")
            elif fi["kind"] == "empty":
                pass
            elif fi["n_sequences"] > conf.MAX_SEQUENCES:
                st.error(f"This file has {fi['n_sequences']:,} sequences. The limit per run on this server is "
                         f"{conf.MAX_SEQUENCES:,}. Please split the file into smaller parts.")
            elif fi["kind"] == "protein" and input_type != pl.INPUT_PROTEIN:
                st.error("This file looks like **protein** sequences, but a DNA type is selected. "
                         "Change the sequence type in the sidebar.")
            elif fi["kind"] == "dna" and input_type == pl.INPUT_PROTEIN:
                st.error("This file looks like **DNA**, but 'Protein sequences' is selected. "
                         "Change the sequence type in the sidebar.")
            else:
                ready = True
                if input_type == pl.INPUT_TRANSCRIPT and fi["max_len"] > 50_000:
                    st.warning("Some sequences are very long (genome-like). For transcripts only the longest ORF "
                               "of each sequence is used. For bacterial or viral genomes choose the bacterial / "
                               "archaeal / viral option instead.")
        with st.container(border=True):
            st.markdown("##### 3 · Run")
            run_name = st.text_input("Name for this run", default_name, help=H["run_name"])
            with st.expander("eggNOG database location"):
                db_dir = st.text_input("Folder containing eggnog.db and eggnog_proteins.dmnd",
                                       str(conf.DB_DIR), help=H["db_dir"])
            env_problems = pl.check_environment(input_type or pl.INPUT_PROKARYOTE, db_dir)
            for p_ in env_problems:
                st.error(p_)
            st.info("⏱ The first run after the server starts is the slowest (about 20 minutes) because the eggNOG "
                    "database must be loaded. Keep this page open while the analysis runs.")
            if st.button("🚀 Run analysis", type="primary",
                         disabled=not ready or bool(env_problems) or not run_name.strip(),
                         help="Starts the analysis with the settings in the sidebar."):
                lines: list = []
                with st.status("Running analysis...", expanded=True) as status:
                    box = st.empty()

                    def log(msg):
                        lines.append(msg)
                        box.code("\n".join(lines[-14:]), language="text")
                    try:
                        run_dir = pl.run_full_pipeline(input_path, run_name.strip(), settings,
                                                       data_dir=db_dir, runs_dir=conf.RUNS_DIR, log=log)
                        status.update(label="Finished", state="complete")
                        st.session_state["active_run"] = str(run_dir)
                        st.session_state["my_runs"].append(str(run_dir))
                        st.session_state.pop("enrichment", None)
                        tables = load_active_tables(th, show_obs)
                        st.success("Done! Open the **③ Results** tab.")
                        st.markdown(f"**Your Run ID** (keep it to reopen these results later):")
                        st.code(Path(run_dir).name, language="text")
                    except pl.PipelineError as exc:
                        status.update(label="Could not finish", state="error")
                        st.error(str(exc))
                    except Exception as exc:  # unexpected
                        status.update(label="Could not finish", state="error")
                        st.error(f"Unexpected error: {exc}")

# ---------------------------------------------------------------------------- ③ RESULTS
with tab_res:
    if tables is None:
        st.info("No results to show yet. Run an analysis in step ②, or open a saved run from the sidebar.")
    else:
        c = tables.counts
        run_name_shown = tables.info.get("run_name", "")
        st.markdown(f"#### Results: {run_name_shown}")
        st.caption(f"Run ID: `{Path(st.session_state['active_run']).name}`")
        for n in tables.notes:
            st.warning(n)
        n = max(c["proteins_analysed"], 1)
        mh = ht.METRIC_HELP
        m = st.columns(6)
        m[0].metric("Proteins analysed", f"{c['proteins_analysed']:,}", help=mh["proteins"])
        m[1].metric("With ortholog hit", f"{c['with_ortholog_hit']:,}", f"{100*c['with_ortholog_hit']/n:.0f}%",
                    delta_color="off", help=mh["hit"])
        m[2].metric("Annotated (GO/KEGG)", f"{c['annotated']:,}", f"{100*c['annotated']/n:.0f}%",
                    delta_color="off", help=mh["annotated"])
        m[3].metric("High confidence", f"{c['high']:,}", help=mh["high"])
        m[4].metric("With GO terms", f"{c['with_GO']:,}", help=mh["go"])
        m[5].metric("With KEGG pathway", f"{c['with_pathway']:,}", help=mh["kegg"])

        ch1, ch2 = st.columns(2)
        with ch1:
            st.markdown("**Confidence levels**")
            conf_df = pd.DataFrame({"level": ["High", "Medium", "Low", "No hit"],
                                    "sequences": [c["high"], c["medium"], c["low"], c["no_hit"]]})
            bar(conf_df, "level", "sequences", color="level")
        with ch2:
            st.markdown("**Share of proteins with each annotation (%)**")
            cov = pd.DataFrame({"annotation": ["Ortholog", "GO MF", "GO BP", "GO CC", "KEGG KO", "KEGG pathway"],
                                "percent": [100 * c[k] / n for k in ("with_ortholog_hit", "with_GO_MF", "with_GO_BP",
                                                                      "with_GO_CC", "with_KO", "with_pathway")]})
            bar(cov, "annotation", "percent")

        detail = st.radio("Table columns", ["Compact", "Full"], horizontal=True, key="detail", help=H["detail"])
        compact_view = detail == "Compact"

        def pick(df, cols):
            return (rs.compact(df, cols) if compact_view else df).copy()

        t_seq, t_go, t_kegg, t_info = st.tabs(["Sequences", "GO terms", "KEGG", "Run info & downloads"])

        with t_seq:
            f1, f2, f3 = st.columns([2, 2, 3])
            levels = f1.multiselect("Confidence", ["High", "Medium", "Low", "No hit"],
                                    default=["High", "Medium", "Low", "No hit"], help=H["conf_filter"])
            status_opts = sorted(tables.summary["annotation_status"].unique())
            status_f = f2.multiselect("Status", status_opts, default=status_opts, help=H["status_filter"])
            q = f3.text_input("Search (ID, gene name, description...)", "", help=H["search_box"])
            v = tables.summary[tables.summary["confidence"].isin(levels) &
                               tables.summary["annotation_status"].isin(status_f)]
            if q.strip():
                mask = v.astype(str).apply(lambda col: col.str.contains(q.strip(), case=False, regex=False, na=False)).any(axis=1)
                v = v[mask]
            st.caption(f"Showing {len(v):,} of {len(tables.summary):,} sequences")
            vv = pick(v, rs.COMPACT_SUMMARY)
            vv["confidence"] = vv["confidence"].map(BADGE)
            show_df(vv, height=520, column_config=stat_cfg(vv))

        with t_go:
            asp = st.radio("Category", ["All", "Molecular Function", "Biological Process", "Cellular Component"],
                           horizontal=True, help=H["aspect"])
            gq = st.text_input("Search GO ID or term", "", key="goq", help=H["search_box"])
            gl = tables.go_long
            if asp != "All":
                gl = gl[gl["aspect"] == asp]
            if gq.strip():
                gl = gl[gl["GO_ID"].str.contains(gq.strip(), case=False, regex=False, na=False) |
                        gl["GO_term"].str.contains(gq.strip(), case=False, regex=False, na=False)]
            st.caption(f"{len(gl):,} annotations (one row per sequence and GO term, no duplicates)")
            g2 = pick(gl, rs.COMPACT_GO)
            g2["confidence"] = g2["confidence"].map(BADGE)
            show_df(g2, height=420, column_config=stat_cfg(g2))
            with st.expander("📖 GO term dictionary - every term found in this run, with its description", expanded=True):
                gd = tables.go_dict if asp == "All" else tables.go_dict[tables.go_dict["aspect"] == asp]
                show_df(gd, height=300)
            with st.expander("🔎 Look up any GO ID"):
                gid = st.text_input("GO ID, e.g. GO:0005737", "", key="golook", help=H["go_lookup"]).strip().upper()
                if gid:
                    dag, err = go_dag()
                    if err:
                        st.warning(err)
                    else:
                        i = rs.go_info(dag, gid)
                        if i["aspect"] == "Unknown":
                            st.warning("This ID was not found in the loaded ontology.")
                        else:
                            st.markdown(f"**{gid} - {i['name']}**  \n*{i['aspect']}*")
                            st.write(i["definition"] or "(no definition)")
                            par = rs.go_parents(dag, gid)
                            if par:
                                st.caption("Parent terms: " + "; ".join(par))
                            hit = tables.go_dict[tables.go_dict["GO_ID"] == gid]
                            st.caption(f"Sequences in this run with this term: "
                                       f"{int(hit['n_sequences'].iloc[0]) if len(hit) else 0}")

        with t_kegg:
            if tables.pathway_long.empty and tables.ko_long.empty:
                st.info("No KEGG annotations in this run.")
            else:
                k1, k2 = st.columns([3, 2])
                with k1:
                    st.markdown("**KEGG pathways found** (number of sequences)")
                    show_df(tables.pathway_dict, height=300)
                with k2:
                    if len(tables.pathway_dict):
                        top = tables.pathway_dict.head(12).copy()
                        top["label"] = top["pathway_name"].where(top["pathway_name"] != rs.NAME_NA, top["pathway_ID"])
                        bar(top, "label", "n_sequences")
                with st.expander("Pathway per sequence", expanded=False):
                    kp = pick(tables.pathway_long, rs.COMPACT_PATHWAY)
                    kp["confidence"] = kp["confidence"].map(BADGE)
                    show_df(kp, height=300, column_config=stat_cfg(kp))
                with st.expander("KEGG orthology (KO) per sequence", expanded=False):
                    ko = pick(tables.ko_long, rs.COMPACT_KO)
                    ko["confidence"] = ko["confidence"].map(BADGE)
                    show_df(ko, height=300, column_config=stat_cfg(ko))
                st.caption("KEGG names come from KEGG's web service (free for academic use).")

        with t_info:
            i1, i2 = st.columns(2)
            with i1:
                st.markdown("**Run details**")
                info = tables.info
                s = info.get("settings", {})
                rows = [("Run ID", info.get("run_id")), ("Run name", info.get("run_name")),
                        ("Date", info.get("created")), ("Input file", info.get("input_file")),
                        ("Sequence type", info.get("input_type_label")),
                        ("Software version", info.get("app_version", "")),
                        ("eggNOG-mapper", info.get("emapper_version", "")),
                        ("Prodigal", info.get("prodigal_version", "")), ("Run time (s)", info.get("seconds", ""))]
                rows += [(k, v_) for k, v_ in s.items()]
                show_df(pd.DataFrame(rows, columns=["item", "value"]).astype(str))
            with i2:
                st.markdown("**How to read the statistics**")
                st.markdown(
                    "- **E-value** - matches this good expected by chance. Smaller is better.\n"
                    "- **Bit score** - alignment score, comparable between searches. Larger is better.\n"
                    "- **Identity %** - share of identical residues in the alignment.\n"
                    "- **Coverage %** - the smaller of query coverage (share of *your* sequence aligned) and "
                    "subject coverage (share of the matched ortholog aligned).\n"
                    f"- **High** = E-value ≤ {th.high_evalue:g}, identity ≥ {th.high_identity:g}%, coverage ≥ {th.high_coverage:g}%\n"
                    f"- **Medium** = E-value ≤ {th.med_evalue:g}, identity ≥ {th.med_identity:g}%, coverage ≥ {th.med_coverage:g}%\n"
                    "- **Low** = passed the search cut-offs but not the Medium rule.")
            st.markdown("**Downloads**")
            d1, d2, d3 = st.columns(3)
            tag = detail.lower()
            d1.download_button(f"⬇️ Summary (CSV, {tag})", rs.to_csv_bytes(pick(tables.summary, rs.COMPACT_SUMMARY)),
                               f"{run_name_shown}_summary_{tag}.csv", "text/csv")
            d2.download_button(f"⬇️ GO annotations (CSV, {tag})", rs.to_csv_bytes(pick(tables.go_long, rs.COMPACT_GO)),
                               f"{run_name_shown}_GO_{tag}.csv", "text/csv")
            d3.download_button(f"⬇️ KEGG pathways (CSV, {tag})", rs.to_csv_bytes(pick(tables.pathway_long, rs.COMPACT_PATHWAY)),
                               f"{run_name_shown}_KEGG_{tag}.csv", "text/csv")
            st.caption("The Excel workbook has compact sheets plus a 'Summary_full' sheet with every column.")
            if st.button("📦 Prepare the full Excel workbook", help="Builds one workbook with all tables and the run settings."):
                with st.spinner("Building workbook..."):
                    st.session_state["xlsx"] = (st.session_state["active_run"],
                                                rs.export_excel(tables, st.session_state.get("enrichment")))
            xl = st.session_state.get("xlsx")
            if xl and xl[0] == st.session_state.get("active_run"):
                st.download_button("⬇️ Download Excel workbook (all tables + run info)", xl[1],
                                   f"{run_name_shown}_results.xlsx",
                                   "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", type="primary")
                st.caption("The workbook uses the confidence thresholds set when it was prepared - prepare it again after changing them.")
            st.markdown("**Methods text for your report**")
            dag_, _err = go_dag()
            mtxt = rs.methods_text(tables, getattr(dag_, "data_version", "") if dag_ is not None else "",
                                   st.session_state.get("enrichment"), conf.APP_NAME)
            st.text_area("Copy this into the Methods section (it uses the settings of this run).", mtxt, height=230,
                         key=f"methods_{abs(hash(mtxt))}",
                         help="A ready-written Methods paragraph built from the settings, tool versions and "
                              "confidence thresholds of this run. Check it, then paste it into your report.")
            with st.expander("Raw files of this run"):
                st.code(str(Path(st.session_state["active_run"])), language="text")

# ---------------------------------------------------------------------------- ④ ENRICHMENT
with tab_enr:
    st.markdown("#### Enrichment: are some functions over-represented in a group of sequences?")
    st.caption("Pick a group of sequences you care about (for example, genes that respond to a treatment). The test "
               "asks whether a GO term or KEGG pathway is more common in your group than expected by chance, "
               "compared with ALL annotated sequences of this run. See the 📘 Help tab for a simple explanation.")
    if tables is None:
        st.info("Open or run an analysis first.")
    else:
        def fill_example():
            ann = tables.go_long["sequence_id"].drop_duplicates().head(15)
            st.session_state["enr_ids"] = "\n".join(ann)

        with st.container(border=True):
            st.text_area("Sequence IDs of interest (one per line, or separated by commas / spaces)",
                         key="enr_ids", height=140, help=H["enr_ids"])
            upl = st.file_uploader("...or upload a text file with IDs", type=["txt", "csv", "tsv"],
                                   key="enr_file", help=H["enr_file"])
            e1, e2, e3 = st.columns(3)
            min_count = e1.number_input("Minimum sequences per term", 1, 50, 2, help=H["enr_min"])
            fdr_cut = e2.selectbox("FDR cut-off", [0.05, 0.01, 0.1], index=0, help=H["enr_fdr"])
            e3.button("Fill with an example", on_click=fill_example,
                      help="Uses the first 15 sequences that have GO terms, just to try the test.")
            run_enr = st.button("Run enrichment", type="primary", help="Runs the over-representation test.")
        if run_enr:
            raw = st.session_state.get("enr_ids", "")
            if upl is not None:
                raw += "\n" + upl.getvalue().decode("utf-8", errors="replace")
            ids = [x for x in raw.replace(",", " ").replace(";", " ").split() if x]
            if not ids:
                st.warning("Please give at least one sequence ID.")
            else:
                with st.spinner("Testing..."):
                    st.session_state["enrichment"] = rs.run_enrichment(tables, go_dag()[0], ids,
                                                                       int(min_count), float(fdr_cut))
        er = st.session_state.get("enrichment")
        if er is not None:
            st.markdown(f"**{er.study_matched}** of {er.study_requested} IDs matched sequences in this run.")
            if er.not_found:
                st.caption("Not found: " + ", ".join(er.not_found[:15]) + (" ..." if len(er.not_found) > 15 else ""))
            for n_ in er.notes:
                st.warning(n_)
            sig_go = int(er.go["significant"].sum()) if len(er.go) else 0
            sig_k = int(er.kegg["significant"].sum()) if len(er.kegg) else 0
            st.success(f"Significant at FDR ≤ {er.fdr_cutoff}: {sig_go} GO terms, {sig_k} KEGG pathways.")
            tg, tk = st.tabs(["GO terms", "KEGG pathways"])
            fmt = {"p_value": st.column_config.NumberColumn("p-value", format="%.2e"),
                   "FDR": st.column_config.NumberColumn("FDR", format="%.2e"),
                   "fold_enrichment": st.column_config.NumberColumn("Fold enrichment", format="%.2f")}
            with tg:
                if len(er.go):
                    show_df(er.go, height=360, column_config=fmt)
                    top = er.go[er.go["significant"]].head(12).copy()
                    if len(top):
                        top["neg_log10_FDR"] = -np.log10(top["FDR"].clip(lower=1e-300))
                        st.markdown("**Most significant GO terms** (higher bar = stronger evidence)")
                        bar(top, "GO_term", "neg_log10_FDR")
                    st.download_button("⬇️ GO enrichment (CSV)", rs.to_csv_bytes(er.go), "enrichment_GO.csv", "text/csv")
                else:
                    st.info("No GO terms reached the minimum count.")
            with tk:
                if len(er.kegg):
                    show_df(er.kegg, height=360, column_config=fmt)
                    st.download_button("⬇️ KEGG enrichment (CSV)", rs.to_csv_bytes(er.kegg), "enrichment_KEGG.csv", "text/csv")
                else:
                    st.info("No KEGG pathways reached the minimum count.")
            st.caption("Method: one-sided hypergeometric test (equivalent to Fisher's exact test), annotations "
                       "propagated to parent GO terms, Benjamini-Hochberg FDR within each GO category.")

# ---------------------------------------------------------------------------- HELP
with tab_help:
    st.markdown("#### How to use this tool")
    st.markdown(ht.QUICK_START)
    with st.expander("Which sequence type should I choose?", expanded=True):
        show_df(ht.SEQUENCE_TYPES)
        st.caption("Not supported directly: raw genomic DNA of plants, animals and fungi (genes with introns). "
                   "Upload predicted proteins or transcripts instead.")
    with st.expander("Every setting explained"):
        st.caption("The same explanations appear as ? hover tips next to each control.")
        show_df(ht.glossary_df(), height=420)
    with st.expander("What do the result columns mean?"):
        show_df(ht.columns_df(), height=420)
    with st.expander("What do the confidence levels mean?"):
        st.markdown(ht.CONFIDENCE_TEXT)
    with st.expander("What is enrichment?"):
        st.markdown(ht.ENRICHMENT_TEXT)
    with st.expander("Limitations - please read"):
        st.markdown(ht.LIMITATIONS_TEXT)
    with st.expander("Frequently asked questions"):
        for q_, a_ in ht.FAQ:
            st.markdown(f"**{q_}**")
            st.write(a_)
    with st.expander("Methods and how to cite"):
        st.markdown("This tool combines the following software and resources. Please cite them in your work "
                    "(check volume and page details against the journal pages before publishing):")
        for ref in ht.CITATIONS:
            st.markdown(f"- {ref}")
        st.caption("The 'Run info & downloads' tab writes a Methods paragraph with the exact settings of each run.")
    with st.expander("About this service"):
        st.markdown(f"**{conf.APP_NAME}** version {conf.APP_VERSION}" + (f" - {conf.INSTITUTE}" if conf.INSTITUTE else ""))
        st.markdown("Predictions are computational and should be validated experimentally. "
                    + (f"Questions or problems: **{conf.CONTACT_EMAIL}**" if conf.CONTACT_EMAIL else
                       "Contact the administrator of this service for questions or problems."))

render_footer()