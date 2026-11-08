"""
helptext.py  -  every explanation shown in the app (hover tips, glossary, FAQ, citations).
Edit the wording here; nothing else needs to change.
"""
import pandas as pd

import appconfig

_DAYS = appconfig.RETENTION_DAYS
_KEPT = (f"Saved runs and uploaded files are deleted from the server automatically {_DAYS} day"
         f"{'s' if _DAYS != 1 else ''} after the run." if _DAYS > 0
         else "Saved runs stay on the server until the administrator removes them.")

# (key, label, default, explanation).  The explanation is used BOTH as the "?" hover tip and in the Help tab.
SETTINGS = [
    ("sequence_type", "Sequence type", "you choose",
     "How your sequences are turned into proteins. Transcripts / mRNA / coding sequences: the longest open "
     "reading frame of each sequence is used (any organism). Bacterial, archaeal or viral DNA: genes are predicted "
     "with Prodigal. Protein: used as they are. A wrong choice gives meaningless results, so you must pick one."),
    ("cpu", "CPU threads", "6 (or all cores if fewer)",
     "Number of processor cores used for the search. More cores finish sooner; leave one or two free if other "
     "people use this computer."),
    ("sensmode", "Search sensitivity", "fast",
     "How thoroughly DIAMOND searches. 'fast' is quickest; 'sensitive' and above find more distant relatives but "
     "take longer. Use a more sensitive mode for final results if time allows."),
    ("evalue", "E-value cut-off", "0.001",
     "Matches with an E-value above this are ignored. The E-value is the number of matches this good you would "
     "expect by chance, so smaller is stricter. 0.001 is eggNOG-mapper's standard."),
    ("min_score", "Minimum bit score", "60",
     "Matches with a bit score below this are ignored. Larger scores mean better alignments. 60 is "
     "eggNOG-mapper's standard."),
    ("min_pident", "Minimum identity %", "0 (off)",
     "Ignore matches whose percent identity (identical residues in the alignment) is lower than this. "
     "0 switches the filter off, which is the standard."),
    ("min_qcov", "Minimum query coverage %", "0 (off)",
     "Ignore matches that cover less than this percentage of YOUR sequence. 0 = off (standard)."),
    ("min_scov", "Minimum subject coverage %", "0 (off)",
     "Ignore matches that cover less than this percentage of the matched eggNOG sequence. 0 = off (standard)."),
    ("tax_scope", "Taxonomic scope", "Auto",
     "Restricts which evolutionary groups annotations may come from. 'Auto' lets eggNOG pick the best level for "
     "each sequence and works for any organism. Narrow it (for example to Plants) only if you want annotations "
     "from that group alone."),
    ("target_orthologs", "Ortholog types used", "all",
     "Which kinds of ortholog relationship may be used to transfer function. 'all' is the standard; 'one2one' is "
     "the strictest and gives fewer but more reliable annotations."),
    ("go_evidence", "GO evidence", "Non-electronic",
     "Which GO terms may be transferred. Non-electronic (standard) uses terms backed by non-electronic evidence; "
     "Experimental only is strictest; All also includes electronically inferred (IEA) terms: more GO terms, "
     "lower reliability."),
    ("min_orf", "Minimum ORF length (aa)", "100",
     "Transcript input only. Open reading frames shorter than this many amino acids are ignored. Lower it for "
     "short genes."),
    ("genetic_code", "Genetic code", "1",
     "Transcript input only. NCBI translation table: 1 = standard (most organisms), 11 = bacterial, "
     "2 = vertebrate mitochondrial."),
    ("both_strands", "Search both strands", "on",
     "Transcript input only. Keep on if the orientation of your sequences is unknown (assembled transcripts). "
     "Turn off if all are mRNA / CDS in the forward (sense) orientation."),
    ("th_he", "High: E-value", "1e-30",
     "A sequence is High confidence if it passes all three High rules. Here: its E-value must be at most 10 to "
     "this power (-30 means 1e-30). More negative = stricter."),
    ("th_hi", "High: identity %", "60",
     "High confidence needs at least this percent identity to the matched ortholog."),
    ("th_hc", "High: coverage %", "70",
     "High confidence needs at least this coverage. Coverage is the smaller of query and subject coverage."),
    ("th_me", "Medium: E-value", "1e-10",
     "A sequence is Medium if it is not High but passes all three Medium rules. Here: E-value at most 10 to "
     "this power."),
    ("th_mi", "Medium: identity %", "40", "Medium confidence needs at least this percent identity."),
    ("th_mc", "Medium: coverage %", "50", "Medium confidence needs at least this coverage (smaller of query and subject)."),
    ("show_obsolete", "Show obsolete GO terms", "off",
     "eggNOG's database is older than the current Gene Ontology, so some of its GO terms are now obsolete. They are "
     "hidden by default because obsolete terms are not valid annotations."),
    ("detail", "Table columns", "Compact",
     "Compact shows the columns most people need; Full shows every column. You can also hide or show single "
     "columns with the eye icon above any table."),
    ("run_name", "Name for this run", "file name",
     "A label for this run. It names the saved results so you can find them later."),
    ("db_dir", "eggNOG database folder", "set by the administrator",
     "Folder on the server that holds eggnog.db and eggnog_proteins.dmnd. Normally you do not need to change it."),
    ("ncbi_term", "NCBI search term", "-",
     "The same kind of query you would type on the NCBI website, for example  Mus musculus[Organism] AND "
     "biomol_mrna[PROP]. If a search finds 0 sequences, simplify it."),
    ("ncbi_db", "NCBI database", "nucleotide",
     "'nucleotide' downloads DNA / RNA sequences. 'protein' downloads protein sequences (then choose 'Protein "
     "sequences' as the sequence type)."),
    ("ncbi_max", "Max sequences to fetch", "300", "Upper limit on how many sequences are downloaded."),
    ("ncbi_name", "Save as", "ncbi_fetched.fasta", "File name for the downloaded sequences."),
    ("ncbi_email", "Your e-mail", "-",
     "NCBI asks every automated user to identify themselves with an e-mail address. It is only sent to NCBI."),
    ("conf_filter", "Confidence filter", "all",
     "Show only sequences with the chosen confidence levels."),
    ("status_filter", "Status filter", "all",
     "Annotated = has GO or KEGG information. 'Ortholog hit, no GO/KEGG' = a match was found but it carries no "
     "GO/KEGG data. 'No ortholog hit' = nothing was found with the chosen cut-offs."),
    ("search_box", "Search", "-", "Finds rows containing this text in any column (ID, gene name, description ...)."),
    ("aspect", "GO category", "All",
     "Molecular Function = what the protein does at the molecular level. Biological Process = the larger process "
     "it takes part in. Cellular Component = where in the cell it is found."),
    ("go_lookup", "Look up a GO ID", "-", "Type a GO ID such as GO:0005737 to see its name, definition and parent terms."),
    ("enr_ids", "Sequence IDs of interest", "-",
     "The group of sequences you want to test, for example genes that respond to a treatment. One per line, or "
     "separated by commas or spaces. Use IDs from the 'sequence_id' or 'source_sequence' column."),
    ("enr_file", "ID file", "-", "A plain text file with one ID per line, as an alternative to pasting."),
    ("enr_min", "Minimum sequences per term", "2",
     "Terms found in fewer sequences of your group than this are not tested. Raising it removes weak, one-off hits."),
    ("enr_fdr", "FDR cut-off", "0.05",
     "False discovery rate: the share of 'significant' results you accept to be wrong. Because hundreds of terms are "
     "tested at once, p-values are corrected; 0.05 is the usual choice, 0.01 is stricter."),
    ("source", "Input", "-", "Where the sequences come from: a file on your computer, pasted text, or a file you downloaded from NCBI."),
    ("upload", "FASTA file", "-", "A FASTA file: each sequence starts with a '>' line holding its name."),
    ("paste", "FASTA text", "-", "Paste sequences in FASTA format (a '>' name line followed by the sequence)."),
]
HELP = {k: text for k, _label, _default, text in SETTINGS}

# Hover tips for the menu and for the controls of "Simple prediction" and "NCBI fetch".
MENU_HELP = {
    "menu": "1 = Simple prediction: upload a protein FASTA file and run it with the default settings. "
            "2 = Advanced analysis: all settings, results tables and enrichment. "
            "3 = NCBI fetch: download sequences from NCBI.",
    "simple_upload": "A protein FASTA file: each sequence starts with a '>' line holding its name. "
                     "DNA files are not accepted here - use Advanced analysis for DNA.",
    "simple_run": "Starts the prediction with the default settings listed above.",
    "simple_dl_go": "All GO annotations of this run: one row per sequence and GO term.",
    "simple_dl_kegg": "All KEGG pathway assignments of this run: one row per sequence and pathway.",
    "simple_dl_summary": "One row per sequence with every column: confidence, statistics, description, GO, KEGG and more.",
    "simple_detail": "Opens this run in Advanced analysis, where you can change the confidence thresholds, "
                     "browse the tables and run enrichment.",
    "ncbi_download": "Saves the fetched sequences as a FASTA file on your computer.",
    "ncbi_send_simple": "Uses the fetched file in Simple prediction, so you do not have to upload it again. "
                        "Only available for protein sequences.",
    "ncbi_send_adv": "Uses the fetched file in Advanced analysis, so you do not have to upload it again.",
}

METRIC_HELP = {
    "proteins": "Number of protein sequences analysed (one per input sequence for transcripts and proteins; "
                "possibly several per contig for bacterial DNA).",
    "hit": "Proteins for which eggNOG found a matching ortholog with the chosen cut-offs.",
    "annotated": "Proteins that received GO and / or KEGG information.",
    "high": "Proteins whose best match passes all High-confidence rules (see the sidebar thresholds).",
    "go": "Proteins with at least one valid GO term.",
    "kegg": "Proteins assigned to at least one KEGG pathway.",
}

QUICK_START = """
1. **Get sequences** - download them from NCBI (menu option 3, NCBI fetch) or skip this and upload your own FASTA file.
2. **Choose the sequence type** in the sidebar (DNA transcripts / bacterial DNA / protein).
3. **Run analysis** (tab ②): upload or paste the sequences, give the run a name and press *Run analysis*.
   Keep the page open while it runs. The first run after the server starts is the slowest.
4. **Read the results** (tab ③): sequences with statistics and confidence, GO terms with descriptions, KEGG pathways,
   and downloads (CSV or one Excel workbook).
5. **Optional - enrichment** (tab ④): test whether some functions are over-represented in a group of sequences.
6. **Keep your Run ID.** It lets you reopen the results later from the sidebar.
"""

SEQUENCE_TYPES = pd.DataFrame([
    ("mRNA, cDNA, transcripts or coding sequences from any organism (plants, animals, fungi, bacteria)",
     "DNA - transcripts, mRNA or coding sequences"),
    ("Bacterial, archaeal or viral genomes, contigs or genes", "DNA - bacterial, archaeal or viral ..."),
    ("Already translated proteins", "Protein sequences"),
], columns=["Your sequences", "Choose"])

CONFIDENCE_TEXT = """
Each protein gets a confidence level from the **best ortholog match** found for it:

- **High** - passes all three High rules (E-value, identity, coverage).
- **Medium** - not High, but passes all three Medium rules.
- **Low** - a match passed the search cut-offs but not the Medium rules (for example a partial match).
- **No hit** - no match with the chosen cut-offs.

The rules are set in the sidebar. They are a transparent rule of thumb for how closely your sequence resembles a
sequence of known function. They are **not** a validated probability that the prediction is correct, and the
starting values are common practice rather than a universal standard - state the values you used in your methods.
"""

ENRICHMENT_TEXT = """
Imagine 5% of all genes in a genome are involved in glycolysis. If 12 of the 40 genes you are interested in are,
that is far more than the 2 you would expect by chance, so glycolysis is **enriched** in your group.

The test gives a p-value for every GO term and KEGG pathway. Because hundreds of terms are tested at once, the
p-values are corrected (**FDR**, Benjamini-Hochberg). Your group is compared with **all annotated sequences of the
same run**, so the result is only meaningful when that run is a sensible background (for example a whole
transcriptome) and your group is a subset of it.
"""

LIMITATIONS_TEXT = """
- Predictions are **orthology-based and computational**. They are hypotheses, not experimental evidence.
- **Raw genomic DNA of plants, animals and fungi is not supported directly**: their genes contain introns and need a
  dedicated gene predictor (for example AUGUSTUS) first. Upload predicted proteins or transcripts instead.
- For transcript input only the **single longest open reading frame** of each sequence is analysed.
- Sequences without a close relative of known function get **No hit** or **Low** confidence. This is expected for
  novel, fast-evolving or fragmentary sequences.
- GO terms are transferred from orthologs. With the default setting only non-electronic evidence is used, which gives
  fewer GO terms than KEGG hits.
- eggNOG's database (5.0.2) is older than the current Gene Ontology and KEGG releases, so a few terms may be obsolete
  or renamed.
"""

FAQ = [
    ("Why is the first run slow?",
     "eggNOG's database (about 48 GB) must first be read from disk. This can take about 20 minutes. Later runs are "
     "faster while the database stays in the computer's memory."),
    ("It says another analysis is running.",
     "To protect the server only one analysis runs at a time. Wait a few minutes and press Run again."),
    ("My file is rejected because it has too many sequences.",
     "There is a limit per run (shown in the message). Split the file into parts and run them one after another."),
    ("Why are there fewer GO terms than KEGG hits?",
     "By default only GO terms with non-electronic evidence are transferred. Choose GO evidence 'All' in the sidebar "
     "for more terms (lower reliability)."),
    ("Why do many sequences have No hit or Low confidence?",
     "No close relative of known function was found, or only a partial match. Check the sequence type, try a more "
     "sensitive search mode, and remember that fragments and novel genes often have no good match."),
    ("How do I get my results back later?",
     "Use the Run ID shown after the run (sidebar of Advanced analysis: 'Open a run by its ID'). " + _KEPT +
     " Download the CSV or Excel files for a permanent copy."),
    ("What happens to my sequences?",
     "They are analysed on this server only. The only things fetched from the internet are your NCBI search (if you use "
     "menu option 3, NCBI fetch) and the public KEGG / Gene Ontology name lists. " + _KEPT),
    ("The NCBI search found 0 sequences.",
     "Simplify the search term, for example remove the part after the last AND, and check the organism spelling."),
]

NCBI_EXAMPLES = {
    "Mouse mRNA (animal)": ("Mus musculus[Organism] AND biomol_mrna[PROP] AND srcdb_refseq[PROP] AND 500:3000[Sequence Length]",
                            "nucleotide", "DNA - transcripts, mRNA or coding sequences"),
    "Fruit fly mRNA (animal)": ("Drosophila melanogaster[Organism] AND biomol_mrna[PROP] AND srcdb_refseq[PROP] AND 500:3000[Sequence Length]",
                                "nucleotide", "DNA - transcripts, mRNA or coding sequences"),
    "Arabidopsis mRNA (plant)": ("Arabidopsis thaliana[Organism] AND biomol_mrna[PROP] AND srcdb_refseq[PROP] AND 500:3000[Sequence Length]",
                                 "nucleotide", "DNA - transcripts, mRNA or coding sequences"),
    "E. coli sequences (bacterium)": ("Escherichia coli[Organism] AND 500:3000[Sequence Length]",
                                      "nucleotide", "DNA - bacterial, archaeal or viral genomes / contigs / genes"),
}

CITATIONS = [
    "Cantalapiedra CP, Hernandez-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation, "
    "orthology assignments, and domain prediction at the metagenomic scale. Mol Biol Evol. 2021;38(12):5825-5829.",
    "Huerta-Cepas J, et al. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource "
    "based on 5090 organisms and 2502 viruses. Nucleic Acids Res. 2019;47(D1):D309-D314.",
    "Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. "
    "Nat Methods. 2021;18:366-368.",
    "Hyatt D, et al. Prodigal: prokaryotic gene recognition and translation initiation site identification. "
    "BMC Bioinformatics. 2010;11:119.",
    "Ashburner M, et al. Gene Ontology: tool for the unification of biology. Nat Genet. 2000;25:25-29.",
    "Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 2000;28:27-30.",
    "Klopfenstein DV, et al. GOATOOLS: a Python library for Gene Ontology analyses. Sci Rep. 2018;8:10872.",
    "Cock PJA, et al. Biopython: freely available Python tools for computational molecular biology and "
    "bioinformatics. Bioinformatics. 2009;25:1422-1423.",
    "Benjamini Y, Hochberg Y. Controlling the false discovery rate. J R Stat Soc B. 1995;57:289-300.",
]

COLUMNS = [
    ("sequence_id", "ID of the protein analysed. Bacterial DNA: contig ID + gene number. Transcripts: sequence ID + _ORF1."),
    ("source_sequence", "ID of the sequence you uploaded that this protein came from."),
    ("confidence", "High / Medium / Low / No hit (see the Help tab and the sidebar thresholds)."),
    ("evalue", "Expected number of matches this good by chance. Smaller = more reliable."),
    ("bitscore", "Alignment score that does not depend on database size. Larger = better."),
    ("identity_pct", "Percent of identical residues in the alignment to the matched ortholog."),
    ("coverage_pct", "The smaller of query and subject coverage - the value used for the confidence level."),
    ("query_cov_pct", "Percent of YOUR sequence covered by the alignment (Full view)."),
    ("subject_cov_pct", "Percent of the matched eggNOG sequence covered by the alignment (Full view)."),
    ("gene_name", "Gene name of the matched ortholog group, when known."),
    ("description", "Functional description of the matched ortholog group."),
    ("GO_molecular_function", "GO terms for what the protein does at the molecular level, with names."),
    ("GO_biological_process", "GO terms for the larger process the protein takes part in, with names."),
    ("GO_cellular_component", "GO terms for where in the cell the protein is found, with names."),
    ("KEGG_KO", "KEGG Orthology IDs (with names) of the matched ortholog group."),
    ("KEGG_pathway", "KEGG pathways (with names) linked to the matched ortholog group."),
    ("KEGG_module", "KEGG modules (functional units inside pathways). Full view."),
    ("cog_category / COG_description", "COG functional category letters and their meaning. Full view."),
    ("EC", "Enzyme Commission number(s). Full view."),
    ("PFAM_domains", "Protein domains found in the matched ortholog group. Full view."),
    ("seed_ortholog", "Identifier of the best eggNOG match. Full view."),
    ("taxonomic_level", "Evolutionary level of the ortholog group used for annotation. Full view."),
    ("annotation_status", "Annotated / Ortholog hit, no GO/KEGG / No ortholog hit. Full view."),
    ("gene_call, start, end, strand", "How the protein was found and where it lies on the source sequence. Full view."),
    ("aspect, GO_ID, GO_term, definition", "GO tables: category, GO identifier, term name and definition."),
    ("n_sequences, n_high_confidence", "Dictionary tables: number of sequences (and High-confidence ones) with that term."),
]


def glossary_df() -> pd.DataFrame:
    return pd.DataFrame([(label, default, text) for _k, label, default, text in SETTINGS
                         if _k not in ("source", "upload", "paste", "ncbi_name", "ncbi_max", "enr_file")],
                        columns=["Setting", "Default", "What it does"])


def columns_df() -> pd.DataFrame:
    return pd.DataFrame(COLUMNS, columns=["Column", "Meaning"])