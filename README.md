# Sequence Function Predictor

A Streamlit-based web application for ortholog-based functional annotation of biological sequences.

## Features
- Protein FASTA, DNA transcript/CDS and genome DNA workflows.
- Functional annotation using eggNOG-mapper.
- Per-protein evidence: E-value, bit score, identity and coverage.
- Heuristic High, Medium and Low confidence categories.
- Gene Ontology (GO) and KEGG annotations.
- GO and KEGG enrichment analysis.
- CSV and Excel exports.
- Automatic Methods paragraph generation.
- NCBI sequence retrieval.

## Requirements
- Python
- eggNOG-mapper
- DIAMOND
- Prodigal
- eggNOG database

Install Python packages:

    pip install -r requirements.txt

## Run

    streamlit run app.py

## Configuration

Application settings are available in scripts/appconfig.py.
The Gene Ontology ontology file is downloaded automatically when needed.

## Scientific note

Confidence categories are heuristic classifications, not probabilities. Functional annotations are computational predictions and should be validated appropriately.

## License

To be determined after institutional authorization.
