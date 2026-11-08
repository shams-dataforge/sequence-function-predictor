
from Bio import Entrez
import time


def download_sequences_batched(
    search_term: str,
    output_path: str,
    database: str = "nucleotide",
    retmax: int = 300,
    batch_size: int = 50,
    email: str = ""
):
    if not email or "@" not in email:
        raise ValueError("Please provide a valid email address for NCBI.")

    Entrez.email = email.strip()

    search_handle = Entrez.esearch(
        db=database,
        term=search_term,
        retmax=retmax
    )

    results = Entrez.read(search_handle)
    search_handle.close()

    id_list = results["IdList"]

    if not id_list:
        return 0

    all_records = ""

    for i in range(0, len(id_list), batch_size):
        batch_ids = id_list[i:i + batch_size]

        fetch_handle = Entrez.efetch(
            db=database,
            id=batch_ids,
            rettype="fasta",
            retmode="text"
        )

        all_records += fetch_handle.read()
        fetch_handle.close()

        time.sleep(0.4)

    with open(output_path, "w") as f:
        f.write(all_records)

    return len(id_list)
