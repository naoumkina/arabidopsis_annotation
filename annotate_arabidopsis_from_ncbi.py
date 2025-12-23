import time
import re
import requests
import pandas as pd
from Bio import Entrez
from lxml import etree

# ==========================
# CONFIGURATION
# ==========================

# NCBI requires an email, and recommends an API key if you have one.
Entrez.email = "marina.naoumkina@usda.gov"  # CHANGE THIS
# Entrez.api_key = "YOUR_NCBI_API_KEY"   # Optional but strongly recommended

# Rate limiting: be polite to NCBI
NCBI_DELAY = 0.4  # seconds between requests; adjust if you hit limits

# Stress-related keyword list (case-insensitive)
STRESS_KEYWORDS = [
    "stress", "drought", "salt", "salinity", "heat", "cold", "chilling",
    "freezing", "osmotic", "oxidative", "ros", "abscisic", "aba",
    "desiccation", "flood", "hypoxia", "wounding", "pathogen",
    "biotic", "abiotic"
]

# Maximum number of PubMed articles to pull titles for (per gene)
MAX_PUBMED_TITLES = 20


# ==========================
# HELPER FUNCTIONS
# ==========================

def search_ncbi_gene_id(agi):
    """
    Given an AGI (e.g., AT5G54510), use Entrez.esearch to find the corresponding
    NCBI Gene ID in Arabidopsis thaliana.
    """
    query = f"{agi}[Gene Name] AND Arabidopsis thaliana[Organism]"
    handle = Entrez.esearch(db="gene", term=query)
    record = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_DELAY)

    id_list = record.get("IdList", [])
    if not id_list:
        return None
    # Usually the first hit is the one you want for AGIs
    return id_list[0]


def get_gene_summary(ncbi_gene_id):
    handle = Entrez.esummary(db="gene", id=ncbi_gene_id)
    record = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_DELAY)

    # NCBI now returns a dict with nested structure
    docs = []

    if isinstance(record, dict):
        docs = record.get("DocumentSummarySet", {}).get("DocumentSummary", [])
    elif isinstance(record, list):
        docs = record
    else:
        docs = []

    # If no documents found, return empty fields
    if not docs:
        return {"symbol": None, "description": None, "summary": None}

    doc = docs[0]

    symbol = doc.get("NomenclatureSymbol") or doc.get("Name")
    description = doc.get("Description")
    summary = doc.get("Summary")

    return {
        "symbol": symbol,
        "description": description,
        "summary": summary
    }


def get_pubmed_links(ncbi_gene_id, retries=3):
    """
    Get PubMed IDs linked to this gene using Entrez.elink (gene -> pubmed).
    Includes retry logic to handle NCBI connection drops.
    """
    for attempt in range(retries):
        try:
            handle = Entrez.elink(dbfrom="gene", db="pubmed", id=ncbi_gene_id)
            record = Entrez.read(handle)
            handle.close()
            time.sleep(NCBI_DELAY)
            break  # success
        except Exception as e:
            print(f"    PubMed link fetch failed (attempt {attempt+1}/{retries}): {e}")
            time.sleep(1.5)
            record = None

    # If still no record, return empty list
    if not record:
        return []

    pubmed_ids = []
    try:
        linkset = record[0]
        for link in linkset.get("LinkSetDb", []):
            if link.get("DbTo") == "pubmed":
                pubmed_ids.extend([l["Id"] for l in link.get("Link", [])])
    except Exception:
        return []

    # Deduplicate while preserving order
    seen = set()
    unique_ids = []
    for pid in pubmed_ids:
        if pid not in seen:
            seen.add(pid)
            unique_ids.append(pid)

    return unique_ids


def fetch_pubmed_titles(pubmed_ids, max_titles=20):
    """
    Given a list of PubMed IDs, fetch the article titles using Entrez.efetch.
    Returns a dict {pmid: title}.
    """
    if not pubmed_ids:
        return {}

    ids_to_fetch = pubmed_ids[:max_titles]
    id_str = ",".join(ids_to_fetch)
    handle = Entrez.efetch(db="pubmed", id=id_str, rettype="xml", retmode="xml")
    records = Entrez.read(handle)
    handle.close()
    time.sleep(NCBI_DELAY)

    pmid_to_title = {}
    for article in records["PubmedArticle"]:
        try:
            pmid = article["MedlineCitation"]["PMID"]
            title = article["MedlineCitation"]["Article"]["ArticleTitle"]
            pmid_to_title[str(pmid)] = str(title)
        except KeyError:
            continue

    return pmid_to_title


def get_go_terms_from_gene(ncbi_gene_id):
    """
    Retrieve GO terms from NCBI Gene record via Entrez.efetch (XML) and parse.
    Returns a list of dicts: [{"go_id": "...", "category": "...", "term": "..."}, ...]
    """
    handle = Entrez.efetch(db="gene", id=ncbi_gene_id, rettype="xml", retmode="xml")
    xml_data = handle.read()
    handle.close()
    time.sleep(NCBI_DELAY)

    # Parse the XML using lxml for easier XPath-like operations
    root = etree.fromstring(xml_data.encode() if isinstance(xml_data, str) else xml_data)

    ns = {
        "e": "http://www.ncbi.nlm.nih.gov"
    }

    go_terms = []

    # NCBI Gene XML structure is a bit nested; GO annotations are in
    # Gene-commentary elements with db="GO".
    for gc in root.xpath("//e:Entrezgene//e:Gene-commentary", namespaces=ns):
        db = gc.find("e:Gene-commentary_source/e:Other-source/e:Other-source_src/e:Dbtag/e:Dbtag_db", namespaces=ns)
        if db is not None and db.text == "GO":
            # GO ID
            go_id_tag = gc.find("e:Gene-commentary_source/e:Other-source/e:Other-source_src/e:Dbtag/e:Dbtag_tag/e:Object-id/e:Object-id_id", namespaces=ns)
            go_id = f"GO:{go_id_tag.text.zfill(7)}" if go_id_tag is not None else None

            # Term name
            go_term_tag = gc.find("e:Gene-commentary_text", namespaces=ns)
            go_term = go_term_tag.text if go_term_tag is not None else None

            # Category (BP/CC/MF) is stored in comment type or label
            category = None
            heading = gc.find("e:Gene-commentary_heading", namespaces=ns)
            if heading is not None and heading.text:
                category = heading.text  # e.g. "Function", "Process", "Component"

            if go_id and go_term:
                go_terms.append({
                    "go_id": go_id,
                    "category": category,
                    "term": go_term
                })

    return go_terms


def flag_stress_related(summary, go_terms, pubmed_titles):
    """
    Heuristic: mark a gene as 'stress-related' if any stress keyword appears
    in gene summary, GO term names, or PubMed titles (case-insensitive).
    Returns (flag_bool, evidence_string).
    """
    text_blobs = []

    if summary:
        text_blobs.append(summary.lower())
    if go_terms:
        text_blobs.extend([gt["term"].lower() for gt in go_terms if gt.get("term")])
    if pubmed_titles:
        text_blobs.extend([t.lower() for t in pubmed_titles.values()])

    combined_text = " ".join(text_blobs)

    hits = []
    for kw in STRESS_KEYWORDS:
        if re.search(rf"\b{kw}\b", combined_text):
            hits.append(kw)

    if hits:
        evidence = "Keywords: " + ", ".join(sorted(set(hits)))
        return True, evidence
    else:
        return False, ""


# ==========================
# MAIN PIPELINE
# ==========================

def annotate_genes(agi_file, output_csv):
    """
    Main driver:
    - Reads AGIs from agi_file
    - For each AGI:
        * find NCBI Gene ID
        * get summary
        * get GO terms
        * get PubMed IDs + titles
        * flag stress-related
    - Writes a clean CSV to output_csv
    """
    with open(agi_file) as f:
        agi_list = [line.strip() for line in f if line.strip()]

    records = []

    for i, agi in enumerate(agi_list, start=1):
        print(f"[{i}/{len(agi_list)}] Processing {agi}...")

        # 1) Find NCBI Gene ID
        ncbi_gene_id = search_ncbi_gene_id(agi)
        if not ncbi_gene_id:
            print(f"  No NCBI Gene ID found for {agi}")
            records.append({
                "AGI": agi,
                "NCBI_GeneID": None,
                "Symbol": None,
                "Description": None,
                "Summary": None,
                "GO_Terms": None,
                "PubMed_IDs": None,
                "PubMed_Titles": None,
                "Stress_Related": False,
                "Stress_Evidence": ""
            })
            continue

        # 2) Gene summary / basic info
        gene_info = get_gene_summary(ncbi_gene_id)
        symbol = gene_info["symbol"]
        description = gene_info["description"]
        summary = gene_info["summary"]

        # 3) GO terms
        go_terms = get_go_terms_from_gene(ncbi_gene_id)
        go_str = "; ".join(
            [f"{g['go_id']} ({g.get('category')}): {g['term']}" for g in go_terms]
        ) if go_terms else None

        # 4) PubMed links + titles
        pubmed_ids = get_pubmed_links(ncbi_gene_id)
        pubmed_id_str = ",".join(pubmed_ids) if pubmed_ids else None

        pmid_to_title = fetch_pubmed_titles(pubmed_ids, max_titles=MAX_PUBMED_TITLES) if pubmed_ids else {}
        pubmed_title_str = "; ".join(
            [f"{pid}: {title}" for pid, title in pmid_to_title.items()]
        ) if pmid_to_title else None

        # 5) Stress flag
        stress_flag, stress_evidence = flag_stress_related(summary, go_terms, pmid_to_title)

        records.append({
            "AGI": agi,
            "NCBI_GeneID": ncbi_gene_id,
            "Symbol": symbol,
            "Description": description,
            "Summary": summary,
            "GO_Terms": go_str,
            "PubMed_IDs": pubmed_id_str,
            "PubMed_Titles": pubmed_title_str,
            "Stress_Related": stress_flag,
            "Stress_Evidence": stress_evidence
        })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"\nAnnotation complete. Wrote {len(df)} records to {output_csv}")


if __name__ == "__main__":
    # Adjust filenames as needed
    annotate_genes("agi_list.txt", "ncbi_annotation_for_arabidopsis_genes.csv")
