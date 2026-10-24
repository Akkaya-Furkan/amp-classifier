"""Download public source data into data/raw (APD3, DBAASP, UniProt Swiss-Prot)."""
import json
import sys
import time
from pathlib import Path

import requests
from concurrent.futures import ThreadPoolExecutor

RAW = Path(__file__).resolve().parents[1] / "data" / "raw"
APD_URL = "https://aps.unmc.edu/assets/sequences/APD_sequence_release_09142020.fasta"
DBAASP_URL = "https://dbaasp.org/peptides"
PAGE = 200
UNIPROT_STREAM = "https://rest.uniprot.org/uniprotkb/stream"

# Exclusion terms keep annotated antimicrobial/toxin-like proteins out of the negatives.
EXCLUDE = (
    "NOT (keyword:KW-0929) NOT (keyword:KW-0044) NOT (keyword:KW-0211) "
    "NOT (keyword:KW-0078) NOT (keyword:KW-0295) NOT (keyword:KW-0800) "
    "NOT (keyword:KW-0081) NOT (keyword:KW-0930) "
    "NOT antimicrobial NOT antibacterial NOT antifungal NOT bacteriocin "
    "NOT defensin NOT antibiotic NOT \"host defense\" NOT toxin NOT hemolysin "
    "NOT cathelicidin NOT lectin NOT microcin NOT lantibiotic"
)
UNIPROT_QUERIES = {
    # secreted short Swiss-Prot entries: hardest realistic negatives
    "uniprot_secreted.fasta": f"(reviewed:true) AND (length:[10 TO 100]) AND (keyword:KW-0964) {EXCLUDE}",
    # secreted Swiss-Prot proteins: random windows of these are the length-unconstrained negatives
    "uniprot_secreted_proteins.fasta": f"(reviewed:true) AND (length:[100 TO 800]) AND (keyword:KW-0964) {EXCLUDE}",
}


def get(url, **kw):
    for attempt in range(5):
        try:
            r = requests.get(url, timeout=120, **kw)
            r.raise_for_status()
            return r
        except requests.RequestException as e:
            print("retry", attempt, e, file=sys.stderr)
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"failed: {url}")


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    (RAW / "apd3.fasta").write_bytes(get(APD_URL).content)

    # DBAASP paginates with offset/limit (a `page` parameter is silently ignored).
    hdr = {"accept": "application/json"}
    first = get(DBAASP_URL, params={"offset": 0, "limit": PAGE}, headers=hdr).json()
    total = first["totalCount"]
    n_pages = -(-total // PAGE)

    def page(i):
        return get(DBAASP_URL, params={"offset": i * PAGE, "limit": PAGE}, headers=hdr).json()["data"]

    with ThreadPoolExecutor(8) as ex:
        rows = [r for chunk in ex.map(page, range(n_pages)) for r in chunk]
    ids = [r["id"] for r in rows]
    assert len(set(ids)) == len(ids), "duplicate DBAASP ids: pagination broken"
    d = {"totalCount": total}
    (RAW / "dbaasp.json").write_text(json.dumps(rows))
    print("dbaasp records:", len(rows), "of", d["totalCount"])

    for name, q in UNIPROT_QUERIES.items():
        r = get(UNIPROT_STREAM, params={"query": q, "format": "fasta"})
        (RAW / name).write_bytes(r.content)
        print(name, r.text.count(">"))


if __name__ == "__main__":
    main()
