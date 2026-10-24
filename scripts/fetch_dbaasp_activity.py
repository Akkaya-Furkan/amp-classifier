"""Fetch per-peptide antimicrobial activity records from DBAASP for the candidate positives.

Output: data/raw/dbaasp_activity.json  {dbaasp_id: [[measure, group, concentration, unit], ...]}
Resumable: re-running only fetches ids that are missing.
"""
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from amp.data import candidate_dbaasp  # noqa: E402

RAW = ROOT / "data" / "raw"
OUT = RAW / "dbaasp_activity.json"


def fetch(pid):
    for _ in range(4):
        try:
            r = requests.get(f"https://dbaasp.org/peptides/{pid}", headers={"accept": "application/json"}, timeout=90)
            r.raise_for_status()
            d = r.json()
            groups = [g.get("name") for g in d.get("targetGroups") or []]
            acts = [[a.get("activityMeasureValue"), (a.get("targetSpecies") or {}).get("name"),
                     a.get("concentration"), (a.get("unit") or {}).get("name")]
                    for a in d.get("targetActivities") or []]
            return pid, {"groups": groups, "acts": acts}
        except Exception:
            pass
    return pid, None


def main():
    ids = [r["id"] for r in candidate_dbaasp(RAW)]
    done = json.loads(OUT.read_text()) if OUT.exists() else {}
    todo = [i for i in ids if str(i) not in done]
    print("candidates", len(ids), "todo", len(todo), flush=True)
    with ThreadPoolExecutor(24) as ex:
        for n, (pid, res) in enumerate(ex.map(fetch, todo), 1):
            if res is not None:
                done[str(pid)] = res
            if n % 500 == 0:
                OUT.write_text(json.dumps(done)); print(n, flush=True)
    OUT.write_text(json.dumps(done))
    print("fetched", len(done), "of", len(ids))


if __name__ == "__main__":
    main()
