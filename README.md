# Antimicrobial peptide classifier

**Question.** Can antimicrobial activity be predicted from sequence features alone, and how much of the
apparent accuracy survives once train and test peptides are properly separated by sequence similarity?

## Data (all public)
| Role | Source | Selection |
|---|---|---|
| Positives | [APD3](https://aps.unmc.edu/) release 09-14-2020 + [DBAASP](https://dbaasp.org/) monomers | 20 standard amino acids, length 10–100, no N-terminal modification, exact-duplicate removal. DBAASP records kept only if DBAASP lists a MIC ≤ 32 µg/ml (µM converted at 110 Da/residue); APD3 entries are taken as curated AMPs. **10,131** unique positives. |
| Negatives A: *secreted peptides* | UniProt Swiss-Prot, reviewed, secreted (KW-0964), length 10–100 | Entries annotated antimicrobial/antibiotic/defensin/bacteriocin/fungicide/toxin etc. excluded (query in `scripts/download_data.py`); entries starting with Met dropped (mostly unprocessed precursors); sampled so the length histogram equals the positives' **at every exact length** (1,279 vs 1,279). |
| Negatives B: *secreted protein fragments* | Windows cut from the same kind of Swiss-Prot secreted proteins (100–800 aa) | One window per positive, exact same length, first 30 residues skipped (signal peptides), each from a different parent protein (10,130 windows). |

Negatives are *presumed* inactive, not experimentally tested — see Limitations.

## Approach
1. `scripts/download_data.py`, `scripts/fetch_dbaasp_activity.py` — download raw data (`data/raw/`, not committed).
2. Features (`amp/features.py`): amino-acid composition, length, net charge at pH 7.4, pI, mean hydropathy,
   hydrophobic moment, class fractions; dipeptide composition for the strongest model.
3. Models (`amp/models.py`): logistic regression, random forest, histogram gradient boosting, plus **length-only controls** (linear and non-linear).
4. **Similarity-aware evaluation** (`amp/split.py`): all sequences are clustered with CD-HIT (40 % identity, word size 2); 5-fold stratified
   cross-validation holds out whole clusters. CD-HIT alone proved leaky (below), so every fold's test rows that have a
   strong local-alignment relative in the training part (shared 5-mer and BLOSUM62 alignment identity ≥ 0.8 of the shorter or ≥ 0.5 of the
   longer sequence) are **purged**. A random 5-fold split is reported alongside for contrast.
5. Metrics on pooled out-of-fold predictions: ROC-AUC and PR-AUC (average precision), 95 % CIs from a bootstrap over CD-HIT clusters.

## Results
Full table: [`results/evaluation_table.md`](results/evaluation_table.md) (CSV alongside); figure: [`results/evaluation.png`](results/evaluation.png).

Best model (dipeptide+AAC+physchem, gradient boosting), ROC-AUC / PR-AUC with 95 % CI:

| Split | Negatives | n test | Positive fraction | ROC-AUC | PR-AUC |
|---|---|---|---|---|---|
| Random 5-fold | secreted peptides | 2,558 | 0.50 | 0.971 [0.964–0.977] | 0.975 [0.968–0.981] |
| **CD-HIT + purge** | secreted peptides | 1,035 | 0.38 | **0.906 [0.881–0.926]** | **0.866 [0.817–0.906]** |
| Random 5-fold | protein fragments | 20,261 | 0.50 | 0.986 [0.983–0.988] | 0.988 [0.985–0.991] |
| **CD-HIT + purge** | protein fragments | 10,357 | 0.21 | **0.947 [0.939–0.954]** | **0.884 [0.862–0.905]** |
| CD-HIT + purge, train fragments → test peptides | | 1,035 | 0.38 | 0.900 [0.876–0.923] | 0.862 [0.821–0.900] |
| CD-HIT + purge, train peptides → test fragments | | 10,357 | 0.21 | 0.921 [0.911–0.930] | 0.809 [0.779–0.836] |

- Random splitting overstates ROC-AUC by about 0.065 (secreted peptides) and 0.04 (fragments) compared to the purged split.
- Simpler models are close: with secreted-peptide negatives, physicochemical features with logistic regression reach 0.884 ROC-AUC; adding composition 0.893; random forest 0.902.
- Length-only controls are at or below chance (ROC-AUC 0.42 linear, 0.27 gradient boosting on the purged peptide set), so length matching removes the
  trivial shortcut. The below-0.5 values come from distribution shift after purging, not from signal.
- PR-AUC depends on the class balance of the surviving test rows (purging removes more positives, which come in analogue families, than negatives);
  compare PR-AUC only between rows with the same positive fraction.

### How well is the split separated?
`results/split_audit.json` (`scripts/audit_split.py`): nearest-neighbour identity of 100 sampled test sequences per fold to the training
part, mean over 5 folds (identity / longer sequence, median): random split 0.56, CD-HIT only 0.49, **CD-HIT + purge 0.34**, shuffled-sequence null 0.31–0.34.
Share of test sequences with identity ≥ 0.6: 0.48 / 0.41 / 0.03, null 0.00–0.06. CD-HIT at 40 % therefore removes little on its own for short peptides
(it only compares sequences with cluster representatives); after purging, test sequences look as unrelated to the training part as shuffled sequences do.
Single-linkage grouping was tried first and collapsed 90 % of the data into one component, hence purging.

## Limitations
- **Negative set choice inflates performance and cannot be removed.** Both negative sets are *assumed* inactive (non-antimicrobial annotation in Swiss-Prot),
  not assayed. They differ from AMPs in origin and processing, so a classifier may partly learn "database/annotation type". Fragments are not real peptides;
  secreted peptides are few (1,279 after matching) and hard to match on every property. Performance against experimentally confirmed inactive peptides is **not measured**
  (not done).
- Peptide and fragment settings use different positive subsets (1,279 vs 10,131 before purging), so their scores differ for more than the negatives.
- Positives: APD3 is not re-validated; DBAASP activity uses a single MIC ≤ 32 µg/ml rule over mixed targets and conditions. N-terminal-modified records are dropped; C-terminal
  modifications (e.g. amidation) are ignored in features and in duplicate removal.
- Purging removes about half of each fold's test rows (≈2,100 of ≈4,300), preferentially analogue-rich positives, so the surviving test sets are biased toward singleton sequences.
  Candidate pairs require a shared exact 5-mer; heavily diverged near-duplicates sharing none are missed.
- Bootstrap CIs resample CD-HIT clusters on fixed out-of-fold predictions; they ignore training and split variability (one split seed).
- Fragment negatives skip the first 30 residues and come from proteins ≥ window+40 aa, a possible positional/parent-length shortcut. Not checked.
- Metrics use roughly balanced classes; precision at real screening prevalence (AMPs ≪ 1 %) would be far lower. No wet-lab validation, no hemolysis/toxicity prediction.

## Reproduce
```bash
conda env create -f environment.yml && conda activate amp-classifier
python scripts/download_data.py            # ~5 min
python scripts/fetch_dbaasp_activity.py    # ~35 min, resumable
python scripts/run_experiment.py           # builds data/processed/, writes results/
python scripts/audit_split.py
pytest
```
`data/processed/dataset.csv` (sequences, labels, source, CD-HIT cluster) and `similarity_edges.npy` are committed; raw downloads are not.
Source databases have their own terms of use; check them before redistributing the raw data.

## Layout
`amp/` features, data assembly, splits, models, audit, evaluation · `scripts/` pipeline · `tests/` unit tests · `results/` tables, figure, split audit.
