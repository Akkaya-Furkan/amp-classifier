# Antimicrobial peptide classifier

**Question.** Can antimicrobial activity be predicted from sequence features?

## Data
Public AMP databases (e.g. APD3, DBAASP) and a documented negative set.

## Approach
1. Data collection and de-duplication
2. Composition and physicochemical features (charge, hydrophobicity)
3. Sequence-similarity-aware split (CD-HIT clusters)
4. Model training and honest evaluation

## Deliverables
Model, evaluation table, figure, write-up.

## Limitations
- The choice of negatives can inflate performance; how they were chosen is stated.

## Status
Planned — part of a 30-project computational biology and scientific software portfolio.
