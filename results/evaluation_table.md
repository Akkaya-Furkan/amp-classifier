| split                        | negatives                        | model                        |   n_test |   prevalence | ROC-AUC [95% CI]    | PR-AUC (avg. precision) [95% CI]   |
|:-----------------------------|:---------------------------------|:-----------------------------|---------:|-------------:|:--------------------|:-----------------------------------|
| cluster (CD-HIT 40% + purge) | secreted peptides                | length_only (LR)             |     1035 |     0.381643 | 0.421 [0.372-0.472] | 0.350 [0.291-0.414]                |
| cluster (CD-HIT 40% + purge) | secreted peptides                | length_only (HGB)            |     1035 |     0.381643 | 0.266 [0.230-0.306] | 0.272 [0.236-0.323]                |
| cluster (CD-HIT 40% + purge) | secreted peptides                | physchem (LR)                |     1035 |     0.381643 | 0.884 [0.857-0.908] | 0.851 [0.807-0.886]                |
| cluster (CD-HIT 40% + purge) | secreted peptides                | AAC+physchem (LR)            |     1035 |     0.381643 | 0.893 [0.869-0.914] | 0.854 [0.815-0.889]                |
| cluster (CD-HIT 40% + purge) | secreted peptides                | AAC+physchem (RF)            |     1035 |     0.381643 | 0.902 [0.876-0.923] | 0.864 [0.817-0.904]                |
| cluster (CD-HIT 40% + purge) | secreted peptides                | dipeptide+AAC+physchem (HGB) |     1035 |     0.381643 | 0.906 [0.881-0.926] | 0.866 [0.817-0.906]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | length_only (LR)             |    10357 |     0.205175 | 0.503 [0.468-0.542] | 0.236 [0.209-0.274]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | length_only (HGB)            |    10357 |     0.205175 | 0.404 [0.369-0.432] | 0.181 [0.165-0.194]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | physchem (LR)                |    10357 |     0.205175 | 0.918 [0.906-0.927] | 0.811 [0.784-0.835]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | AAC+physchem (LR)            |    10357 |     0.205175 | 0.926 [0.917-0.935] | 0.822 [0.800-0.842]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | AAC+physchem (RF)            |    10357 |     0.205175 | 0.943 [0.935-0.951] | 0.866 [0.842-0.890]                |
| cluster (CD-HIT 40% + purge) | secreted protein fragments       | dipeptide+AAC+physchem (HGB) |    10357 |     0.205175 | 0.947 [0.939-0.954] | 0.884 [0.862-0.905]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | length_only (LR)             |     1035 |     0.381643 | 0.527 [0.475-0.583] | 0.443 [0.383-0.521]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | length_only (HGB)            |     1035 |     0.381643 | 0.418 [0.372-0.462] | 0.367 [0.329-0.419]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | physchem (LR)                |     1035 |     0.381643 | 0.885 [0.859-0.905] | 0.858 [0.818-0.889]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | AAC+physchem (LR)            |     1035 |     0.381643 | 0.894 [0.870-0.912] | 0.856 [0.821-0.884]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | AAC+physchem (RF)            |     1035 |     0.381643 | 0.907 [0.884-0.926] | 0.878 [0.843-0.905]                |
| cluster (CD-HIT 40% + purge) | train fragments -> test peptides | dipeptide+AAC+physchem (HGB) |     1035 |     0.381643 | 0.900 [0.876-0.923] | 0.862 [0.821-0.900]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | length_only (LR)             |    10357 |     0.205175 | 0.475 [0.442-0.508] | 0.221 [0.193-0.253]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | length_only (HGB)            |    10357 |     0.205175 | 0.443 [0.420-0.470] | 0.200 [0.178-0.228]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | physchem (LR)                |    10357 |     0.205175 | 0.910 [0.897-0.920] | 0.787 [0.758-0.813]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | AAC+physchem (LR)            |    10357 |     0.205175 | 0.909 [0.897-0.920] | 0.773 [0.747-0.798]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | AAC+physchem (RF)            |    10357 |     0.205175 | 0.918 [0.907-0.928] | 0.801 [0.772-0.828]                |
| cluster (CD-HIT 40% + purge) | train peptides -> test fragments | dipeptide+AAC+physchem (HGB) |    10357 |     0.205175 | 0.921 [0.911-0.930] | 0.809 [0.779-0.836]                |
| random                       | secreted peptides                | length_only (LR)             |     2558 |     0.5      | 0.477 [0.452-0.497] | 0.484 [0.441-0.528]                |
| random                       | secreted peptides                | length_only (HGB)            |     2558 |     0.5      | 0.325 [0.305-0.343] | 0.389 [0.353-0.425]                |
| random                       | secreted peptides                | physchem (LR)                |     2558 |     0.5      | 0.935 [0.921-0.947] | 0.945 [0.930-0.957]                |
| random                       | secreted peptides                | AAC+physchem (LR)            |     2558 |     0.5      | 0.947 [0.938-0.957] | 0.954 [0.944-0.963]                |
| random                       | secreted peptides                | AAC+physchem (RF)            |     2558 |     0.5      | 0.965 [0.957-0.971] | 0.969 [0.960-0.975]                |
| random                       | secreted peptides                | dipeptide+AAC+physchem (HGB) |     2558 |     0.5      | 0.971 [0.964-0.977] | 0.975 [0.968-0.981]                |
| random                       | secreted protein fragments       | length_only (LR)             |    20261 |     0.500025 | 0.494 [0.480-0.506] | 0.495 [0.464-0.524]                |
| random                       | secreted protein fragments       | length_only (HGB)            |    20261 |     0.500025 | 0.459 [0.451-0.466] | 0.466 [0.438-0.497]                |
| random                       | secreted protein fragments       | physchem (LR)                |    20261 |     0.500025 | 0.959 [0.954-0.963] | 0.963 [0.956-0.969]                |
| random                       | secreted protein fragments       | AAC+physchem (LR)            |    20261 |     0.500025 | 0.968 [0.965-0.971] | 0.971 [0.966-0.976]                |
| random                       | secreted protein fragments       | AAC+physchem (RF)            |    20261 |     0.500025 | 0.983 [0.981-0.986] | 0.986 [0.983-0.989]                |
| random                       | secreted protein fragments       | dipeptide+AAC+physchem (HGB) |    20261 |     0.500025 | 0.986 [0.983-0.988] | 0.988 [0.985-0.991]                |