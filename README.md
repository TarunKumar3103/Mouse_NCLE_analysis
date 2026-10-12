# Mouse NCLE–NED association analysis

Tests whether mouse proteins with a native non-covalent lasso entanglement (NCLE) are more likely to show non-exponential degradation (NED vs ED in McShane et al., *Cell* 2016, NIH-3T3 cells). The repo contains the full pipeline, frozen input snapshots, results, and two independent verifications.

## Results

| Analysis | N | NCLE odds ratio (length-adjusted) | 95% CI | P |
|---|---|---|---|---|
| Primary: NED vs ED | 788 | 1.49 | 0.95–2.32 | 0.082 |
| Sensitivity: NED vs ED + UN | 1,302 | 1.38 | 0.89–2.12 | 0.146 |

The point estimates are positive, but every confidence interval includes 1. No statistically significant association was detected; a moderate positive effect isn't ruled out either. Crude Fisher tests and a High-quality-detection-only sensitivity check give the same picture (see `RESULTS_SUMMARY.txt`).

**Methods in brief.** Model 1 is logistic regression `NED ~ NCLE + z(protein length)`. The crude analysis uses a Fisher exact P with a Woolf CI. Cohort filters:

- exact UniProt match to the native entanglement archive
- binary NCLE call, with any covalent lasso excluded
- not knotted
- no UniProt Membrane (KW-0472) or Secreted (KW-0964) keyword
- whole-protein mean pLDDT ≥ 70, computed from AlphaFold v6 per-residue confidence

Full details are in `README.txt`.

## Reproduce

Python 3.12 and the pinned packages in `requirements.txt` (statsmodels is needed only for the independent check). Everything runs offline from the frozen files in `raw_sources/`.

~~~
python -m pip install -r requirements.txt statsmodels

python scripts/prepare_mouse.py \
  --native raw_sources/Mouse_Native_HQ_GE.zip \
  --mcshane raw_sources/McShane_mouse_mmc1.xlsx \
  --out results/local_preparation
python scripts/prepare_mouse_confidence_inputs.py \
  --prepared results/local_preparation/mcshane_native_exact_merge.csv \
  --native-annotations results/local_preparation/native_protein_annotations.csv \
  --knots raw_sources/knot_all.txt \
  --alphalasso raw_sources/AlphaLasso.zip \
  --keywords raw_sources/uniprot_mouse_keyword_snapshot_2026_10_10.tsv \
  --out results/local_confidence_preparation
python scripts/analyze_mouse.py \
  --prepared results/local_preparation/mcshane_native_exact_merge.csv \
  --metadata raw_sources/alphafold_metadata.csv \
  --confidence raw_sources/residue_confidence_means.csv \
  --knots raw_sources/knot_all.txt \
  --alphalasso raw_sources/AlphaLasso.zip \
  --keywords raw_sources/uniprot_mouse_keyword_snapshot_2026_10_10.tsv \
  --out results/local_replay
python scripts/verify_mouse.py --root . \
  --results-dir results/local_replay --verification-dir verification/local_replay
python scripts/independent_checks.py --mouse .
~~~

Expected output:

- `verify_mouse.py` ends with `"Numerical_failures": 0` (88 fields).
- `independent_checks.py` reports **12/12 independent checks passed**.

`independent_checks.py` rebuilds the cohort from the raw files without importing any package code and refits the models with statsmodels. It also confirms that your rerun matches the shipped tables. Last-digit float rounding differences between platforms are reported but accepted. To rerun, delete the `local_*` output folders first; the scripts refuse to overwrite them.

## Layout

| Folder | Contents |
|---|---|
| `scripts/` | Pipeline, network fetchers (not needed to reproduce), verifiers |
| `raw_sources/` | Frozen inputs: McShane supplement, native NCLE archive, knots, AlphaLasso, UniProt keyword snapshot, AlphaFold v6 metadata and confidence JSON |
| `preparation/` | Intermediate tables |
| `results/` | Shipped result tables |
| `verification/` | Shipped verification outputs |
| `Mouse_NCLE_NED_association_results.xlsx` | Formatted workbook |

`SOURCE_MANIFEST.json` lists SHA-256 hashes for every file.

## Limitations

- **Isoform IDs:** exact accession matching drops all 373 McShane isoform IDs; no remapping is done.
- **Keyword snapshot:** the UniProt keywords are a 2026_03 snapshot, not the annotations McShane used.
- **Covalent lasso calls:** the native entanglement detector itself was not rerun; its supplied annotations were audited.
- **Interpretation:** NCLE status doesn't directly measure misfolding, and the ORs are associations, not causal effects.

## Sources

McShane et al., *Cell* 167, 803–815 (2016), https://doi.org/10.1016/j.cell.2016.09.015 · AlphaFold DB, https://alphafold.ebi.ac.uk · UniProt, https://www.uniprot.org
