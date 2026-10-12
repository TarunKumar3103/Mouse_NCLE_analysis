VERIFIED RESULTS
Primary N=788; sensitivity N=1,302. Main length-adjusted NCLE OR: 1.485521
(95% CI 0.951341-2.319644; P=0.081757). See RESULTS_SUMMARY.txt.

MOUSE NCLE-NED ASSOCIATION: SOURCES, METHODS AND REPRODUCTION

Purpose
Test the association between detected native non-covalent lasso entanglement
(NCLE) and McShane's NED degradation category in mouse proteins. Run the crude
association and Model 1 only: NED ~ NCLE + z(protein length).

Analyses
1. Primary: NED versus ED; UN is excluded.
2. Sensitivity: NED versus ED+UN; UN is not relabeled as exponential decay.
3. Conservative local-quality check: repeat both comparisons requiring a
   non-covalent detection with Quality=High for a positive NCLE call. Proteins
   with only Low detections are excluded from this check rather than assigned
   to the negative group. This is an additional source-quality sensitivity
   analysis, not a change to the primary human-consistent coding rule.

Native coding and source checks
- There are 21,452 unique F1/version-6 model files in the supplied archive.
- Native CSVs are pipe-delimited. All detected rows have ENT=True.
- Any CCbond=True native row excludes that protein from binary NCLE analyses.
- Otherwise any CCbond=False row gives the main NCLE=1 call.
- A header-only file is interpreted as no detected NCLE (NCLE=0). Missing
  files remain unmatched. Detector completion logs were not supplied.
- The HQ folder contains both High and Low annotation rows. The filename
  alone does not establish a whole-protein mean pLDDT of at least 70.
- The main rule matches the earlier human CCbond aggregation. The conservative
  check makes the local-confidence ambiguity visible without treating uncertain
  Low-only detections as confirmed negatives.

Protein quality, length and exclusions
- Whole-protein mean pLDDT is calculated as the arithmetic mean of the exact
  version-6 confidenceScore array, one score per residue. Positions must be
  consecutive 1 through N. The threshold is >=70, including exactly 70.
- The residue-array length is checked against the same model's full sequence
  in the official AlphaFold API. The model ID, version 6, F1 fragment and mouse
  taxonomy 10090 must match. Length is the model's full sequence length.
- The API summary pLDDT is retained for comparison but is not the cutoff field:
  a test found that it can differ from the arithmetic residue mean.
- Exclude supplied knot IDs, and covalent-lasso IDs from any native covalent
  record or AlphaLasso table (Organism=Mouse).
- Exclude UniProt Membrane (KW-0472) and Secreted (KW-0964) keywords. The mouse
  snapshot was retrieved on 2026-10-10 from UniProt release 2026_03.
- Retain meaningful isoform suffixes. Only the terminal AlphaFold -F<number>
  suffix in AlphaLasso records is removed. No canonical/isoform remapping is used.
- Neither MobiDB disorder nor complex membership is used as a filter or covariate.
- Topology and keyword exclusions precede residue-score retrieval to reduce
  downloads. These intersections give the same final cohort as filtering on
  pLDDT first. Source losses and their ordering are recorded in cohort_flow.csv.

Statistical definitions
- Outcome: NED=1; ED=0 in primary. ED or UN=0 in sensitivity.
- Crude OR: odds of NED among NCLE=1 divided by odds among NCLE=0.
- Crude P-value: two-sided Fisher exact test. CI: Woolf log-OR Wald interval.
- Model 1: logistic maximum likelihood, intercept + NCLE + standardized length.
- z(length) uses the mean and sample standard deviation (ddof=1) of each exact
  fitted cohort. NCLE OR is adjusted only for protein length. Length OR is per
  one sample SD. Regression CIs and P-values use Wald normal approximations.
- ORs are associations and do not establish misfolding or causal degradation.
- A non-significant result is not proof of no effect. Evaluate the CI as well.
- The quality checks and UN sensitivity are reported without multiplicity
  adjustment, alongside the designated main analysis.

Where to run the code
Extract this ZIP. Open its Mouse_NCLE_analysis folder in VS Code.
Use Terminal > New Terminal. All commands below are terminal commands, not
Python expressions for pasting directly into zsh or notebook cells.
Run them while your terminal is in the Mouse_NCLE_analysis folder.

Your existing human_ncle_repro environment uses the same analysis dependencies:

conda activate human_ncle_repro
python --version
python -m pip install -r requirements.txt

Alternatively create a separate environment:

conda create -n mouse_ncle_repro python=3.12 -y
conda activate mouse_ncle_repro
python -m pip install -r requirements.txt

Offline reproduction: no online queries are needed
Use the frozen source files already included in raw_sources.

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

python scripts/verify_mouse.py \
  --root . \
  --results-dir results/local_replay \
  --verification-dir verification/local_replay

The independent check reconstructs cohort IDs and source fields from the
original archive, McShane workbook, exclusion sources and raw public responses.
It solves logistic-score equations independently of the analysis fitter and
enumerates the fixed-margin hypergeometric distribution for crude P-values.
It reports the number of comparisons and any numerical failures.

Use a new result directory if rerunning analyze_mouse.py; it does not overwrite
an existing result directory. Open the resulting CSVs in Excel to see the same
scientific tables. The supplied XLSX is the formatted presentation of them.
The included build_workbook.mjs requires the ChatGPT spreadsheet runtime; it
is not required to reproduce any numerical result in a local Python setup.

Code-to-table map
prepare_mouse.py:
  native_protein_annotations.csv: one row per supplied native model.
  mcshane_native_exact_merge.csv: original experimental rows plus exact matches.
  archive_profile_overlap.csv: degradation labels versus native-file availability.
  unmatched_degradation_records.csv: exact accession losses; includes isoforms.

prepare_mouse_confidence_inputs.py:
  confidence_candidate_ids.csv: all 1,653 eligible candidates before the mean
    pLDDT filter, ordered by UniProt_ID; used for residue-confidence retrieval.
  alphalasso_mouse_exact_evidence.csv: 1,128 Mouse-labelled AlphaLasso records
    whose accessions match the full supplied native archive (888 unique IDs).
    These are not limited to the experimental overlap. analyze_mouse.py reads
    the original AlphaLasso archive independently for its final exclusions.
  confidence_preparation_summary.json: candidate rules, counts and source hashes.
  This step is offline. It is not necessary to download AlphaFold again.

fetch_mouse_metadata.py:
  alphafold_metadata.csv and raw_prediction_json/: exact-model identity,
  full sequence length, summary confidence, source body hashes and sequences.

fetch_mouse_confidence.py:
  residue_confidence_means.csv and raw_confidence_json/: exact-v6 residue arrays,
  whole-protein arithmetic means, coordinate and length information.

fetch_mouse_keywords.py:
  uniprot_mouse_keyword_snapshot_2026_10_10.tsv and query metadata: taxonomy,
  current UniProt length, Membrane and Secreted keyword coverage.

analyze_mouse.py + stats_model1.py:
  association_results.csv: crude, NCLE, length and intercept regression estimates.
  contingency_counts.csv: counts, totals and NED fractions for every fitted cohort.
  cohort_flow.csv: sequential exclusions for primary and UN sensitivity.
  primary_any_ncle_cohort.csv / sensitivity_any_ncle_cohort.csv: main cohorts.
  primary_high_only_cohort.csv / sensitivity_high_only_cohort.csv: conservative
    local-quality cohorts, with Low-only detections excluded.
  all_experimental_record_audit.csv: all 3,292 source records and exclusion flags.
  model_diagnostics.csv: convergence, standardization and information checks.
  covalent_lasso_evidence.csv: matching source rows and table names.
  local_quality_definition_differences.csv: Low-only detections for review.
  sequence_length_version_differences.csv: current UniProt versus AF sequence length.
  pLDDT_summary_comparison.csv: API summary versus actual residue means.
  analysis_manifest.json: methods, source hashes, versions and limitations.

verify_mouse.py:
  independent_cohort_checks.csv: exact membership and field comparisons.
  independent_numerical_checks.csv: independent versus reported numerical fields.
  verification_summary.json: verification outcome and number of comparisons.

build_workbook.mjs:
  Results, Cohort flow and Protein records sheets. Workbook count and crude-OR
  formulas are checked against Python. Regression estimates and P-values are
  saved Python fits; rerun the Python scripts after any data edit.

Sources
McShane et al., Cell 2016: https://doi.org/10.1016/j.cell.2016.09.015
Official AlphaFold DB: https://alphafold.ebi.ac.uk/faq
Official AlphaFold API: https://alphafold.ebi.ac.uk/api-docs
UniProt API: https://rest.uniprot.org
Native NCLE, AlphaLasso and knot sources: the supplied original attachments.

PACKAGING REVISION — 11 OCTOBER 2026
The previously undocumented confidence-candidate and native-overlap AlphaLasso
preparation tables now have an executable generator and commands above.
The candidate file is now in ascending UniProt_ID order. All 1,653 ID/model
records match the historical file; only its row order changed. The 1,128-row
AlphaLasso evidence file reproduces byte-for-byte. Scientific analysis scripts,
raw scientific inputs, result tables and the formatted workbook are preserved.
The new generator was replayed from a freshly reconstructed native/experimental
merge, and its output matches the packaged preparation tables.

INTERPRETATION FOR COMPARING SPECIES
Mouse and human use McShane's fitted-model ED/NED classifications. The yeast
analysis uses Quyen's operational R2 rule, which measures poor single-exponential
fit and can reflect noise or model mismatch. Report yeast separately rather
than pooling it as the same biological NED endpoint.
Mouse has a positive point estimate but its CI includes 1. Report that no
statistically significant association was detected, not that there is no effect.
The species use different structural confidence sources and keyword snapshots.
In the human workflow, the membrane export filename is dated 2025-03-19 and
its secreted export filename is dated 2026-09-21. They are not both March 2025
exports. Export dates and UniProt release versions are separate provenance.

Revision verification: all 13 scientific CSVs and the workbook are unchanged.
Independent refit: 88 numerical comparisons passed; no failures.
