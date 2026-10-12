#!/usr/bin/env python3
"""Crude NCLE–NED association and length-adjusted logistic regression.

Input: a curated CSV/TSV, one protein per row, with columns
UniProt_ID, NCLE, Degradation_profile, Protein_length.
No ED/NED labels are inferred from rates, half-lives, or R-squared.
Only Model 1 is fitted. Complex membership and disorder are never used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy import optimize, special, stats


REQUIRED = ["UniProt_ID", "NCLE", "Degradation_profile", "Protein_length"]


def validate_input(path: Path) -> pd.DataFrame:
    if path.suffix.lower() not in {".csv", ".tsv", ".txt"}:
        raise ValueError("Provide a curated CSV/TSV. A raw workbook needs source-specific reconstruction first.")
    data = pd.read_csv(path, sep="\t" if path.suffix.lower() != ".csv" else ",", dtype=str)
    data.columns = data.columns.str.strip()
    if data.columns.duplicated().any():
        raise ValueError("Duplicate column names after whitespace normalization")
    if not set(REQUIRED).issubset(data.columns):
        raise ValueError(f"Missing columns: {sorted(set(REQUIRED) - set(data.columns))}")
    if data.empty:
        raise ValueError("The input has headers but no protein records. An empty template cannot be analyzed.")
    data = data.copy()
    for column in ["UniProt_ID", "Degradation_profile"]:
        data[column] = data[column].astype("string").str.strip()
        data[column] = data[column].replace("", pd.NA)
    if data["UniProt_ID"].isna().any():
        raise ValueError("Every input row must have a UniProt ID")
    if data["UniProt_ID"].duplicated().any():
        ids = data.loc[data["UniProt_ID"].duplicated(False), "UniProt_ID"].unique().tolist()
        raise ValueError(f"Duplicate IDs must be resolved before analysis: {ids[:20]}")
    if data["UniProt_ID"].str.contains(r"[;|]", regex=True).any():
        raise ValueError("Grouped accessions/FASTA headers need a documented mapping before analysis")
    data["Degradation_profile"] = data["Degradation_profile"].str.upper()
    invalid_profiles = set(data["Degradation_profile"].dropna().unique()) - {"ED", "NED", "UN"}
    if invalid_profiles:
        raise ValueError(f"Unexpected degradation labels: {sorted(invalid_profiles)}")
    for column in ["NCLE", "Protein_length"]:
        raw = data[column].astype("string").str.strip().replace("", pd.NA)
        numeric = pd.to_numeric(raw, errors="coerce")
        invalid = raw.notna() & (numeric.isna() | ~np.isfinite(numeric))
        if invalid.any():
            raise ValueError(f"Non-numeric/non-finite {column} values: {raw[invalid].unique().tolist()[:10]}")
        data[column] = numeric.astype(float)
    if not set(data["NCLE"].dropna().unique()).issubset({0.0, 1.0}):
        raise ValueError("NCLE must be 1 (native non-covalent entanglement), 0 (confirmed no detected NCLE), or missing")
    lengths = data["Protein_length"].dropna()
    if (lengths <= 0).any() or not np.equal(lengths, np.floor(lengths)).all():
        raise ValueError("Nonmissing Protein_length must be a positive integer residue count")
    return data


def crude(cohort: pd.DataFrame) -> tuple[dict, list[dict]]:
    a = int(((cohort["NCLE"] == 1) & (cohort["NED_binary"] == 1)).sum())
    b = int(((cohort["NCLE"] == 1) & (cohort["NED_binary"] == 0)).sum())
    c = int(((cohort["NCLE"] == 0) & (cohort["NED_binary"] == 1)).sum())
    d = int(((cohort["NCLE"] == 0) & (cohort["NED_binary"] == 0)).sum())
    table = np.array([[a, b], [c, d]], dtype=int)
    rows = [
        {"NCLE": 1, "NED": a, "Non_NED": b, "Total": a + b},
        {"NCLE": 0, "NED": c, "Non_NED": d, "Total": c + d},
    ]
    if len(cohort) == 0 or (table.sum(axis=0) == 0).any() or (table.sum(axis=1) == 0).any():
        return {"N": len(cohort), "Status": "not_estimable_missing_comparison_group"}, rows
    estimate, p = stats.fisher_exact(table, alternative="two-sided")
    if np.all(table > 0):
        se = np.sqrt(np.sum(1.0 / table))
        lo, hi = np.exp(np.log(estimate) + stats.norm.ppf(0.975) * np.array([-se, se]))
        ci_method = "Woolf log-OR Wald"
    else:
        interval = stats.contingency.odds_ratio(table, kind="conditional").confidence_interval(0.95)
        lo, hi = interval.low, interval.high
        ci_method = "Exact conditional interval (uncorrected sample OR reported)"
    return {
        "N": len(cohort), "OR": float(estimate), "CI_lower": float(lo), "CI_upper": float(hi),
        "P_value": float(p), "P_method": "Two-sided Fisher exact", "CI_method": ci_method,
        "Status": "ok" if np.all(table > 0) else "zero_cell",
    }, rows


def fit_model1(cohort: pd.DataFrame) -> tuple[list[dict], dict]:
    n = len(cohort)
    if n == 0 or cohort["NCLE"].nunique() != 2 or cohort["NED_binary"].nunique() != 2:
        return [], {"N": n, "Status": "not_estimable_missing_comparison_group"}
    mean = float(cohort["Protein_length"].mean())
    sd = float(cohort["Protein_length"].std(ddof=1))
    if not np.isfinite(sd) or sd <= 0:
        return [], {"N": n, "Status": "not_estimable_constant_length"}
    cohort["z_Protein_length"] = (cohort["Protein_length"] - mean) / sd
    x = np.column_stack([np.ones(n), cohort["NCLE"], cohort["z_Protein_length"]])
    y = cohort["NED_binary"].to_numpy(dtype=float)
    if np.linalg.matrix_rank(x) != 3:
        return [], {"N": n, "Status": "not_estimable_rank_deficient"}
    signed = (2 * y - 1)[:, None] * x
    separation = optimize.linprog(
        np.zeros(3), A_ub=np.vstack([-signed, -signed.sum(axis=0)]),
        b_ub=np.r_[np.zeros(n), -1.0], bounds=[(None, None)] * 3, method="highs",
    )
    if separation.status == 0:
        return [], {"N": n, "Status": "not_estimable_complete_or_quasi_separation"}
    if separation.status != 2:
        return [], {"N": n, "Status": "separation_diagnostic_failed", "Message": separation.message}

    def loss(beta):
        eta = x @ beta
        return float(np.sum(np.logaddexp(0, eta) - y * eta))

    def gradient(beta):
        return x.T @ (special.expit(x @ beta) - y)

    def information(beta):
        p = special.expit(x @ beta)
        return x.T @ ((p * (1 - p))[:, None] * x)

    initial = np.array([special.logit(y.mean()), 0.0, 0.0])
    fitted = optimize.minimize(loss, initial, method="trust-exact", jac=gradient, hess=information,
                               options={"gtol": 1e-8, "maxiter": 1000})
    score_max = float(np.max(np.abs(gradient(fitted.x))))
    if not np.isfinite(fitted.fun) or score_max > 1e-6 or np.max(np.abs(fitted.x)) > 30:
        return [], {"N": n, "Status": "logistic_fit_not_reliable", "Message": str(fitted.message),
                    "Max_abs_score": score_max}
    info = information(fitted.x)
    covariance = np.linalg.inv(info)
    errors = np.sqrt(np.diag(covariance))
    if not np.isfinite(errors).all():
        return [], {"N": n, "Status": "nonfinite_standard_errors"}
    rows = []
    for term, beta, se in zip(["Intercept", "NCLE", "z_Protein_length"], fitted.x, errors):
        z = beta / se
        radius = stats.norm.ppf(0.975) * se
        rows.append({"Term": term, "N": n, "Coefficient": float(beta), "SE": float(se),
                     "OR": float(np.exp(beta)), "CI_lower": float(np.exp(beta - radius)),
                     "CI_upper": float(np.exp(beta + radius)), "P_value": float(2 * stats.norm.sf(abs(z))),
                     "P_method": "Two-sided Wald z", "CI_method": "Logistic Wald", "Status": "ok"})
    return rows, {"N": n, "Status": "ok", "Mean_length": mean, "SD_length_ddof1": sd,
                  "Log_likelihood": -float(fitted.fun), "Max_abs_score": score_max,
                  "Information_condition_number": float(np.linalg.cond(info)), "Iterations": int(fitted.nit)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, metavar="SPECIES=CURATED_CSV")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--include-un-sensitivity", action="store_true",
                        help="Optional: also compare NED with ED+UN. UN is not asserted to be ED.")
    args = parser.parse_args()
    if args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("Output directory must be empty or new, to prevent overwriting prior results")
    all_results, all_counts, audits = [], [], []
    export_cohorts = []
    names = set()
    for argument in args.input:
        if "=" not in argument:
            raise ValueError("Each --input must have the form SPECIES=PATH")
        species, value = argument.split("=", 1)
        if not species or any(c not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for c in species):
            raise ValueError("Use a species label containing letters, digits, hyphens, or underscores")
        if species in names:
            raise ValueError(f"Duplicate species label: {species}")
        names.add(species)
        path = Path(value)
        data = validate_input(path)
        modes = [("Primary_NED_vs_ED", {"ED", "NED"})]
        if args.include_un_sensitivity:
            modes.append(("Sensitivity_NED_vs_ED_UN", {"ED", "NED", "UN"}))
        for mode, labels in modes:
            selected = data.loc[data["Degradation_profile"].isin(labels)].copy()
            selected = selected.loc[selected["NCLE"].notna()].copy()
            selected["NED_binary"] = selected["Degradation_profile"].eq("NED").astype(int)
            model_cohort = selected.loc[selected["Protein_length"].notna()].copy()
            crude_result, count_rows = crude(selected)
            all_results.append({"Species": species, "Analysis": mode, "Method": "Crude",
                                "Term": "NCLE", **crude_result})
            all_counts.extend({"Species": species, "Analysis": mode, "Cohort": "Crude", **r} for r in count_rows)
            if len(model_cohort) != len(selected):
                same_result, same_counts = crude(model_cohort)
                all_results.append({"Species": species, "Analysis": mode, "Method": "Crude_Model1_cohort",
                                    "Term": "NCLE", **same_result})
                all_counts.extend({"Species": species, "Analysis": mode, "Cohort": "Model1", **r} for r in same_counts)
            model_rows, diagnostic = fit_model1(model_cohort)
            if model_rows:
                all_results.extend({"Species": species, "Analysis": mode, "Method": "Model1", **r} for r in model_rows)
            else:
                all_results.append({"Species": species, "Analysis": mode, "Method": "Model1",
                                    "Term": "NCLE", **diagnostic})
            audits.append({"Species": species, "Analysis": mode, "Input_file": path.name,
                           "Input_SHA256": hashlib.sha256(path.read_bytes()).hexdigest(),
                           "Input_rows": len(data), "Input_profiles": data["Degradation_profile"].value_counts(dropna=False).to_dict(),
                           "Selected_profile_rows": int(data["Degradation_profile"].isin(labels).sum()),
                           "Excluded_missing_NCLE": int(data.loc[data["Degradation_profile"].isin(labels), "NCLE"].isna().sum()),
                           "Crude_N": len(selected), "Excluded_missing_length_from_Model1": len(selected) - len(model_cohort),
                           "Model1_diagnostic": diagnostic})
            model_cohort.insert(0, "Species", species)
            export_cohorts.append((f"{species}_{mode}_cohort.csv", model_cohort))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_results).to_csv(args.output_dir / "association_results.csv", index=False)
    pd.DataFrame(all_counts).to_csv(args.output_dir / "contingency_counts.csv", index=False)
    for filename, cohort in export_cohorts:
        cohort.to_csv(args.output_dir / filename, index=False)
    metadata = {"Model": "NED ~ NCLE + z(Protein_length)", "Length_standardization": "Within each species and model cohort, sample SD ddof=1",
                "Outcome": "NED=1, ED=0 (UN=0 only in explicitly requested sensitivity)",
                "Population": "Input must already contain the agreed eligible proteins. No additional pLDDT, knot, membrane, secreted, or disorder filtering occurs in this runner.",
                "Runtime": {"Python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
                "Audits": audits}
    (args.output_dir / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    view = pd.DataFrame(all_results)
    print(view.loc[view["Term"].eq("NCLE")].to_string(index=False))


if __name__ == "__main__":
    main()
