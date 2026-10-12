#!/usr/bin/env python3
"""Rebuild mouse confidence candidates and AlphaLasso evidence for native models.

Run after prepare_mouse.py. No network access, pLDDT cutoff, or fitted outcome
is used here: topology and keyword filters define candidates for subsequent
exact-v6 residue-confidence retrieval. The analysis applies pLDDT >=70 later.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import zipfile
from pathlib import Path

import pandas as pd


def fingerprint(path: Path) -> dict:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return {"filename": path.name, "bytes": path.stat().st_size, "sha256": digest.hexdigest()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", required=True, type=Path)
    parser.add_argument("--native-annotations", required=True, type=Path)
    parser.add_argument("--knots", required=True, type=Path)
    parser.add_argument("--alphalasso", required=True, type=Path)
    parser.add_argument("--keywords", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists() and any(args.out.iterdir()):
        raise ValueError("Use a new or empty output directory")
    args.out.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(args.prepared, dtype={"UniProt_ID": str})
    required = {"UniProt_ID", "Model_ID", "Model_version", "NCLE_any", "Covalent_native_records"}
    if not required.issubset(data.columns):
        raise ValueError("Unexpected prepared mouse schema")
    if data.UniProt_ID.isna().any() or not data.UniProt_ID.is_unique:
        raise ValueError("Prepared accessions must be unique and nonmissing")
    matched = data.loc[data.Model_ID.notna()]
    if not matched.Model_version.eq(6).all():
        raise ValueError("Confidence candidates require the supplied model version 6")
    if not matched.Model_ID.eq("AF-" + matched.UniProt_ID + "-F1-model_v6").all():
        raise ValueError("Prepared model ID does not match its accession and F1/version 6")
    if not matched.Covalent_native_records.isin([0, 1]).all():
        raise ValueError("Matched native covalent flags must be binary")
    if not data.NCLE_any.dropna().isin([0, 1]).all():
        raise ValueError("Nonmissing native NCLE calls must be binary")
    native = pd.read_csv(args.native_annotations, dtype={"UniProt_ID": str})
    if not {"UniProt_ID", "Model_ID"}.issubset(native.columns):
        raise ValueError("Unexpected native annotation schema")
    if native.UniProt_ID.isna().any() or not native.UniProt_ID.is_unique:
        raise ValueError("Native accessions must be unique and nonmissing")
    native_ids = set(native.UniProt_ID)
    expected_models = matched.merge(native[["UniProt_ID", "Model_ID"]],
                                   on="UniProt_ID", how="left", validate="one_to_one",
                                   suffixes=("_prepared", "_native"))
    if not expected_models.Model_ID_prepared.eq(expected_models.Model_ID_native).all():
        raise ValueError("Prepared native matches disagree with native annotation models")

    knots = pd.read_csv(args.knots, sep=";", usecols=[0], dtype=str)
    knot_ids = set(knots.iloc[:, 0].dropna().str.strip())
    evidence = []
    with zipfile.ZipFile(args.alphalasso) as archive:
        for name in archive.namelist():
            if not name.endswith(".tsv") or "__MACOSX" in name:
                continue
            with archive.open(name) as handle:
                for chunk in pd.read_csv(handle, sep="\t", dtype=str, chunksize=100000):
                    if not {"Organism", "Uniprot"}.issubset(chunk.columns):
                        raise ValueError(f"Unexpected AlphaLasso schema in {name}")
                    rows = chunk.loc[chunk.Organism.str.strip().eq("Mouse")].copy()
                    if rows.empty:
                        continue
                    rows["UniProt_ID"] = rows.Uniprot.str.strip().str.replace(r"-F\d+$", "", regex=True)
                    rows["AlphaLasso_table"] = name
                    evidence.append(rows)
    evidence = pd.concat(evidence, ignore_index=True) if evidence else pd.DataFrame(columns=["UniProt_ID"])
    lasso_ids = set(evidence.UniProt_ID.dropna())

    keyword = pd.read_csv(args.keywords, sep="\t", dtype=str)
    if not {"Entry", "Organism (ID)", "Keyword ID"}.issubset(keyword.columns):
        raise ValueError("Unexpected frozen keyword snapshot schema")
    if keyword.Entry.isna().any() or not keyword.Entry.is_unique or not keyword["Organism (ID)"].eq("10090").all():
        raise ValueError("Keyword accessions must be unique mouse entries")
    kw = keyword.set_index("Entry")["Keyword ID"].fillna("").str.split(";").map(
        lambda values: {value.strip() for value in values})
    covered = data.UniProt_ID.isin(kw.index)
    excluded_keywords = data.UniProt_ID.map(kw).map(
        lambda values: isinstance(values, set) and bool({"KW-0472", "KW-0964"} & values))

    masks = [
        ("Exact native accession match", data.Model_ID.notna()),
        ("Binary native NCLE call", data.NCLE_any.notna()),
        ("Exclude supplied knots", ~data.UniProt_ID.isin(knot_ids)),
        ("Exclude native and AlphaLasso covalent records",
         data.Covalent_native_records.eq(0) & ~data.UniProt_ID.isin(lasso_ids)),
        ("Require frozen keyword coverage", covered),
        ("Exclude Membrane and Secreted keywords", ~excluded_keywords),
    ]
    keep = pd.Series(True, index=data.index)
    flow = [{"Step": "All prepared experimental records", "Remaining": len(data), "Removed": 0}]
    for step, condition in masks:
        before = int(keep.sum())
        keep &= condition.fillna(False)
        flow.append({"Step": step, "Remaining": int(keep.sum()), "Removed": before - int(keep.sum())})
    candidates = data.loc[keep, ["UniProt_ID", "Model_ID", "Model_version"]].sort_values("UniProt_ID")
    candidates.to_csv(args.out / "confidence_candidate_ids.csv", index=False)
    # Reproduce the historical evidence scope: Mouse-labelled archive rows
    # whose exact accession occurs in the FULL supplied native model archive.
    native_evidence = evidence.loc[evidence.UniProt_ID.isin(native_ids)].copy()
    native_evidence.to_csv(args.out / "alphalasso_mouse_exact_evidence.csv", index=False)
    summary = {
        "Confidence_candidate_count": len(candidates),
        "Mouse_AlphaLasso_rows": len(evidence),
        "Mouse_AlphaLasso_unique_accessions": len(lasso_ids),
        "Native_overlap_AlphaLasso_rows": len(native_evidence),
        "Native_overlap_AlphaLasso_unique_accessions": native_evidence.UniProt_ID.nunique(),
        "AlphaLasso_output_scope": "Mouse-labelled archive rows with an exact accession in the full supplied native annotation archive",
        "Candidate_output_order": "Ascending UniProt_ID; row order does not enter the analysis",
        "Candidate_definition": "Exact native match, binary NCLE, no knot/covalent flag, keyword coverage, no Membrane or Secreted keyword",
        "Protein_pLDDT_filter_applied": False,
        "Degradation_category_filter_applied": False,
        "Network_access_used": False,
        "Flow": flow,
        "Sources": [fingerprint(path) for path in [args.prepared, args.native_annotations, args.knots, args.alphalasso, args.keywords]],
    }
    (args.out / "confidence_preparation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
