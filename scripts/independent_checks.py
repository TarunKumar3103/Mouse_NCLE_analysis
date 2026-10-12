#!/usr/bin/env python3
"""Independent checks of the human, mouse and yeast NCLE-NED packages.

Rebuilds every cohort from the raw source files without importing any package
code, refits every model with statsmodels, and compares against the numbers
each package reports. Also confirms that your fresh replays are byte-identical
to the shipped result tables.

Usage (from the folder that holds the three extracted packages):
    python independent_checks.py --human Human_NCLE_reproduction \
        --mouse Mouse_NCLE_analysis --yeast Yeast_NCLE_ED_NED_analysis
"""
import argparse
import filecmp
import io
import json
import sys
import warnings
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

warnings.filterwarnings("ignore")
RESULTS = []


def check(package, name, passed, detail=""):
    RESULTS.append((package, name, bool(passed)))
    print(f"  [{'PASS' if passed else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def crude(ncle, y):
    a = int(((ncle == 1) & (y == 1)).sum()); b = int(((ncle == 1) & (y == 0)).sum())
    c = int(((ncle == 0) & (y == 1)).sum()); d = int(((ncle == 0) & (y == 0)).sum())
    odds = a * d / (b * c); se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d); q = stats.norm.ppf(0.975)
    return dict(OR=odds, lo=odds * np.exp(-q * se), hi=odds * np.exp(q * se),
                p=stats.fisher_exact([[a, b], [c, d]])[1], table=(a, b, c, d))


def logit(df, y, terms):
    X = sm.add_constant(df[terms].astype(float))
    fit = sm.Logit(df[y].astype(float), X).fit(disp=0, method="newton", tol=1e-12, maxiter=200)
    ci = fit.conf_int()
    return {t: dict(OR=np.exp(fit.params[t]), lo=np.exp(ci.loc[t, 0]), hi=np.exp(ci.loc[t, 1]),
                    p=fit.pvalues[t]) for t in terms}, fit


def close(mine, theirs, ci_rtol=1e-6):
    """OR and P to 1e-6; CI to ci_rtol (human code uses 1.96, not 1.959964)."""
    return (np.isclose(mine["OR"], theirs[0], rtol=1e-6) and np.isclose(mine["p"], theirs[3], rtol=1e-5, atol=1e-10)
            and np.isclose(mine["lo"], theirs[1], rtol=ci_rtol) and np.isclose(mine["hi"], theirs[2], rtol=ci_rtol))


def fasta_lengths(path):
    lengths, current = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            current = line[1:].split("|")[1]; lengths[current] = 0
        elif current:
            lengths[current] += len(line)
    return lengths


def same_tables(package, shipped_dir, replay_dir, names):
    replay_dir = Path(replay_dir)
    if not replay_dir.is_dir():
        check(package, f"fresh replay folder exists: {replay_dir.name}", False, "run the package commands first")
        return
    byte_diff, problems, worst = [], [], 0.0
    for n in names:
        a_path, b_path = Path(shipped_dir) / n, replay_dir / n
        if filecmp.cmp(a_path, b_path, shallow=False):
            continue
        byte_diff.append(n)
        # Different numpy/pandas/BLAS builds can change the last printed digit
        # of a float. Accept that only if every value agrees to 1e-9 relative.
        a, b = pd.read_csv(a_path, low_memory=False), pd.read_csv(b_path, low_memory=False)
        if list(a.columns) != list(b.columns) or a.shape != b.shape:
            problems.append(f"{n}: columns/shape differ"); continue
        for col in a.columns:
            x, y = a[col], b[col]
            if pd.api.types.is_numeric_dtype(x) and pd.api.types.is_numeric_dtype(y):
                x, y = x.to_numpy(float), y.to_numpy(float)
                if not np.array_equal(np.isnan(x), np.isnan(y)):
                    problems.append(f"{n}:{col} missing values differ"); continue
                ok = ~np.isnan(x)
                if ok.any():
                    rel = np.abs(x[ok] - y[ok]) / np.maximum(np.abs(x[ok]), 1e-300)
                    worst = max(worst, float(np.max(np.where(np.abs(x[ok] - y[ok]) <= 1e-12, 0, rel))))
                    if not np.allclose(x[ok], y[ok], rtol=1e-9, atol=1e-12):
                        problems.append(f"{n}:{col} values differ")
            elif not x.astype(str).equals(y.astype(str)):
                problems.append(f"{n}:{col} text differs")
    if not byte_diff:
        check(package, f"fresh replay byte-identical to shipped ({len(names)} tables)", True)
    else:
        check(package, f"fresh replay matches shipped ({len(names)} tables; {len(byte_diff)} differ only in float "
              f"rounding, max relative diff {worst:.1e})", not problems, "; ".join(problems))


# ----------------------------------------------------------------------------- human
def human(root):
    print("\nHUMAN")
    rs = root / "raw_sources"; wb = root / "data/reference_human_analysis.xlsx"
    mc = pd.read_excel(rs / "McShane_human_mmc4.xlsx", sheet_name="FinalTable_RPE_")[
        ["Protein IDs (Uniprot)", "Degradation profile"]]
    mc.columns = ["id", "prof"]
    pl = pd.read_csv(rs / "Human_Avg_pLDDTs.csv")
    ent = pd.read_csv(rs / "Human_AF_combined_20250614.csv", sep="|", dtype=str, usecols=["gene", "CCBond"])
    cc = ent.CCBond.map({"True": True, "1.0": True, "False": False, "0.0": False})
    g = pd.DataFrame({"gene": ent.gene, "cov": cc.eq(True), "nc": cc.eq(False)}).groupby("gene").agg(
        cov=("cov", "any"), nc=("nc", "any"))
    g["NCLE"] = np.where(g["cov"], np.nan, np.where(g["nc"], 1.0, 0.0))
    prot = fasta_lengths(rs / "UP000005640_9606.fasta")
    mem = set(fasta_lengths(rs / "uniprotkb_human_keywords_membrane_2025_03_19.fasta"))
    kn = pd.read_csv(rs / "AF4_knots.dat", sep="\t"); kn.columns = [c.lstrip("#") for c in kn.columns]
    knots = set(kn.loc[kn.Organism.eq("Human"), "Uniprot"].str.replace(r"-F\d+$", "", regex=True))
    m = mc.merge(pl.rename(columns={"gene": "id"}), on="id").merge(g[["NCLE"]], left_on="id", right_index=True, how="left")
    m = m[m["<pLDDT>"] >= 70]
    reported = pd.read_excel(wb, sheet_name="Association Results")
    ann = pd.read_csv(root / "data/complex_annotations.csv")
    expected_flow = {"Primary": [1043, 1034, 1023, 1022, 725], "Sensitivity": [1836, 1810, 1793, 1792, 1230]}
    for name, profiles, label in [("Primary", ["ED", "NED"], "Primary: NED vs ED"),
                                  ("Sensitivity", ["ED", "NED", "UN"], "Sensitivity: NED vs ED + UN")]:
        d = m[m.prof.isin(profiles)]; flow = [len(d)]
        d = d[d.NCLE.notna()]; flow.append(len(d))
        d = d[~d.id.isin(knots)]; flow.append(len(d))
        d = d[d.id.isin(prot)]; flow.append(len(d))
        d = d[~d.id.isin(mem)]; flow.append(len(d))
        check("human", f"{name}: cohort flow rebuilt from raw = {expected_flow[name]}", flow == expected_flow[name], str(flow))
        final = pd.read_excel(wb, sheet_name=f"{name} Cohort")
        check("human", f"{name}: final {len(final)} proteins are a subset of the rebuilt {flow[-1]}",
              set(final.UniProt_ID) <= set(d.id), f"secreted step removed {flow[-1] - len(final)}")
        j = final.merge(d.assign(L=d.id.map(prot)), left_on="UniProt_ID", right_on="id")
        mism = int((j.NCLE_x != j.NCLE_y).sum() + (j.Protein_length != j.L).sum()
                   + (~np.isclose(j.Average_pLDDT, j["<pLDDT>"])).sum() + (j.Degradation_profile != j.prof).sum())
        check("human", f"{name}: NCLE, length, pLDDT, profile match raw sources", mism == 0, f"{mism} mismatches")
        c = final.drop(columns="Complex_member").merge(ann, on="UniProt_ID", how="left")
        c["z"] = (c.Protein_length - c.Protein_length.mean()) / c.Protein_length.std(ddof=1)
        rep = reported[reported.Analysis.eq(label)]

        def ref(method_prefix, term):
            r = rep[rep.Method.str.startswith(method_prefix) & rep.Term.eq(term)].iloc[0]
            return r.Odds_ratio, r.CI_lower, r.CI_upper, r.P_value
        cr = crude(c.NCLE, c.NED_binary)
        check("human", f"{name}: crude OR {cr['OR']:.3f}, Fisher P {cr['p']:.3f}", close(cr, ref("Crude", "NCLE"), ci_rtol=1e-4))
        m1, _ = logit(c, "NED_binary", ["NCLE", "z"])
        for t, rt in [("NCLE", "NCLE"), ("z", "z_Protein_length")]:
            check("human", f"{name}: Model 1 {rt} OR {m1[t]['OR']:.3f}", close(m1[t], ref("Model 1", rt), ci_rtol=1e-4))
        cc2 = c.dropna(subset=["Complex_member"])
        m2, f2 = logit(cc2, "NED_binary", ["NCLE", "Complex_member", "z"])
        for t, rt in [("NCLE", "NCLE"), ("Complex_member", "Complex_member"), ("z", "z_Protein_length")]:
            check("human", f"{name}: Model 2 {rt} OR {m2[t]['OR']:.3f} (N={len(cc2)})",
                  close(m2[t], ref("Model 2", rt), ci_rtol=1e-4))
        _, f0 = logit(cc2, "NED_binary", ["NCLE", "z"])
        lr_p = stats.chi2.sf(2 * (f2.llf - f0.llf), 1)
        print(f"  [INFO] {name}: complex term Wald P {m2['Complex_member']['p']:.3f} vs likelihood-ratio P {lr_p:.4f}")
    md = pd.read_csv(root / "data/processed_MobiDB_master.csv")
    ok = md[md.Disorder_query_status.eq("ok")].copy()
    ok["fasta"] = ok.UniProt_ID.map(prot)
    n_mis = int((ok.fasta.notna() & ok.MobiDB_length.ne(ok.fasta)).sum())
    print(f"  [INFO] MobiDB rows whose UniProt length differs from the proteome FASTA: {n_mis} (expected 33)")


# ----------------------------------------------------------------------------- yeast
def yeast(root):
    print("\nYEAST")
    rs = root / "raw_sources"
    src = pd.read_excel(rs / "Christiano_yeast_mmc2.xlsx", sheet_name=0)
    nat = pd.read_pickle(rs / "AF_Yeast_entanglement.pkl")
    kn = {l.split(";")[0].strip() for l in open(rs / "knot_all.txt")}
    las = set()
    with zipfile.ZipFile(rs / "AlphaLasso.zip") as z:
        for n in z.namelist():
            if n.endswith(".tsv") and "__MACOSX" not in n:
                las |= set(pd.read_csv(z.open(n), sep="\t", usecols=["Uniprot"], dtype=str)
                           .Uniprot.str.replace(r"-F\d+$", "", regex=True))
    kw = pd.read_csv(rs / "uniprot_S_cerevisiae_keyword_snapshot_2026_10_10.tsv", sep="\t", dtype=str
                     ).set_index("Entry")["Keyword ID"].fillna("")
    src["id"] = src.UNIPROT.astype(str).str.strip()
    d = src[~src.id.str.contains(";")]; flow = [len(src), len(d)]
    d = d.merge(nat, left_on="id", right_on="Uniprot"); flow.append(len(d))
    d = d[d.plddt >= 70]; flow.append(len(d))
    d = d[~d.id.isin(kn)]; flow.append(len(d))
    d = d[~d.id.isin(las)]; flow.append(len(d))
    flow.append(len(d))
    d = d[~d.id.map(kw).fillna("").str.contains("KW-0472")]; flow.append(len(d))
    d = d[~d.id.map(kw).fillna("").str.contains("KW-0964")]; flow.append(len(d))
    shipped_flow = pd.read_csv(root / "results/cohort_flow.csv").Remaining.tolist()
    check("yeast", f"cohort flow rebuilt from raw = {shipped_flow}", flow == shipped_flow, str(flow))
    coh = pd.read_csv(root / "results/eligible_yeast_cohort.csv")
    check("yeast", f"same {len(coh)} protein IDs as the package", set(coh.UniProt_ID) == set(d.id))
    d = d.assign(NCLE=d.entanglement_AF.eq("Yes").astype(int), R2=d["R2 (quality of curve fitting)"])
    d["z"] = (d.length_AF - d.length_AF.mean()) / d.length_AF.std(ddof=1)
    rep = pd.read_csv(root / "results/association_results.csv")
    for cut in [0.90, 0.75, 0.85, 0.95]:
        d["y"] = (d.R2 < cut).astype(int)
        r = rep[np.isclose(rep.R2_cutoff, cut)]
        cr = crude(d.NCLE, d.y); row = r[r.Method.eq("Crude")].iloc[0]
        check("yeast", f"R2<{cut:.2f}: crude OR {cr['OR']:.3f}, Fisher P {cr['p']:.3f} ({int(d.y.sum())} NED)",
              close(cr, (row.OR, row.CI_lower, row.CI_upper, row.P_value)))
        m1, _ = logit(d, "y", ["NCLE", "z"])
        for t, rt in [("NCLE", "NCLE"), ("z", "z_Protein_length")]:
            row = r[r.Method.eq("Model 1") & r.Term.eq(rt)].iloc[0]
            check("yeast", f"R2<{cut:.2f}: Model 1 {rt} OR {m1[t]['OR']:.3f}",
                  close(m1[t], (row.OR, row.CI_lower, row.CI_upper, row.P_value)))
    same_tables("yeast", root / "results", root / "results/local_reproduction",
                ["association_results.csv", "contingency_counts.csv", "cohort_flow.csv", "eligible_yeast_cohort.csv",
                 "all_turnover_record_audit.csv", "model_diagnostics.csv", "covalent_lasso_evidence.csv",
                 "grouped_accession_review.csv", "native_annotation_export.csv"])


# ----------------------------------------------------------------------------- mouse
def mouse(root):
    print("\nMOUSE")
    rs = root / "raw_sources"
    mc = pd.read_excel(rs / "McShane_mouse_mmc1.xlsx", sheet_name="FinalTable_PostNorm_3T3_")[
        ["Protein IDs (Uniprot)", "Degradation profile"]]
    mc.columns = ["id", "prof"]
    native = {}
    with zipfile.ZipFile(rs / "Mouse_Native_HQ_GE.zip") as z:
        for n in z.namelist():
            if n.endswith("/"):
                continue
            acc = Path(n).name.split("-F1-")[0][3:]
            df = pd.read_csv(io.BytesIO(z.read(n)), sep="|", dtype=str)
            cov = (df.CCbond == "True").any(); nc = ((df.CCbond == "False") & (df.ENT == "True")).any()
            native[acc] = np.nan if cov else float(nc)
    kn = {l.split(";")[0].strip() for l in open(rs / "knot_all.txt")}
    las = set()
    with zipfile.ZipFile(rs / "AlphaLasso.zip") as z:
        for n in z.namelist():
            if n.endswith(".tsv") and "__MACOSX" not in n:
                df = pd.read_csv(z.open(n), sep="\t", usecols=["Uniprot", "Organism"], dtype=str)
                las |= set(df.Uniprot[df.Organism.str.strip().eq("Mouse")].str.replace(r"-F\d+$", "", regex=True))
    kw = pd.read_csv(rs / "uniprot_mouse_keyword_snapshot_2026_10_10.tsv", sep="\t", dtype=str
                     ).set_index("Entry")["Keyword ID"].fillna("")
    conf = {}
    for f in (rs / "raw_confidence_json").glob("*_confidence_v6.json"):
        scores = json.loads(f.read_text())["confidenceScore"]
        conf[f.name.split("_")[0]] = (float(np.mean(scores)), len(scores))
    mc["NCLE"] = mc.id.map(native)
    reported = pd.read_csv(root / "results/association_results.csv")
    shipped_flow = pd.read_csv(root / "results/cohort_flow.csv")
    for name, profiles, label in [("Primary", ["ED", "NED"], "Primary: NED vs ED"),
                                  ("Sensitivity", ["ED", "NED", "UN"], "Sensitivity: NED vs ED + UN")]:
        d = mc[mc.prof.isin(profiles)]; flow = [len(d)]
        d = d[d.id.isin(native)]; flow.append(len(d))
        d = d[d.NCLE.notna()]; flow.append(len(d))
        d = d[~d.id.isin(kn)]; flow.append(len(d))
        d = d[~d.id.isin(las)]; flow.append(len(d)); flow.append(len(d))
        d = d[~d.id.map(kw).str.contains("KW-0472")]; flow.append(len(d))
        d = d[~d.id.map(kw).str.contains("KW-0964")]; flow.append(len(d)); flow.append(len(d))
        d = d.assign(pl=d.id.map(lambda a: conf[a][0]), L=d.id.map(lambda a: conf[a][1]))
        d = d[d.pl >= 70]; flow.append(len(d))
        expected = shipped_flow[shipped_flow.Analysis.eq(label)].Total.tolist()
        check("mouse", f"{name}: cohort flow rebuilt from raw = {expected}", flow == expected, str(flow))
        ship = pd.read_csv(root / f"results/{name.lower()}_any_ncle_cohort.csv")
        check("mouse", f"{name}: same {len(ship)} protein IDs as the package", set(ship.UniProt_ID) == set(d.id))
        d["y"] = d.prof.eq("NED").astype(int)
        d["z"] = (d.L - d.L.mean()) / d.L.std(ddof=1)
        r = reported[reported.Analysis.eq(label) & reported.NCLE_definition.eq("Any detected NCLE")]
        cr = crude(d.NCLE, d.y); row = r[r.Method.str.startswith("Crude")].iloc[0]
        check("mouse", f"{name}: crude OR {cr['OR']:.3f}, Fisher P {cr['p']:.3f}",
              close(cr, (row.OR, row.CI_lower, row.CI_upper, row.P_value)))
        m1, _ = logit(d, "y", ["NCLE", "z"])
        for t, rt in [("NCLE", "NCLE"), ("z", "z_Protein_length")]:
            row = r[r.Method.str.startswith("Model 1") & r.Term.eq(rt)].iloc[0]
            check("mouse", f"{name}: Model 1 {rt} OR {m1[t]['OR']:.3f}",
                  close(m1[t], (row.OR, row.CI_lower, row.CI_upper, row.P_value)))
    same_tables("mouse", root / "results", root / "results/local_replay",
                ["association_results.csv", "contingency_counts.csv", "cohort_flow.csv", "model_diagnostics.csv",
                 "primary_any_ncle_cohort.csv", "primary_high_only_cohort.csv", "sensitivity_any_ncle_cohort.csv",
                 "sensitivity_high_only_cohort.csv", "all_experimental_record_audit.csv",
                 "covalent_lasso_evidence.csv", "local_quality_definition_differences.csv",
                 "sequence_length_version_differences.csv", "pLDDT_summary_comparison.csv"])
    same_tables("mouse", root / "preparation", root / "results/local_confidence_preparation",
                ["confidence_candidate_ids.csv", "alphalasso_mouse_exact_evidence.csv"])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--human", type=Path); p.add_argument("--mouse", type=Path); p.add_argument("--yeast", type=Path)
    a = p.parse_args()
    if not (a.human or a.mouse or a.yeast):
        p.error("give at least one of --human, --mouse, --yeast")
    for fn, path in [(human, a.human), (mouse, a.mouse), (yeast, a.yeast)]:
        if path:
            if not path.is_dir():
                sys.exit(f"Folder not found: {path}")
            fn(path)
    failed = [r for r in RESULTS if not r[2]]
    print(f"\n{len(RESULTS) - len(failed)}/{len(RESULTS)} independent checks passed.")
    print("ALL INDEPENDENT CHECKS PASSED" if not failed else "SOME CHECKS FAILED - see [FAIL] lines above")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
