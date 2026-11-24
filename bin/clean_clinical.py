#!/usr/bin/env python
"""
CLEAN_CLINICAL
Turn the audit's patient table into a model-ready table for one cancer.

Steps:
  1. group detailed AJCC stages (Stage IA, Stage IIIB, ...) into I, II, III, IV
  2. convert survival time from days to years
  3. drop patients missing sex, age, stage or survival (and say how many, per sex)
  4. sanity checks: stop with a clear error if the cleaned data look wrong

Cancers without AJCC staging (e.g. leukaemia, brain tumours) are handled by
--min-staged: if fewer than this fraction of patients have a stage, stage is
not required and the Cox model adjusts for age only (stage_group = "unstaged").

Input : <CANCER>_patients.csv from representation_audit.py
Output: <CANCER>_clean.csv and <CANCER>_dropped.csv (how many were lost, per sex)

Usage:
  clean_clinical.py --cancer STAD --patients results/audit/STAD_patients.csv \
      --outdir results/clean
"""
import argparse
import re
from pathlib import Path

import pandas as pd


def group_stage(stage):
    """'Stage IIIB' -> 'III'; anything unparseable -> missing."""
    if pd.isna(stage):
        return pd.NA
    m = re.match(r"^\s*stage\s+(IV|III|II|I)", str(stage), flags=re.IGNORECASE)
    return m.group(1).upper() if m else pd.NA


def sanity_check(clean, cancer):
    """Stop the pipeline loudly instead of modelling broken data."""
    problems = []
    if clean.empty:
        problems.append("no patients left after cleaning")
    if clean["_PATIENT"].duplicated().any():
        problems.append(f"{clean['_PATIENT'].duplicated().sum()} duplicated patients")
    # impossible ages stop the run; young patients are real (e.g. paediatric
    # brain tumours) so they only get a warning
    bad_age = clean[(clean["age"] < 0) | (clean["age"] > 110)]
    if len(bad_age):
        problems.append(f"{len(bad_age)} impossible ages (<0 or >110)")
    young = int((clean["age"] < 18).sum())
    if young:
        print(f"WARNING: TCGA-{cancer} includes {young} patient(s) under 18 "
              f"(ages {sorted(clean.loc[clean['age'] < 18, 'age'].tolist())})")
    if (clean["time_years"] > 30).any():
        problems.append(f"{(clean['time_years'] > 30).sum()} follow-up times over 30 years")
    if not set(clean["sex"]).issubset({"female", "male"}):
        problems.append(f"unexpected sex values: {set(clean['sex'])}")
    if problems:
        raise SystemExit(f"[CLEAN_CLINICAL] TCGA-{cancer} failed sanity checks: "
                         + "; ".join(problems))
    print(f"Sanity checks passed for TCGA-{cancer}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cancer", required=True)
    p.add_argument("--patients", required=True)
    p.add_argument("--outdir", default=".")
    p.add_argument("--min-staged", type=float, default=0.5,
                   help="below this fraction of staged patients, ignore stage")
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(a.patients)
    start = df["sex"].value_counts()

    df["stage_group"] = df["stage"].map(group_stage)
    frac_staged = df["stage_group"].notna().mean()
    use_stage = frac_staged >= a.min_staged
    if not use_stage:
        print(f"TCGA-{a.cancer}: only {frac_staged:.0%} of patients staged; "
              "modelling without stage")
        df["stage_group"] = "unstaged"
    df["time_years"] = df["OS.time"] / 365.25
    df["event"] = df["OS"]

    keep = (
        df["sex"].isin(["female", "male"])
        & df["age"].notna()
        & df["stage_group"].notna()
        & df["time_years"].notna() & (df["time_years"] > 0)
        & df["event"].isin([0, 1])
    )
    clean = df.loc[keep, ["sample", "_PATIENT", "sex", "age",
                          "stage_group", "time_years", "event"]].copy()
    clean["event"] = clean["event"].astype(int)
    clean.insert(0, "cancer", a.cancer)

    end = clean["sex"].value_counts()
    dropped = pd.DataFrame({
        "cancer": a.cancer,
        "sex": ["female", "male"],
        "n_before": [int(start.get(s, 0)) for s in ["female", "male"]],
        "n_after": [int(end.get(s, 0)) for s in ["female", "male"]],
    })
    dropped["n_dropped"] = dropped["n_before"] - dropped["n_after"]
    dropped["n_deaths_after"] = [int(clean.loc[clean["sex"] == s, "event"].sum())
                                 for s in ["female", "male"]]

    sanity_check(clean, a.cancer)
    clean.to_csv(out / f"{a.cancer}_clean.csv", index=False)
    dropped.to_csv(out / f"{a.cancer}_dropped.csv", index=False)

    print(f"TCGA-{a.cancer}: {len(clean)} patients kept")
    print("Stage groups:", clean["stage_group"].value_counts().sort_index().to_dict())
    print(dropped.to_string(index=False))


if __name__ == "__main__":
    main()
