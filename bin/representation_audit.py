#!/usr/bin/env python

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no screen needed (works inside Nextflow/Docker)
import matplotlib.pyplot as plt
import pandas as pd

SEX_COL = "gender.demographic"
AGE_COL = "age_at_index.demographic"
STAGE_COL = "ajcc_pathologic_stage.diagnoses"


def load_primary_tumours(clinical_path, survival_path):
    clin = pd.read_csv(clinical_path, sep="\t", low_memory=False)
    surv = pd.read_csv(survival_path, sep="\t")

    # first column of the clinical table is the sample ID
    clin = clin.rename(columns={clin.columns[0]: "sample"})

    # not every cohort has every column (e.g. lymphoma has no AJCC stage)
    if SEX_COL not in clin:
        raise SystemExit(f"[AUDIT] {clinical_path}: no '{SEX_COL}' column")
    if AGE_COL not in clin:
        alt = "age_at_diagnosis.diagnoses"          # in days
        clin[AGE_COL] = clin[alt] / 365.25 if alt in clin else pd.NA
    if STAGE_COL not in clin:
        print(f"note: no '{STAGE_COL}' column; stage treated as missing")
        clin[STAGE_COL] = pd.NA

    df = surv.merge(clin[["sample", SEX_COL, AGE_COL, STAGE_COL]],
                    on="sample", how="left")


    df = df[df["sample"].str[13:15].isin(["01", "03", "09"])]

    # one sample per patient (sort so the choice is reproducible)
    df = df.sort_values("sample").drop_duplicates("_PATIENT", keep="first")

    df = df.rename(columns={SEX_COL: "sex", AGE_COL: "age", STAGE_COL: "stage"})
    df["sex"] = df["sex"].str.lower()
    # treat "not reported" / "unknown" / blanks as missing
    df["stage"] = df["stage"].replace(
        {"Not Reported": pd.NA, "not reported": pd.NA, "Unknown": pd.NA, "": pd.NA})
    return df

def summarise(df, cancer):
    rows = []
    n_total = len(df)
    for sex in ["female", "male"]:
        d = df[df["sex"] == sex]
        n = len(d)
        rows.append({
            "cancer": cancer,
            "sex": sex,
            "n_patients": n,
            "pct_of_cohort": round(100 * n / n_total, 1) if n_total else 0,
            "n_deaths": int(d["OS"].sum()),          # events drive Cox power
            "median_age": d["age"].median(),
            "pct_missing_stage": round(100 * d["stage"].isna().mean(), 1) if n else 0,
            "pct_missing_age": round(100 * d["age"].isna().mean(), 1) if n else 0,
            "pct_missing_survival": round(
                100 * d[["OS", "OS.time"]].isna().any(axis=1).mean(), 1) if n else 0,
        })
    n_unknown_sex = int((~df["sex"].isin(["female", "male"])).sum())
    return pd.DataFrame(rows), n_unknown_sex


def plot(summary, cancer, out_png):
    fig, axes = plt.subplots(1, 2, figsize=(8, 3.5))
    colors = {"female": "#C2185B", "male": "#1565C0"}
    for ax, col, title in [(axes[0], "n_patients", "Patients"),
                           (axes[1], "n_deaths", "Deaths (events)")]:
        ax.bar(summary["sex"], summary[col],
               color=[colors[s] for s in summary["sex"]])
        for i, v in enumerate(summary[col]):
            ax.text(i, v, str(v), ha="center", va="bottom")
        ax.set_title(title)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(f"TCGA-{cancer}: who is in the data?")
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cancer", required=True, help="TCGA code, e.g. STAD")
    p.add_argument("--clinical", required=True)
    p.add_argument("--survival", required=True)
    p.add_argument("--outdir", default=".")
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    df = load_primary_tumours(a.clinical, a.survival)
    summary, n_unknown = summarise(df, a.cancer)

    summary.to_csv(out / f"{a.cancer}_representation.csv", index=False)
    df.to_csv(out / f"{a.cancer}_patients.csv", index=False)  # used by later steps
    plot(summary, a.cancer, out / f"{a.cancer}_representation.png")

    print(f"TCGA-{a.cancer}: {len(df)} patients (primary tumour, 1 per patient); "
          f"{n_unknown} with sex not recorded")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
