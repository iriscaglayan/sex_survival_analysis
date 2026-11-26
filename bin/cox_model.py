#!/usr/bin/env python
import argparse
from pathlib import Path
import pandas as pd
from lifelines import CoxPHFitter


def design_matrix(d, include_sex):
    X = pd.DataFrame({
        "time_years": d["time_years"],
        "event": d["event"],
        "age_per10y": d["age"] / 10,
    })
    if include_sex:
        X["female"] = (d["sex"] == "female").astype(int)
    stages = pd.get_dummies(d["stage_group"], prefix="stage", dtype=int)
    stages = stages.drop(columns=["stage_unstaged"], errors="ignore")
    # reference = lowest stage present; order I < II < III < IV
    order = [c for c in ["stage_I", "stage_II", "stage_III", "stage_IV"] if c in stages]
    X = pd.concat([X, stages[order[1:]]], axis=1)
    # drop terms with no variation (e.g. a stage nobody in this subgroup has)
    const = [c for c in X.columns if c not in ("time_years", "event") and X[c].nunique() < 2]
    return X.drop(columns=const)


def fit(d, model, include_sex, cancer, min_events):
    n, n_events = len(d), int(d["event"].sum())
    base = {"cancer": cancer, "model": model, "n": n, "n_events": n_events,
            "low_events_flag": n_events < min_events,
            "stage_adjusted": bool((d["stage_group"] != "unstaged").any())}
    if n_events < 5:
        return [{**base, "term": None, "note": "too few deaths to fit"}]

    X = design_matrix(d, include_sex)
    cph = CoxPHFitter()
    try:
        cph.fit(X, duration_col="time_years", event_col="event")
    except Exception as e:  # convergence trouble: refit with a tiny penalty
        cph = CoxPHFitter(penalizer=0.01)
        cph.fit(X, duration_col="time_years", event_col="event")
        base["note"] = f"refit with penalizer=0.01 ({type(e).__name__})"

    s = cph.summary
    rows = []
    for term in s.index:
        rows.append({**base,
                     "term": term,
                     "HR": s.loc[term, "exp(coef)"],
                     "CI_low": s.loc[term, "exp(coef) lower 95%"],
                     "CI_high": s.loc[term, "exp(coef) upper 95%"],
                     "p": s.loc[term, "p"],
                     "concordance": cph.concordance_index_})
    return rows


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cancer", required=True)
    p.add_argument("--clean", required=True)
    p.add_argument("--outdir", default=".")
    p.add_argument("--min-events", type=int, default=40)
    a = p.parse_args()

    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(a.clean)

    rows = []
    rows += fit(d, "all", True, a.cancer, a.min_events)
    rows += fit(d[d["sex"] == "female"], "female", False, a.cancer, a.min_events)
    rows += fit(d[d["sex"] == "male"], "male", False, a.cancer, a.min_events)

    res = pd.DataFrame(rows)
    res.to_csv(out / f"{a.cancer}_cox.csv", index=False)

    pd.set_option("display.float_format", lambda v: f"{v:.3f}")
    print(res[["model", "term", "HR", "CI_low", "CI_high", "p",
               "n", "n_events", "low_events_flag"]].to_string(index=False))

    sex = res[(res["model"] == "all") & (res["term"] == "female")]
    if not sex.empty:
        r = sex.iloc[0]
        print(f"\nTCGA-{a.cancer}: women vs men, adjusted for age and stage: "
              f"HR = {r.HR:.2f} (95% CI {r.CI_low:.2f}-{r.CI_high:.2f}), p = {r.p:.3f}")


if __name__ == "__main__":
    main()
