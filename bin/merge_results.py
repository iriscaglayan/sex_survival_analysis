#!/usr/bin/env python

import argparse
from pathlib import Path

import pandas as pd
from statsmodels.stats.multitest import multipletests


def stack(paths):
    return pd.concat([pd.read_csv(p) for p in paths], ignore_index=True) \
             .sort_values("cancer", kind="stable")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--representation", nargs="+", required=True)
    p.add_argument("--dropped", nargs="+", required=True)
    p.add_argument("--cox", nargs="+", required=True)
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    rep = stack(a.representation)
    drop = stack(a.dropped)
    cox = stack(a.cox)
    rep.to_csv(out / "all_representation.csv", index=False)
    drop.to_csv(out / "all_dropped.csv", index=False)
    cox.to_csv(out / "all_cox.csv", index=False)

    # headline table: one row per cancer
    sex = cox[(cox["model"] == "all") & (cox["term"] == "female")]
    sex = sex[["cancer", "HR", "CI_low", "CI_high", "p", "n", "n_events",
               "stage_adjusted"]]
    pct_f = rep[rep["sex"] == "female"][["cancer", "pct_of_cohort"]] \
              .rename(columns={"pct_of_cohort": "pct_female"})
    deaths = drop.pivot(index="cancer", columns="sex", values="n_deaths_after") \
                 .add_prefix("deaths_").reset_index()
    headline = sex.merge(pct_f, on="cancer").merge(deaths, on="cancer")
    # many cancers tested at once -> Benjamini-Hochberg false discovery rate
    ok = headline["p"].notna()
    headline["q_fdr"] = pd.NA
    if ok.any():
        headline.loc[ok, "q_fdr"] = multipletests(headline.loc[ok, "p"],
                                                  method="fdr_bh")[1]
    headline.to_csv(out / "sex_effect.csv", index=False)

    pd.set_option("display.float_format", lambda v: f"{v:.2f}")
    print(headline.to_string(index=False))


if __name__ == "__main__":
    main()
