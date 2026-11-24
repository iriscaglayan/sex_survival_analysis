#!/usr/bin/env python
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, NullFormatter, NullLocator
import pandas as pd
from lifelines import KaplanMeierFitter

COL = {"female": "#eb6834", "male": "#2a78d6"}
GREY = "#8a8985"
INK, INK2 = "#0b0b0b", "#52514e"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False,
    "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight",
})


def log_ticks(ax, ticks):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(FixedLocator(ticks))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())


def estimable(r):
    return np.isfinite(r["CI_low"]) and np.isfinite(r["CI_high"]) \
        and r["CI_low"] > 0 and r["CI_high"] / max(r["CI_low"], 1e-12) < 1e4


# figure 1
def fig_representation(rep, drop, min_events, out):
    order = rep[rep.sex == "female"].sort_values("pct_of_cohort")["cancer"].tolist()
    height = 0.6 * len(order) + 1.6
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, height),
                                 sharey=True, gridspec_kw={"wspace": 0.35})
    # keep the title close to the plots 
    fig.subplots_adjust(top=1 - 0.75 / height, bottom=0.9 / height)
    y = np.arange(len(order))

    pf = rep[rep.sex == "female"].set_index("cancer").loc[order, "pct_of_cohort"]
    a1.barh(y, pf, color=COL["female"], height=0.55)
    a1.barh(y, 100 - pf, left=pf, color=COL["male"], height=0.55)
    a1.axvline(50, color=INK2, lw=1, ls="--")
    for i, v in enumerate(pf):
        inside = v >= 28                     
        a1.text(v - 1.5 if inside else v + 1.5, i, f"{v:.0f}% women", va="center",
                ha="right" if inside else "left", color="white",
                fontsize=9, fontweight="bold")
    a1.set_yticks(y, order)
    a1.set_xlim(0, 100)
    a1.set_xlabel("% of patients")
    a1.set_title("Share of patients", loc="left", color=INK)

    d = drop.pivot(index="cancer", columns="sex", values="n_deaths_after").loc[order]
    h = 0.36
    for k, (sex, off) in enumerate([("female", h / 2), ("male", -h / 2)]):
        a2.barh(y + off, d[sex], height=h - 0.04, color=COL[sex], label=sex.capitalize())
        for i, v in enumerate(d[sex]):
            a2.text(v + 1, i + off, str(int(v)), va="center", fontsize=8, color=INK2)
    a2.axvline(min_events, color=INK2, lw=1, ls="--")
    a2.text(min_events, -0.75, f" {min_events}-death threshold",
            fontsize=8, color=INK2, va="center")
    a2.set_ylim(-1, len(order) - 0.4)
    a2.set_xlabel("Deaths available for the model")
    a2.set_title("Deaths per sex (what a survival model learns from)",
                 loc="left", color=INK)
   
    a2.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))

    fig.suptitle("Sex distribution in each cancer", x=0.06, ha="left", fontweight="bold",
                 y=1 - 0.2 / height)
    fig.savefig(out / "fig1_representation.png")
    plt.close(fig)


# figure 2
def fig_sex_effect(sex, drop, min_events, min_sex_events, out):
    sex = sex.sort_values("HR").reset_index(drop=True)
    d = drop.pivot(index="cancer", columns="sex", values="n_deaths_after")
    fig, ax = plt.subplots(figsize=(8, 0.55 * len(sex) + 1.5))
    for i, r in sex.iterrows():
        f_d, m_d = int(d.loc[r.cancer, "female"]), int(d.loc[r.cancer, "male"])
        if r["n_events"] < min_events:
            reason = "too few deaths"
        elif min(f_d, m_d) < min_sex_events:
            reason = f"only {min(f_d, m_d)} {'female' if f_d < m_d else 'male'} deaths"
        else:
            reason = ""
        low = bool(reason)
        c = GREY if low else INK
        if estimable(r):
            ax.plot([r.CI_low, r.CI_high], [i, i], color=c, lw=2, solid_capstyle="round")
            ax.plot(r.HR, i, "s", ms=8, mfc="white" if low else c, mec=c, mew=2)
        else:
            ax.text(1, i, "not estimable", ha="center", va="center", color=GREY)
        note = f"{f_d} / {m_d}" + (f"  ({reason})" if low else "")
        ax.annotate(note, (1.02, i), xycoords=("axes fraction", "data"),
                    va="center", fontsize=9, color=GREY if low else INK2)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    log_ticks(ax, [0.25, 0.5, 1, 2, 4])
    ax.set_xlim(0.15, 6)
    ax.set_yticks(range(len(sex)),
                  [c + ("" if s else "*") for c, s in
                   zip(sex["cancer"], sex.get("stage_adjusted", [True] * len(sex)))])
    ax.set_ylim(-0.7, len(sex) - 0.3)
    ax.set_xlabel("Hazard ratio for death, women vs men (adjusted for age and stage)")
    ax.text(0.17, len(sex) - 0.35, "← lower risk for women", fontsize=8, color=INK2)
    ax.text(5.8, len(sex) - 0.35, "higher risk for women →", fontsize=8,
            color=INK2, ha="right")
    ax.annotate("deaths F / M", (1.02, len(sex) - 0.35),
                xycoords=("axes fraction", "data"), fontsize=8, color=INK2)
    ax.set_title("The effect of gender on survival prediction", loc="left", fontweight="bold")
    if "stage_adjusted" in sex and not sex["stage_adjusted"].all():
        # fixed distance (in points) below the axis label, whatever the figure height
        ax.annotate("* no AJCC staging: adjusted for age only", xy=(0, 0),
                    xycoords="axes fraction", xytext=(0, -42),
                    textcoords="offset points", fontsize=8, color=INK2, va="top")
    fig.savefig(out / "fig2_sex_effect.png")
    plt.close(fig)


# figure 3
def fig_precision_gap(cox, min_events, out, term="age_per10y"):
    t = cox[(cox.term == term) & cox.model.isin(["female", "male"])]
    sizes = cox.groupby(["cancer", "model"])["n_events"].first()
    order = sorted(cox.cancer.unique())
    fig, ax = plt.subplots(figsize=(8, 0.75 * len(order) + 1.5))
    for i, c in enumerate(order):
        for sex, off in [("female", 0.17), ("male", -0.17)]:
            r = t[(t.cancer == c) & (t.model == sex)]
            if r.empty:             
                n_ev = sizes.get((c, sex), 0)
                ax.text(1, i + off, f"{sex}: not estimable ({int(n_ev)} deaths)",
                        va="center", fontsize=8, color=GREY)
                continue
            r = r.iloc[0]
            low = r["n_events"] < min_events
            col = COL[sex]
            if estimable(r):
                ax.plot([r.CI_low, r.CI_high], [i + off] * 2, color=col, lw=2,
                        alpha=0.45 if low else 1)
                ax.plot(r.HR, i + off, "o", ms=7, mfc="white" if low else col,
                        mec=col, mew=2)
                ax.text(r.CI_high * 1.06, i + off,
                        f"{int(r.n_events)} deaths" + (" ⚠" if low else ""),
                        va="center", fontsize=8, color=INK2)
            else:
                ax.text(1, i + off, f"{sex}: not estimable ({int(r.n_events)} deaths)",
                        va="center", fontsize=8, color=GREY)
    ax.axvline(1, color=INK2, lw=1, ls="--")
    log_ticks(ax, [0.25, 0.5, 1, 2, 4, 8, 16])
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("Hazard ratio per 10 years of age (95% CI), log scale")
    ax.plot([], [], "o-", color=COL["female"], label="Women-only model")
    ax.plot([], [], "o-", color=COL["male"], label="Men-only model")
    ax.plot([], [], "o", mfc="white", mec=GREY, mew=2, ls="",
            label=f"Hollow = fewer than {min_events} deaths")
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.set_title("The effect of age in each sex",
                 loc="left", fontweight="bold")
    fig.savefig(out / "fig3_precision_gap.png")
    plt.close(fig)


# figure 4
def fig_km(clean_files, out):
    data = pd.concat([pd.read_csv(f) for f in clean_files]).sort_values("cancer")
    cancers = data.cancer.unique().tolist()
    ncol = min(3 if len(cancers) <= 6 else 5, len(cancers))
    nrow = int(np.ceil(len(cancers) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.3 * nrow),
                             squeeze=False, sharey=True)
    for ax, c in zip(axes.flat, cancers):
        d = data[data.cancer == c]
        for sex in ["female", "male"]:
            s = d[d.sex == sex]
            if len(s) == 0:
                continue
            km = KaplanMeierFitter().fit(s.time_years, s.event,
                                         label=f"{sex.capitalize()} (n={len(s)}, "
                                               f"{int(s.event.sum())} deaths)")
            km.plot_survival_function(ax=ax, color=COL[sex], lw=2, ci_alpha=0.15)
        ax.set_title(c, loc="left", fontweight="bold")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1.02)
        ax.set_xlabel("Years since diagnosis")
        ax.set_ylabel("Proportion alive")
        ax.legend(frameon=False, fontsize=7 if len(cancers) > 6 else 8,
                  loc="lower left")
    for ax in list(axes.flat)[len(cancers):]:
        ax.set_visible(False)
    fig.suptitle("Survival by sex (shaded = 95% CI)",
                 x=0.02, ha="left", fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 1 - 0.45 / fig.get_figheight()])
    fig.savefig(out / "fig4_km_by_sex.png")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--representation", required=True)
    p.add_argument("--dropped", required=True)
    p.add_argument("--cox", required=True)
    p.add_argument("--sex-effect", required=True)
    p.add_argument("--clean", nargs="+", required=True)
    p.add_argument("--min-events", type=int, default=40)
    p.add_argument("--min-sex-events", type=int, default=10,
                   help="flag the sex effect if either sex has fewer deaths than this")
    p.add_argument("--outdir", default=".")
    a = p.parse_args()
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)

    rep = pd.read_csv(a.representation)
    drop = pd.read_csv(a.dropped)
    cox = pd.read_csv(a.cox)
    sex = pd.read_csv(a.sex_effect)

    fig_representation(rep, drop, a.min_events, out)
    fig_sex_effect(sex, drop, a.min_events, a.min_sex_events, out)
    fig_precision_gap(cox, a.min_events, out)
    fig_km(a.clean, out)
    print("Wrote 4 figures to", out)


if __name__ == "__main__":
    main()
