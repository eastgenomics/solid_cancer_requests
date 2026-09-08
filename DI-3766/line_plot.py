"""
Faceted line charts: supporting reads (junction + spanning) vs downsample depth, one panel
per sample+fusion, one line per caller (Arriba, STAR-Fusion).

Produces one PNG per checked TSV listed in INPUTS.
"""

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from summary import EXCLUDE   # non-fusion / QC-fail names to skip

# also kept out of the plots (but counted in summary.py's rate): borderline
# fusions whose per-depth line is just noise
PLOT_EXCLUDE = EXCLUDE | {"EWSR1::FLI1", "FGFR1::PLAG1"}

sns.set_theme(style="ticks")

# (input checked TSV, output PNG, plot title)
INPUTS = [
    ("25RESVal1_260819_RES_Validation_checked.tsv", "line_plot_RESVal1.png", "RESVal1"),
    ("26RESVal2_260819_RES_Validation_checked.tsv", "line_plot_RESVal2.png", "RESVal2"),
    ("26RESVal3_260819_RES_Validation_checked.tsv", "line_plot_RESVal3.png", "RESVal3"),
    ("26RESVal4_260819_RES_Validation_checked.tsv", "line_plot_RESVal4.png", "RESVal4"),
    ("26RESVal5_260819_RES_Validation_checked.tsv", "line_plot_RESVal5.png", "RESVal5"),
    ("seraseq_RESVal1_checked.tsv", "line_plot_seraseq_RESVal1.png", "Seraseq - RESVal1"),
    ("seraseq_RESVal3_checked.tsv", "line_plot_seraseq_RESVal3.png", "Seraseq - RESVal3"),
    ("seraseq_RESVal4_checked.tsv", "line_plot_seraseq_RESVal4.png", "Seraseq - RESVal4"),
    ("seraseq_RESVal5_checked.tsv", "line_plot_seraseq_RESVal5.png", "Seraseq - RESVal5"),
]

ORDER = ["STAR-Fusion", "Arriba"]


def make_plot(input_tsv: str, output_png: str, title: str) -> None:
    df = pd.read_csv(input_tsv, sep="\t", na_values=["."])

    # "ds10M-25062S0042" -> depth 10, sample "25062S0042"
    parts = df["sample_id"].str.extract(r"ds(\d+)M-(.+)")
    df["depth"] = parts[0].astype(int)
    df["sample"] = parts[1]

    # Total supporting reads per caller = junction/split reads + spanning pairs.
    #   Arriba      = split_reads1 + split_reads2 + discordant_mates
    #   STAR-Fusion = JunctionReadCount + SpanningFragCount
    # min_count=1 -> NaN only when EVERY part is missing (a real "no call"),
    # not when just one part happens to be blank.
    df["arriba_support"] = df[
        ["Arriba_split_reads1", "Arriba_split_reads2", "Arriba_discordant_mates"]
    ].sum(axis=1, min_count=1)
    df["star_support"] = df[
        ["StarFusion_JunctionReadCount", "StarFusion_SpanningFragCount"]
    ].sum(axis=1, min_count=1)

    # drop non-fusion / QC-fail rows (same list as the summary.py rate)
    df = df[~df["fusion_name"].isin(PLOT_EXCLUDE)]

    # one panel per sample+fusion (same fusion name can occur in >1 sample)
    df["target"] = df["fusion_name"] + "\n" + df["sample"]

    # keep every depth of any target called by either tool at least once, so a
    # depth where neither called it still appears (and becomes a 0 below).
    # Targets never called at any depth
    any_call = df["star_support"].notna() | df["arriba_support"].notna()
    detected_targets = df.loc[any_call, "target"].unique()
    called = df[df["target"].isin(detected_targets)].copy()
    if called.empty:
        print(f"skip {input_tsv}: no called fusions")
        return

    long = called.melt(
        id_vars=["target", "sample", "fusion_name", "depth"],
        value_vars=["arriba_support", "star_support"],
        var_name="caller",
        value_name="support",
    )
    long["caller"] = long["caller"].map({
        "arriba_support": "Arriba",
        "star_support": "STAR-Fusion",
    })

    # a depth where a caller made no call -> plot it as 0 so the line visibly
    # drops to the axis (rather than the point silently vanishing / the line
    # bridging the gap). Every (target, depth, caller) row exists after the
    # melt, so this fills exactly the non-detections.
    long["support"] = long["support"].fillna(0)

    n_panels = long["target"].nunique()
    n_cols = min(4, n_panels)

    # STAR-Fusion drawn first (solid), Arriba last so its dashed line sits on top
    g = sns.relplot(
        data=long,
        x="depth", y="support",
        hue="caller", hue_order=ORDER,
        style="caller", style_order=ORDER,
        dashes={"STAR-Fusion": "", "Arriba": (3, 2)},
        markers={"STAR-Fusion": "X", "Arriba": "o"},
        col="target",
        col_wrap=n_cols,
        kind="line",
        palette={"STAR-Fusion": "#eb6834", "Arriba": "#2a78d6"},
        height=2.4, aspect=1.2,
        facet_kws=dict(sharey=False, sharex=False),
    )

    # open circles for Arriba so a coincident STAR-Fusion 'X' stays visible
    for ax in g.axes.flat:
        for line in ax.lines:
            if line.get_color() == "#2a78d6":
                line.set_markerfacecolor("none")
                line.set_markeredgewidth(1.4)
                line.set_linewidth(1.6)

    g.set_titles("{col_name}")
    g.set(xticks=[5, 10, 15, 20, 25])
    for ax in g.axes.flat:
        ax.set_ylim(bottom=0)
        ax.set_xlabel("depth (M reads)", fontsize=9)
        ax.set_ylabel("supporting reads", fontsize=9)
        ax.tick_params(labelbottom=True, labelleft=True)

    g.legend.set_title(None)
    sns.move_legend(g, "lower center", bbox_to_anchor=(0.5, -0.04), ncol=2, frameon=False)
    g.figure.suptitle(
        f"{title} - fusion supporting reads (junction + spanning) vs read depth",
        y=1.02, fontsize=13,
    )
    g.figure.tight_layout()

    g.savefig(output_png, dpi=150, bbox_inches="tight")
    plt.close(g.figure)
    print(f"wrote {output_png}  ({n_panels} panels)")


def main() -> None:
    for input_tsv, output_png, title in INPUTS:
        try:
            make_plot(input_tsv, output_png, title)
        except FileNotFoundError:
            print(f"skip {input_tsv}: not found")


if __name__ == "__main__":
    main()
