import pandas as pd

INPUTS = [
    ("25RESVal1_260819_RES_Validation_checked.tsv", "RESVal1"),
    ("26RESVal2_260819_RES_Validation_checked.tsv", "RESVal2"),
    ("26RESVal3_260819_RES_Validation_checked.tsv", "RESVal3"),
    ("26RESVal4_260819_RES_Validation_checked.tsv", "RESVal4"),
    ("26RESVal5_260819_RES_Validation_checked.tsv", "RESVal5"),
    ("seraseq_RESVal1_checked.tsv", "Seraseq-RESVal1"),
    ("seraseq_RESVal3_checked.tsv", "Seraseq-RESVal3"),
    ("seraseq_RESVal4_checked.tsv", "Seraseq-RESVal4"),
    ("seraseq_RESVal5_checked.tsv", "Seraseq-RESVal5"),
]

DEPTHS = [5, 10, 15, 20, 25]
DEPTH_COLS = ["5M", "10M", "15M", "20M", "25M"]

# "Detected" = the caller reported the fusion for that sample at that depth
# (a non-"." row in the checked TSV). No read-count threshold.

# Rows excluded from the detection-rate denominator: not gene fusions.
# Kept in the rate (genuine fusions, some missed): ALK.
EXCLUDE = {
    "/", "\\", "NEGATIVE CONTROL",   # no fusion
    "SDC4::ROS1",                    # nothing similar in workbook
    "NPV",                           # NPV on WGS
    "TERT",                          # TERTp, no fusion
    "MAP3KB",                        # suspected only
    "PRDM10",                        # ?PRDM10, uncertain
    "NR4A3",                         # overexpression, not a fusion
    "GRM1",                          # overexpression, not a fusion
    "BCOR",                          # BCOR ITD
    "TFE3",                          # single-gene target
    "TRAC",                          # single-gene target (LMO2::TRAC; absent from downsampled workbook)
    "KMT2A", "KMT2A::PTD",           # partial tandem duplication
    "MET", "MET ex 14 skipping",     # exon 14 skipping
    "EGFR", "EGFR Variant III",      # EGFR ITD / vIII deletion
    "Seraseq control",               # placeholder row
    "BCR::ABL1", "TPM::NUP210L",     # sample 25020K0005: 12M unique reads, 11.82M dup
}

def load():
    frames = []
    for path, run in INPUTS:
        try:
            d = pd.read_csv(path, sep="\t", na_values=["."])
        except FileNotFoundError:
            print(f"skip {path}: not found")
            continue
        parts = d["sample_id"].str.extract(r"ds(\d+)M-(.+)")
        d["depth"] = parts[0].astype(int)
        d["sample"] = parts[1]
        d["run"] = run
        d["arriba_det"] = d[
            ["Arriba_split_reads1", "Arriba_split_reads2", "Arriba_discordant_mates"]
        ].notna().any(axis=1)
        d["star_det"] = d[
            ["StarFusion_JunctionReadCount", "StarFusion_SpanningFragCount"]
        ].notna().any(axis=1)
        d["either_det"] = d["arriba_det"] | d["star_det"]
        frames.append(d)
    return pd.concat(frames, ignore_index=True)

def detection_matrix(alld):
    alld = alld[~alld["fusion_name"].isin(EXCLUDE)]
    mat = alld.pivot_table(
        index=["run", "fusion_name", "sample"],
        columns="depth",
        values="either_det",
        aggfunc="first",
    )
    mat = mat.reindex(columns=DEPTHS).replace({True: "Y", False: "N"}).fillna("-")
    mat.columns = DEPTH_COLS
    mat = mat.reset_index().rename(columns={
        "run": "Validation run",
        "fusion_name": "Fusion",
        "sample": "Sample",
    })
    return mat.sort_values(["Validation run", "Fusion", "Sample"]).reset_index(drop=True)

def detection_rate(alld):
    is_real = alld["fusion_name"].notna() & ~alld["fusion_name"].isin(EXCLUDE)
    real = alld[is_real]

    rows = []
    for depth in DEPTHS:
        sub = real[real["depth"] == depth]
        rows.append({
            "Depth": f"{depth}M",
            "Arriba": f"{sub['arriba_det'].mean() * 100:.0f}%",
            "STAR-Fusion": f"{sub['star_det'].mean() * 100:.0f}%",
            "Combined": f"{sub['either_det'].mean() * 100:.0f}%",
            "n_fusions": len(sub),
        })
    return pd.DataFrame(rows)

def main():
    alld = load()

    mat = detection_matrix(alld)
    mat.to_csv("detection_matrix.tsv", sep="\t", index=False)
    print(mat.to_string(index=False))
    print()

    # clinical validation samples and the Seraseq control answer different
    # questions, so rate them separately
    is_seraseq = alld["run"].str.startswith("Seraseq")
    for label, subset, out in [
        ("clinical validation samples", alld[~is_seraseq], "detection_rate_clinical.tsv"),
        ("Seraseq control", alld[is_seraseq], "detection_rate_seraseq.tsv"),
    ]:
        if subset.empty:
            continue
        rate = detection_rate(subset)
        rate.to_csv(out, sep="\t", index=False)
        print(f"--- {label} ---")
        print(rate.to_string(index=False))
        print()


if __name__ == "__main__":
    main()
