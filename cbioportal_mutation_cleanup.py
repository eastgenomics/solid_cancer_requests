"""
Bit of code to do some cleanup of the file downloaded from cbioportal with
mutation data, expects a folder of tsvs downloaded from the genes tables.

Splits out the single annotation and functional impact columns to separate
columns of each field they contain.
"""

from pathlib import Path

import pandas as pd

download_path = "/home/jethro/Projects/solid_cancer_requests/files"
output_path = Path(download_path).parent / "cleaned_files"

output_path.mkdir(parents=True, exist_ok=True)

files = list(Path("").rglob("*.tsv"))
print(f"Found {len(files)} in {download_path}")

for idx, tsv in enumerate(files, 1):
    print(f"[{idx}/{len(files)}] Cleaning {tsv}")
    df = pd.read_csv(tsv, sep="\t")

    # split out annotation to separate fields
    annotation_col_names = [
        x.split(":")[0] for x in df["Annotation"].iloc[0].split(";")
    ]
    df[annotation_col_names] = df["Annotation"].str.split(";", expand=True)

    # remove the field name that is now the column name
    df[annotation_col_names] = df[annotation_col_names].applymap(
        lambda x: x.split(":", 1)[1]
    )

    df = df.join(
        df["OncoKB"].str.split(",", expand=True).add_prefix("OncoKB_")
    )
    df = df.join(df["CIViC"].str.split(",", expand=True).add_prefix("CIViC_"))
    df.drop(["Annotation", "OncoKB", "CIViC"], axis=1, inplace=True)

    # split out annotation to separate fields
    impact_col_names = [
        x.split(":")[0] for x in df["Functional Impact"].iloc[0].split(";")
    ]
    df[impact_col_names] = df["Functional Impact"].str.split(";", expand=True)

    # remove the field name that is now the column name
    df[impact_col_names] = df[impact_col_names].applymap(
        lambda x: x.split(":", 1)[1]
    )

    impact_cols = ["MutationAssessor", "SIFT", "Polyphen-2", "AlphaMissense"]

    for col in impact_cols:
        df = df.join(df[col].str.split(",", expand=True).add_prefix(f"{col}_"))

    df.drop(impact_cols + ["Functional Impact"], axis=1, inplace=True)

    # cols to just rename
    rename_cols = {
        "OncoKB_0": "OncoKB",
        "OncoKB_1": "OncoKB_level",
        "OncoKB_2": "OncoKB_resistance",
    }

    # cols to rename and split on ':'
    split_rename_cols = {
        "CIViC_0": "CIViC_diagnosticCount",
        "CIViC_1": "CIViC_predictiveCount",
        "CIViC_2": "CIViC_prognosticCount",
        "CIViC_3": "CIViC_predisposingCount",
        "CIViC_4": "CIViC_oncogenicCount",
        "CIViC_5": "CIViC_functionalCount",
        "MutationAssessor_0": "MutationAssessor_impact",
        "MutationAssessor_1": "MutationAssessor_score",
        "SIFT_0": "SIFT_impact",
        "SIFT_1": "SIFT_score",
        "Polyphen-2_0": "Polyphen-2_impact",
        "Polyphen-2_1": "Polyphen-2_score",
        "AlphaMissense_0": "AlphaMissense_pathogenicity",
        "AlphaMissense_1": "AlphaMissense_score",
    }

    for old, new in split_rename_cols.items():
        df[old] = (
            df[old]
            .fillna("")
            .apply(lambda x: x.split(":", 1)[1] if ":" in x and x else x)
        )

    df.rename(columns={**rename_cols, **split_rename_cols}, inplace=True)

    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    df.to_csv(
        f"{output_path / tsv.replace('.tsv', '_cleaned.tsv')}",
        index=False,
        sep="\t",
    )
