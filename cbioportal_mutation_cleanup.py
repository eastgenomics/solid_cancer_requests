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

files = sorted(Path(download_path).rglob("*.tsv"))
print(f"Found {len(files)} in {download_path}")

for idx, tsv in enumerate(files, 1):
    print(f"[{idx}/{len(files)}] Cleaning {tsv}")

    df = pd.read_csv(tsv, sep="\t")

    if df.empty:
        print(f"Warning: file empty: {tsv}")
        continue

    if "Annotation" in df.columns:
        # the annotation column is semi-colon separated but reVUE column
        # can have semi_colons in it (:sadpepe:). Therefore use regex to split
        # out the field name plus contents until the next field (regardless of
        # the order)
        annotation_df = (
            df["Annotation"]
            .str.extractall(
                r"((OncoKB|reVUE|CIViC|CancerHotspot|3DHotspot):.*?)(?=(OncoKB|reVUE|CIViC|CancerHotspot|3DHotspot):|$)"
            )[0]
            .reset_index()
            .pivot(index="level_0", columns="match", values=0)
            .reset_index()
        )
        annotation_df = annotation_df.applymap(
            lambda x: x.rstrip(";") if x and isinstance(x, str) else x
        )
        annotation_df.drop(columns=["level_0"], inplace=True)

        for col in annotation_df.columns:
            name = list(
                set([x.split(":", 1)[0] for x in annotation_df[col].to_list()])
            )

            if len(name) != 1:
                print("Something borked in annotation data: {name}")
                exit()

            annotation_df[col] = annotation_df[col].apply(
                lambda x: x.split(":", 1)[1]
            )
            annotation_df.rename(columns={col: name[0]}, inplace=True)

        df = df.join(annotation_df)

        if "OncoKB" in df.columns:
            # split out the multiple OncoKB fields to separate columns
            df = df.join(
                df["OncoKB"].str.split(",", expand=True).add_prefix("OncoKB_")
            )
            df.drop(["OncoKB"], axis=1, inplace=True)

            # do some specific clean up
            df["OncoKB_level"] = df["OncoKB_level"].str.replace("level ", "")
            df["OncoKB_resistance"] = df["OncoKB_resistance"].str.replace(
                "resistance ", ""
            )

        if "CIViC" in df.columns:
            # split out the multiple CIViC fields to separate columns
            df = df.join(
                df["CIViC"].str.split(",", expand=True).add_prefix("CIViC_")
            )
            df.drop(["CIViC"], axis=1, inplace=True)

        df.drop(["Annotation"], axis=1, inplace=True)

    if "Functional Impact" in df.columns:
        # split out functional impact to separate fields
        impact_col_names = [
            x.split(":")[0] for x in df["Functional Impact"].iloc[0].split(";")
        ]
        df[impact_col_names] = df["Functional Impact"].str.split(
            ";", expand=True
        )

        # remove the field name that is now the column name
        df[impact_col_names] = df[impact_col_names].applymap(
            lambda x: x.split(":", 1)[1]
        )

        impact_cols = [
            "MutationAssessor",
            "SIFT",
            "Polyphen-2",
            "AlphaMissense",
        ]

        for col in impact_cols:
            df = df.join(
                df[col].str.split(",", expand=True).add_prefix(f"{col}_")
            )

        df.drop(impact_cols + ["Functional Impact"], axis=1, inplace=True)

    # cols to just rename
    rename_cols = {
        "OncoKB_0": "OncoKB",
        "OncoKB_1": "OncoKB_level",
        "OncoKB_2": "OncoKB_resistance",
    }

    # cols to rename and split on ':' to remove the title from the value
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
        if old in df.columns:
            df[old] = (
                df[old]
                .fillna("")
                .apply(lambda x: x.split(":", 1)[1] if ":" in x and x else x)
            )

    df.rename(columns={**rename_cols, **split_rename_cols}, inplace=True)

    df = df.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    df.to_csv(
        str(output_path / Path(tsv).stem) + ".cleaned.tsv",
        index=False,
        sep="\t",
    )
