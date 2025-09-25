"""
Script to parse tables from the cancer WGS supplementary HTML files.

Aggregates the top 4 case information tables to a single row per case
(~30 columns) and the 5th table of germline and tumour info which are then
written to separate TSV files. Initial ask of this was to quantify some
metrics across a large number of WGS cases we have analysed (i.e turn
around times, tumour QC metrics etc.)
"""

import argparse
from datetime import datetime
from functools import reduce
from pathlib import Path
from typing import List, Tuple

import pandas as pd


def extract_tables(file: Path) -> List[pd.DataFrame]:
    """
    Parameters
    ----------
    file : pathlib.Path
        Path to file to read from

    Returns
    -------
    List[pd.DataFrame]
        List of dataframes read from file
    """
    return pd.read_html(file)


def select_tables(
    tables: List[pd.DataFrame],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Splits list of tables by the first 4 which contain patient / case info
    and the 5th table with germline / tumour sequencing quality info.

    The first 4 tables have one row each and columns are concatenated to
    a single df, the other will contain one row each for germline and
    tumour info.

    Parameters
    ----------
    tables : list
        List of dataframes read from HTML

    Returns
    -------
    pd.DataFrame
        Single row DataFrame of columns from first 4 tables
    pd.DataFrame
        Two row DataFrame from fifth table in HTML
    """
    if len(tables) < 5:
        raise ValueError(f"Expected at least 5 tables, got {len(tables)}")

    # join data from first 4 tables, ensure no duplicate columns names
    # by adding .1 suffix (i.e where sample ID is in 2 tables)
    sample_tables = reduce(
        lambda left, right: pd.merge(
            left,
            right,
            how="inner",
            suffixes=[None, ".1"],
            left_index=True,
            right_index=True,
        ),
        tables[:4],
    )

    qual_table = tables[4]

    # add case / sample identifiers to keep unique rows per file later
    columns_to_add = {
        0: "Referral ID",
        1: "Patient ID",
        2: "Histopathology or SIHMDS LAB ID",
        3: "Sample ID",
        4: "Sample ID.1",
    }

    for position, name in columns_to_add.items():
        qual_table.insert(position, name, [sample_tables.iloc[0][name]] * 2)

    return sample_tables, qual_table


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--paths",
        help="path(s) to search for supplementary HTMLs within",
        required=True,
        nargs="+",
    )
    args = parser.parse_args()

    all_file_paths = {}
    duplicate_files = {}

    for path in args.paths:
        if not Path(path).is_dir():
            raise ValueError(
                f"Path does not exist or is not a directory: {path}"
            )

        files = list(Path(path).rglob("*.supplementary.html"))
        print(f"Found {len(files)} in {path}")

        for html in files:
            if html.name in all_file_paths:
                print(
                    f"Warning: {html.name} already found at"
                    f" {all_file_paths[html.name]}, skipping..."
                )
                duplicate_files[html.name] = html.resolve().parent
            else:
                all_file_paths[html.name] = html.resolve().parent

    sample_detail_dfs = []
    germline_tumour_quality_dfs = []
    error_files = []

    for idx, file in enumerate(all_file_paths.items(), 1):
        print(f"[{idx}/{len(all_file_paths)}] {file[0]}")

        try:
            tables = extract_tables(file[1] / Path(file[0]))
            sample_tables_df, qual_table_df = select_tables(tables)

            sample_detail_dfs.append(sample_tables_df)
            germline_tumour_quality_dfs.append(qual_table_df)
        except Exception as exc:
            print(f"Error processing {file}: {exc}")
            error_files.append(file)

    all_sample_details_df = pd.concat(
        sample_detail_dfs, ignore_index=True, sort=False
    )
    all_germline_tumour_quality_dfs = pd.concat(
        germline_tumour_quality_dfs, ignore_index=True, sort=False
    )

    print(all_sample_details_df)
    print(all_germline_tumour_quality_dfs)

    print(f"Total sample detail rows: {len(all_sample_details_df)}")
    print(
        "Total germline / tumour quality rows:"
        f" {len(all_germline_tumour_quality_dfs)}"
    )

    all_sample_details_df = all_sample_details_df.drop_duplicates()
    all_germline_tumour_quality_dfs = (
        all_germline_tumour_quality_dfs.drop_duplicates()
    )

    print(
        "Total sample detail rows after dropping duplicates:"
        f" {len(all_sample_details_df)}"
    )
    print(
        "Total germline / tumour quality rows after dropping duplicates:"
        f" {len(all_germline_tumour_quality_dfs)}"
    )

    now = datetime.now().strftime("%Y%m%d_%H%M")
    all_sample_details_df.to_csv(
        f"cancer_wgs_html_extract_sample_details_{now}.tsv",
        sep="\t",
        index=False,
    )
    all_germline_tumour_quality_dfs.to_csv(
        f"cancer_wgs_html_extract_germline_tumour_quality_{now}.tsv",
        sep="\t",
        index=False,
    )

    print(
        f"Found {len(all_file_paths)} unique HTML files\nIgnored"
        f" {len(duplicate_files)} duplicate files in different locations"
    )

    with open("cancer_wgs_html_extract_extracted_files.txt", "w") as fh:
        for file, path in all_file_paths.items():
            fh.write(f"{path}/{file}\n")

    if duplicate_files:
        with open("cancer_wgs_html_extract_duplicate_files.txt", "w") as fh:
            for file, path in duplicate_files.items():
                fh.write(f"{path}/{file}\n")

    if error_files:
        errors = "\n\t".join([f"{file[1]}/{file[0]}" for file in error_files])
        print(f"{len(error_files)} files failed to be read from:{errors}")


if __name__ == "__main__":
    main()
