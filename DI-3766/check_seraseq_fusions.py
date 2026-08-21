"""
Cross-reference a CSV of expected/known RNA fusions against the Arriba and
STAR-Fusion sheets of a DNAnexus RNA-fusions workbook, and write a combined
TSV showing what each tool reported for every expected fusion.

Requires an active DNAnexus auth context (e.g. `dx login`) before running.
"""

import argparse
import re
import sys
import openpyxl
import io
import subprocess
import dxpy
import pandas as pd

ARRIBA_SPECIMEN_COLUMN = "file_name"
ARRIBA_GENE1_COLUMN = "#gene1"
ARRIBA_GENE2_COLUMN = "gene2"
ARRIBA_OUTPUT_COLUMNS = [
    "present"
]

STARFUSION_FILENAME_COLUMN = "file_name"
STARFUSION_FUSION_NAME_COLUMN = "#FusionName"
STARFUSION_OUTPUT_COLUMNS = [
    "present"
]

MISSING_VALUE = "."


def parse_arguments() -> argparse.Namespace:
    """
    Parse command-line arguments.

    Returns
    -------
    args: argparse.Namespace
        Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=(
            "Check expected RNA fusions against Arriba and STAR-Fusion "
            "calls in a DNAnexus RNA-fusions workbook."
        )
    )
    parser.add_argument(
        "--workbook_file_id",
        type=str,
        nargs='+',
        required=True,
        help="DNAnexus file ID of the RNA-fusions workbook (e.g. "
        "file-J9jKpfj4yXj1KZ3GK8bvj5Vx).",
    )
    parser.add_argument(
        "--project_id",
        type=str,
        required=False,
        default=None,
        help="DNAnexus project ID to qualify the file ID canonically. "
        "Recommended to avoid resolving the wrong project copy of the "
        "file; if omitted the file's default project context is used.",
    )
    parser.add_argument(
        "--expected_fusions_tsv",
        type=str,
        required=True,
        help="Local TSV of expected fusions, with 'sample_id' and "
        "'fusion_name' columns (one row per expected fusion per sample).",
    )
    parser.add_argument(
        "--arriba_sheet_name",
        type=str,
        required=False,
        default="Arriba",
        help="Name of the Arriba sheet in the workbook (default: 'Arriba').",
    )
    parser.add_argument(
        "--starfusion_sheet_name",
        type=str,
        required=False,
        default="STAR-Fusion",
        help="Name of the STAR-Fusion sheet in the workbook (default: "
        "'STAR-Fusion').",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=False,
        default=None,
        help="Path to write the output TSV to (default: derived from "
        "--expected_fusions_tsv, e.g. 'expected_fusions_checked.tsv').",
    )
    return parser.parse_args()





def load_tool_sheets(
    workbook, arriba_sheet: str, starfusion_sheet: str
) -> tuple:
    """
    Load the Arriba and STAR-Fusion sheets from the workbook.

    Parameters
    ----------

    arriba_sheet : str
        Name of the Arriba sheet
    starfusion_sheet : str
        Name of the STAR-Fusion sheet

    Returns
    -------
    arriba_df : pd.DataFrame
    starfusion_df : pd.DataFrame

    Raises
    ------
    RuntimeError
        If either sheet name is not found in the workbook
    """
    sheets = workbook.sheetnames
    for sheet_name in (arriba_sheet, starfusion_sheet):
        if sheet_name not in sheets:
            raise RuntimeError(
                f"Sheet '{sheet_name}' not found in workbook. Available "
                f"sheets: {list(sheets.keys())}"
            )

    return workbook[arriba_sheet], workbook[starfusion_sheet]


def build_arriba_fusion_names(arriba_worksheet) -> pd.DataFrame:
    """
    Add a 'fusion_name' column to the Arriba dataframe as 'gene1--gene2',
    and strip whitespace from the sample and gene columns used for matching.

    Parameters
    ----------
    arriba_worksheet : pd.DataFrame
        Raw Arriba sheet data

    Returns
    -------
    pd.DataFrame
        Arriba data with an added 'fusion_name' column
    """
    arriba_df = pd.DataFrame(arriba_worksheet.values)
    arriba_df.columns = arriba_df.iloc[0]
    arriba_df = arriba_df[1:].reset_index(drop=True)

    arriba_df[ARRIBA_SPECIMEN_COLUMN] = (
        arriba_df[ARRIBA_SPECIMEN_COLUMN].astype(str).str.strip()
    )
    gene1 = arriba_df[ARRIBA_GENE1_COLUMN].astype(str).str.strip()
    gene2 = arriba_df[ARRIBA_GENE2_COLUMN].astype(str).str.strip()
    arriba_df["fusion_name"] = gene1 + "::" + gene2
    return arriba_df

def build_starfusion_dataframe(starfusion_worksheet) -> pd.DataFrame:
    """
    Add a 'fusion_name' column to the STAR-Fusion dataframe as 'gene1::gene2',
    and strip whitespace from the sample and gene columns used for matching.

    Parameters
    ----------
    starfusion_worksheet : pd.DataFrame
        Raw STAR-Fusion sheet data

    Returns
    -------
    pd.DataFrame
        STAR-Fusion data with an added 'fusion_name' column
    """
    starfusion_df = pd.DataFrame(starfusion_worksheet.values)
    starfusion_df.columns = starfusion_df.iloc[0]
    starfusion_df = starfusion_df[1:].reset_index(drop=True)
    starfusion_df["fusion_name"] =  starfusion_df["#FusionName"].astype(str).str.replace("--", "::")

    return starfusion_df


def match_arriba(arriba_df: pd.DataFrame, fusion_name: str) -> dict:
    """
    Look up an expected fusion in the Arriba sheet for a given sample.

    Parameters
    ----------
    arriba_df : pd.DataFrame
        Arriba data with a 'fusion_name' column (see build_arriba_fusion_names)
    sample_id : str
        Sample identifier to match exactly against the SPECIMEN column
    fusion_name : str
        Expected fusion name (order-sensitive, e.g. 'GENE1::GENE2')

    Returns
    -------
    dict
        Mapping of ARRIBA_OUTPUT_COLUMNS to matched values, or "." for
        every column if no match was found
    """
    matches = arriba_df[
        arriba_df["fusion_name"] == fusion_name
    ]

    if matches.empty:
        matches = arriba_df[
        (arriba_df["#gene1"] == fusion_name) | (arriba_df["gene2"] == fusion_name)
        ]
        if matches.empty:
            return {column: MISSING_VALUE for column in ARRIBA_OUTPUT_COLUMNS}

    if len(matches) > 1:
        print(
            f"Warning: {len(matches)} Arriba matches found for sample "
            f" fusion '{fusion_name}'; using the first match."
        )

    matched_row = matches.iloc[0]
    print(f"Matched Arriba row for sample fusion '{fusion_name}':")
    matched_row["present"] = "yes"
    return {column: matched_row[column] for column in ARRIBA_OUTPUT_COLUMNS}


def match_starfusion(
    starfusion_df: pd.DataFrame, fusion_name: str
) -> dict:
    """
    Look up an expected fusion in the STAR-Fusion sheet for a given sample.

    Parameters
    ----------
    starfusion_df : pd.DataFrame
        Raw STAR-Fusion sheet data
    sample_id : str
        Sample identifier expected as a substring of the file_name column
    fusion_name : str
        Expected fusion name (order-sensitive, matched against #FusionName)

    Returns
    -------
    dict
        Mapping of STARFUSION_OUTPUT_COLUMNS to matched values, or "." for
        every column if no match was found
    """
    matches = starfusion_df[starfusion_df["fusion_name"] == fusion_name]

    if matches.empty:
        matches = starfusion_df[(starfusion_df["#FusionName"].str.split("--").str[0] == fusion_name) | (starfusion_df["#FusionName"].str.split("--").str[1] == fusion_name)]

        if matches.empty:
            return {column: MISSING_VALUE for column in STARFUSION_OUTPUT_COLUMNS}

    if len(matches) > 1:
        print(
            f"Warning: {len(matches)} STAR-Fusion matches found for sample "
            f"' fusion '{fusion_name}'; using the first match."
        )

    matched_row = matches.iloc[0]
    matched_row["present"] = "yes"
    return {column: matched_row[column] for column in STARFUSION_OUTPUT_COLUMNS}


def main() -> None:
    args = parse_arguments()
    output_path = args.output
    output_rows = []
    for workbook_file_id in args.workbook_file_id:

        try:
            cmd = (
            f"dx describe {workbook_file_id} --name "
            )    
            run_name_output = subprocess.run(cmd, shell=True,
                                capture_output=True)
            run_name_list = run_name_output.stdout.decode("utf-8").strip().split("_")[0:5]
            run_name = "_".join(run_name_list)

            print(run_name)

            cmd = (
            f"dx cat {workbook_file_id}"
            )    
            output = subprocess.run(cmd, shell=True,
                                capture_output=True)

            sample_excel = openpyxl.load_workbook(io.BytesIO(output.stdout),
                                            data_only=False) 
            csv_stem = run_name
            output_path = f"control_comparison_checked.tsv"
        except pd.errors.EmptyDataError as error:
            print(f"Error reading file {args.workbook_file_id}: {error}")
            print(f"Check archival status of {args.workbook_file_id}")
            continue

        arriba_worksheet, starfusion_worksheet = load_tool_sheets(
                sample_excel, args.arriba_sheet_name, args.starfusion_sheet_name
        )

        arriba_df = build_arriba_fusion_names(arriba_worksheet)
        starfusion_df = build_starfusion_dataframe(starfusion_worksheet)

        expected_fusions_in_run = pd.read_csv(args.expected_fusions_tsv, sep="\t")
        print(expected_fusions_in_run)

        for _, expected in expected_fusions_in_run.iterrows():
            fusion_name = str(expected["Fusion detected"])
            print(f"Checking sample fusion '{fusion_name}'...")

            arriba_result = match_arriba(arriba_df, fusion_name)
            starfusion_result = match_starfusion(starfusion_df, fusion_name)

            row = {"run_name": run_name, "fusion_name": fusion_name}
            row.update({f"Arriba_{k}": v for k, v in arriba_result.items()})
            row.update({f"StarFusion_{k}": v for k, v in starfusion_result.items()})
            output_rows.append(row)
    
    output_df = pd.DataFrame(output_rows)
    print(output_df)
    output_pivot = pd.pivot_table(output_df,index=['fusion_name'],
                                  values=["Arriba_present", "StarFusion_present"],
                                  aggfunc="first",
                                  columns=['run_name'])
    
    output_pivot.to_csv(output_path, sep="\t")
    print(f"Wrote {len(output_df)} rows to {output_path}")


if __name__ == "__main__":
    main()
