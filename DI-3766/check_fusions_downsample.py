"""
Cross-reference a CSV of expected/known RNA fusions against the Arriba and
STAR-Fusion sheets of a DNAnexus RNA-fusions workbook, and write a combined
TSV showing what each tool reported for every expected fusion.

Requires an active DNAnexus auth context (e.g. `dx login`) before running.
"""

import argparse
import io
import re
import subprocess
import sys

import openpyxl
import pandas as pd

ARRIBA_SPECIMEN_COLUMN = "file_name"
ARRIBA_GENE1_COLUMN = "#gene1"
ARRIBA_GENE2_COLUMN = "gene2"
ARRIBA_OUTPUT_COLUMNS = [
    "breakpoint1",
    "breakpoint2",
    "split_reads1",
    "split_reads2",
    "discordant_mates",
    "coverage1",
    "coverage2",
    "confidence",
]

STARFUSION_FILENAME_COLUMN = "file_name"
STARFUSION_FUSION_NAME_COLUMN = "#FusionName"
STARFUSION_OUTPUT_COLUMNS = [
    "LeftBreakpoint",
    "RightBreakpoint",
    "JunctionReadCount",
    "SpanningFragCount",
]

MISSING_VALUE = "."

# read-support columns used to pick the strongest row when a fusion is reported
# at more than one breakpoint for a sample
ARRIBA_SUPPORT_COLUMNS = ["split_reads1", "split_reads2", "discordant_mates"]
STARFUSION_SUPPORT_COLUMNS = ["JunctionReadCount", "SpanningFragCount"]


def gene_token(cell) -> str:
    """One gene symbol from a cell: stripped and upper-cased, '' if blank."""
    text = str(cell).strip().upper()
    return "" if text in ("", ".", "NAN") else text


def fusion_pair(fusion_name: str) -> tuple:
    """
    Ordered gene tuple for an expected fusion name:
    'GENE1::GENE2' or 'GENE1--GENE2' -> ('GENE1', 'GENE2'); single-gene
    'GENE' -> ('GENE',).
    """
    parts = [gene_token(p) for p in re.split(r"::|--", str(fusion_name))]
    return tuple(p for p in parts if p)


def find_matches(df: pd.DataFrame, sample_id: str, expected: tuple) -> pd.DataFrame:
    """
    Rows for this sample matching `expected`.
    Two-gene expected: the row's genes must be in the SAME order (5' then 3');
    'A--B' does not match an expected 'B::A'.
    Single-gene expected: any row where that gene is either partner.
    """
    in_sample = df["file_name"].str.contains(sample_id, regex=False, na=False)
    if len(expected) == 1:
        gene = expected[0]
        return df[in_sample & df["gene_set"].map(lambda s: gene in s)]
    return df[in_sample & df["ordered_pair"].map(lambda p: p == expected)]


def row_with_most_support(matches: pd.DataFrame, support_columns: list):
    """
    From one or more rows matching the same fusion, return the row whose
    supporting-read columns sum highest. Breakpoint scatter for one fusion
    produces several rows; this keeps the strongest breakpoint instead of
    whichever row happens to be first in the sheet (which varies by depth).
    """
    support = (
        matches[support_columns].apply(pd.to_numeric, errors="coerce").sum(axis=1)
    )
    return matches.loc[support.idxmax()]


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
    parser.add_argument(
        "--run_name",
        type=str,
        required=False,
        default=None,
        help="Explicit run label to select rows of --expected_fusions_tsv, "
        "matched against its 'Run' column (e.g. '26-RESVal2'). Use this when "
        "the workbook filename does not uniquely identify the run (several "
        "downsampling workbooks share the name "
        "'260819_RES_Validation_fusion_workbook.xlsx'). If omitted, the run "
        "name is derived from the workbook filename and matched against "
        "'Run_name'.",
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
                f"Sheet '{sheet_name}' not found in workbook. "
                f"Sheets present: {sheets}"
            )

    return workbook[arriba_sheet], workbook[starfusion_sheet]


def build_arriba_fusion_names(arriba_worksheet) -> pd.DataFrame:
    """
    Turn the Arriba worksheet into a dataframe and add the columns used for
    matching:
      fusion_name  'gene1::gene2' string
      gene_set     {gene1, gene2}

    Parameters
    ----------
    arriba_worksheet : openpyxl worksheet
        The Arriba sheet.

    Returns
    -------
    pd.DataFrame
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
    g1 = [gene_token(g) for g in gene1]
    g2 = [gene_token(g) for g in gene2]
    arriba_df["ordered_pair"] = list(zip(g1, g2))      # (5' gene, 3' gene)
    arriba_df["gene_set"] = [{a, b} - {""} for a, b in zip(g1, g2)]
    return arriba_df

def build_starfusion_dataframe(starfusion_worksheet) -> pd.DataFrame:
    """
    Turn the STAR-Fusion worksheet into a dataframe and add the same
    fusion_name / gene_set columns as build_arriba_fusion_names.

    Parameters
    ----------
    starfusion_worksheet : openpyxl worksheet
        The STAR-Fusion sheet.

    Returns
    -------
    pd.DataFrame
    """
    starfusion_df = pd.DataFrame(starfusion_worksheet.values)
    starfusion_df.columns = starfusion_df.iloc[0]
    starfusion_df = starfusion_df[1:].reset_index(drop=True)
    starfusion_df["fusion_name"] = (
        starfusion_df["#FusionName"].astype(str).str.replace("--", "::")
    )
    starfusion_df["ordered_pair"] = starfusion_df["#FusionName"].map(fusion_pair)
    starfusion_df["gene_set"] = starfusion_df["ordered_pair"].map(set)

    return starfusion_df


def match_fusion(
    df: pd.DataFrame,
    sample_id: str,
    fusion_name: str,
    output_columns: list,
    support_columns: list,
    tool: str,
) -> dict:
    """
    Look up an expected fusion in one tool's sheet for a given sample.

    df               sheet from build_arriba_fusion_names / build_starfusion_dataframe
    sample_id        substring of the file_name column for this sample
    fusion_name      'GENE1::GENE2', 'GENE1--GENE2', or single 'GENE'. Two-gene
                     names are matched in order; a single gene matches either
                     partner.
    output_columns   sheet columns to copy into the result
    support_columns  columns summed to pick the strongest row on a tie
    tool             label for the multi-match warning

    Returns a dict of output_columns -> value plus 'n_matches' (0 if no row
    matched). Values are MISSING_VALUE when nothing matched.
    """
    missing = {column: MISSING_VALUE for column in output_columns}
    missing["n_matches"] = 0

    expected = fusion_pair(fusion_name)
    if not expected:
        return missing

    matches = find_matches(df, sample_id, expected)
    if matches.empty:
        return missing

    if len(matches) > 1:
        print(
            f"Warning: {len(matches)} {tool} matches for sample '{sample_id}' "
            f"fusion '{fusion_name}'; using the one with most support."
        )

    matched_row = row_with_most_support(matches, support_columns)
    result = {column: matched_row[column] for column in output_columns}
    result["n_matches"] = len(matches)
    return result


def match_arriba(arriba_df: pd.DataFrame, sample_id: str, fusion_name: str) -> dict:
    return match_fusion(
        arriba_df, sample_id, fusion_name,
        ARRIBA_OUTPUT_COLUMNS, ARRIBA_SUPPORT_COLUMNS, "Arriba",
    )


def match_starfusion(
    starfusion_df: pd.DataFrame, sample_id: str, fusion_name: str
) -> dict:
    return match_fusion(
        starfusion_df, sample_id, fusion_name,
        STARFUSION_OUTPUT_COLUMNS, STARFUSION_SUPPORT_COLUMNS, "STAR-Fusion",
    )


def main() -> None:
    args = parse_arguments()

    if args.output and len(args.workbook_file_id) > 1:
        sys.exit(
            "--output is a single path; with several --workbook_file_id values "
            "each run would overwrite it. Drop --output (files are named from "
            "the run) or process one workbook at a time."
        )

    for workbook_file_id in args.workbook_file_id:

        try:
            cmd = (
            f"dx describe {workbook_file_id} --name "
            )    
            run_name_output = subprocess.run(cmd, shell=True,
                                capture_output=True)
            run_name_list = run_name_output.stdout.decode("utf-8").strip().split("_")[0:4]
            run_name = args.run_name or "_".join(run_name_list)

            print(run_name)

            cmd = (
            f"dx cat {workbook_file_id}"
            )    
            output = subprocess.run(cmd, shell=True,
                                capture_output=True)

            sample_excel = openpyxl.load_workbook(io.BytesIO(output.stdout),
                                            data_only=False) 
            csv_stem = run_name
            output_path = args.output or f"{csv_stem}_checked.tsv"
        except pd.errors.EmptyDataError as error:
            print(f"Error reading file {args.workbook_file_id}: {error}")
            print(f"Check archival status of {args.workbook_file_id}")
            continue

        arriba_worksheet, starfusion_worksheet = load_tool_sheets(
                sample_excel, args.arriba_sheet_name, args.starfusion_sheet_name
        )

        arriba_df = build_arriba_fusion_names(arriba_worksheet)
        starfusion_df = build_starfusion_dataframe(starfusion_worksheet)

        expected_fusions = pd.read_csv(args.expected_fusions_tsv, sep="\t")
        if args.run_name:
            # explicit label -> exact match on the 'Run' column
            expected_fusions_in_run = expected_fusions[
                expected_fusions["Run"].astype(str).str.strip() == args.run_name
            ]
        else:
            # derived from the workbook filename -> substring match on 'Run_name'
            expected_fusions_in_run = expected_fusions[
                expected_fusions["Run_name"].str.contains(run_name, na=False)
            ]

        output_rows = []
        for _, expected in expected_fusions_in_run.iterrows():
            sample_id = str(expected["Cases (SP)"]).strip()
            fusion_name = str(expected["Fusion detected"]).strip()
            print(f"Checking sample '{sample_id}' fusion '{fusion_name}'...")

            arriba_result = match_arriba(arriba_df, sample_id, fusion_name)
            starfusion_result = match_starfusion(starfusion_df, sample_id, fusion_name)

            row = {"sample_id": sample_id, "fusion_name": fusion_name}
            row.update({f"Arriba_{k}": v for k, v in arriba_result.items()})
            row.update({f"StarFusion_{k}": v for k, v in starfusion_result.items()})
            output_rows.append(row)

        output_df = pd.DataFrame(output_rows)
        output_df.to_csv(output_path, sep="\t", index=False)
        print(f"Wrote {len(output_df)} rows to {output_path}")


if __name__ == "__main__":
    main()
