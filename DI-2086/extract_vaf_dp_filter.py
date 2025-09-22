import pandas as pd
import subprocess
import argparse


def read_features_file(features_file, sample_name):
    """
    Read a sompy features.csv file and remove duplicate entries (duplicate variants
    due to multiple annotations for more than one transcript)

    Parameters
    ----------
    features_file : str
        Path to the input CSV file containing variant features from sompy.
    sample_name : str
        Sample name corresponding to the features_file

    Returns
    -------
    df : pandas.DataFrame
        Table containing unique variants from the features CSV file

    Examples
    --------
    >>> sample_name, records = parse_features_file('sample123.features.csv')
    """
    # Read the features.csv file
    df = pd.read_csv(features_file)
    
    # Check if required columns exist
    required_cols = ['REF', 'ALT', 'REF.truth', 'ALT.truth']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        print(f"Warning: Missing columns: {missing_cols}")
        return df
    
    # Replace NaN values in REF column with REF.truth
    df.loc[df['REF'].isna(), 'REF'] = df.loc[df['REF'].isna(), 'REF.truth']
    
    # Replace NaN values in ALT column with ALT.truth
    df.loc[df['ALT'].isna(), 'ALT'] = df.loc[df['ALT'].isna(), 'ALT.truth']
    
    # Add the sample name to the records dataframe
    df['sample'] = sample_name
    
    # Return dataframe with dropped duplicates
    return df.drop_duplicates(subset=['CHROM','POS', 'REF', 'ALT'])


def extract_variant_data_df(vcf_path):
    """
    Extract variant depth (DP), allele frequency (VAF) from a VCF file.

    Parameters
    ----------
    vcf_path : str
        Path to the VCF file to query


    Returns
    -------
    df : pandas.DataFrame
        Table containing unique variants with corresponding FILTER, DP and VAF
    
    Raises
    ------
    FileNotFoundError
        If vcf_path is None, empty or not found
    """

    # Raise an error if file is not found
    if not vcf_path:
        raise FileNotFoundError(f"VCF file not identified in {vcf_path}")

    # Run the bcftools command
    try:
        result = subprocess.run(
            [
                "bcftools",
                "query",
                "-f",
                "'%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%VF]\t%FILTER\n'",
                vcf_path
            ],
            capture_output=True,
            text=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] bcftools query failed: {e.stderr}")
    
    lines = result.stdout.strip().split("\n'")
    # Filter empty lines
    data = [line.split('\t') for line in lines if line]
    df = pd.DataFrame(data, columns=['CHROM', 'POS', 'REF', 'ALT', 'DP', 'VAF', 'FILTER'])
    df = df.replace("'", "", regex=True)
    
    # Convert POS column into a numeric column
    df['POS'] = pd.to_numeric(df['POS'], errors='raise')
    
    return df.drop_duplicates(subset=['CHROM','POS', 'REF', 'ALT'])


def merge_records_to_query_truth(records_df, query_data_df, truth_data_df):
    """
    Merge records dataframe with the data from query and truth vcfs

    Parameters
    ----------
    records_df : pd.DataFrame
        Dataframe obtained from the output of parse_features_file()
    query_data_df : pd.DataFrame
        Query dataframe obtained from the output of extract_variant_data_df()
        when the query_vcf is provided
    truth_data_df : pd.DataFrame
        Query dataframe obtained from the output of extract_variant_data_df()
        when the truth_vcf is provided
    
    Returns
    -------
    merged_df : pd.DataFrame
        Merged dataframe with query and truth data, NAs are handled on FILTER
        DP and VAF columns
    """

    # Specify the columns to merge
    merge_cols = ['CHROM', 'POS', 'REF', 'ALT']
    
    # Common function to rename the columns
    def rename_columns(df, exclude_cols, string_suffix):
        return df.rename(columns={
            col: f"{col}_{string_suffix}" for col in df.columns if col not in exclude_cols})
    
    # Rename all other columns from the query_df with '_query'
    query_data_df_renamed = rename_columns(query_data_df, merge_cols, "query")
    
    # Rename all other columns from the truth_df with '_truth'
    truth_data_df_renamed = rename_columns(truth_data_df, merge_cols, "truth")
    
    # Merge all dataframes
    merge_df = records_df.merge(query_data_df_renamed, on=merge_cols, how='left')
    merge_df = merge_df.merge(truth_data_df_renamed, on=merge_cols, how='left')

    
    # Manage the NaN on some of the columns if they exist
    nan_handling = {
        'FILTER_query': '.',
        'FILTER_truth': '.',
        'DP_query': 0,
        'DP_truth': 0,
        'VAF_query':0.0,
        'VAF_truth':0
        }
    for col, default_val in nan_handling.items():
        if col in merge_df.columns:
            merge_df[col] = merge_df[col].fillna(default_val)
    
    
    return merge_df


def parse_args() -> argparse.Namespace:
    
    parser = argparse.ArgumentParser(
        description="Extract VAF and DP from VCF files based on features CSV.")
    parser.add_argument("--features_file", help="Input features CSV file")
    parser.add_argument("--query_file", help="Input query VCF file")
    parser.add_argument("--truth_file", help="Input truth VCF file")
    parser.add_argument("--output_file", help="Filename to storing output")
    parser.add_argument("--sample_name", help="Sample identifier")
    
    args = parser.parse_args()
    
    return args

def main():
    """
    Main function to extract VAF and DP data from VCF files based on features
    CSV.
    """
    args = parse_args()
    
    # Read the features.csv file
    records = read_features_file(args.features_file, args.sample_name)
    
    # Declare the query and truth vcf data and sore vcf files in a pd.DataFrame
    query_vcf_data = extract_variant_data_df(args.query_file)
    truth_vcf_data = extract_variant_data_df(args.truth_file)
    
    # Merge records and vcf into a dataframe
    merged_df = merge_records_to_query_truth(records, query_vcf_data,
                                             truth_vcf_data)
    
    # Save merged_df
    merged_df.to_csv(args.output_file, index=False)
    print(f"File output stored in {args.output_file}")


if __name__ == "__main__":
    main()
