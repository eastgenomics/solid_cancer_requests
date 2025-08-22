import csv
import os
import glob
import subprocess
import argparse


def parse_features_file(features_file):
    """
    Parse a features CSV file and extract variant information.

    Parameters
    ----------
    features_file : str
        Path to the input CSV file containing variant features from sompy.

    Returns
    -------
    tuple
        A tuple containing:
        - sample_name : str
            Sample name extracted from the filename
        - records : list of dict
            List of dictionaries containing variant information with keys:
            'sample', 'CHROM', 'POS', 'REF', 'ALT' and all original CSV columns

    Examples
    --------
    >>> sample_name, records = parse_features_file('sample123.features.csv')
    >>> print(sample_name)
    'sample123'
    >>> print(len(records))
    100
    """
    sample_name = os.path.basename(features_file).split(".")[0]
    print(f"Parsing features file: {features_file} for sample: {sample_name}")
    records = []

    with open(features_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            chrom = row["CHROM"]
            pos = row["POS"]
            
            # Use REF.truth if REF is empty or missing
            ref = row.get("REF").strip()
            if not ref or ref == "":
                ref = row.get("REF.truth").strip()
            
            # Use ALT.truth if ALT is empty or missing
            alt = row.get("ALT").strip()
            if not alt or alt == "":
                alt = row.get("ALT.truth").strip()
            
            # Skip rows where we still don't have ref/alt after trying fallbacks
            if not ref or not alt:
                print(f"Warning: Skipping row at {chrom}:{pos} - missing REF/ALT data")
                continue
            row_copy = row.copy()
            row_copy["REF"] = ref
            row_copy["ALT"] = alt
            row_copy["POS"] = pos
            row_copy["sample"] = sample_name
            
            records.append(row_copy)

    return sample_name, records


def find_vcf_files(sample_name):
    """
    Find S1 and S2 VCF files for a given sample name.

    Parameters
    ----------
    sample_name : str
        Full sample name used to construct the base name for file matching.
    Expected format: "instrument_id-SP-additional-info"

    Returns
    -------
    tuple
        A tuple containing:
        - s2_vcf : str or None
            Path to the S2 VCF file if found, None otherwise
        - s1_vcf : str or None
            Path to the S1 VCF file if found, None otherwise

    Notes
    -----
    The function uses the first two hyphen-separated parts of the sample name
    to create a base name for globbing VCF files in the 'VCFs/' directory.
    S2 files are identified by containing 's2' in the filename.

    Examples
    --------
    >>> s2_vcf, s1_vcf = find_vcf_files('137157064-25122S0075-25TSOD55-8471')
    >>> print(s2_vcf)
    'VCFs/137157064-25122S0075-25TSOD55s2-8471_processed.vcf.gz'
    >>> print(s1_vcf)
    'VCFs/137157064-25122S0075-25TSOD55-8471_processed.vcf.gz'
    """
    base_name = "-".join(sample_name.split("-")[:2])
    print(f"Finding VCF files for base name: {base_name}")
    s2_vcf = glob.glob(f"query_filtered_VCFs/S2_VCF/{base_name}*s2*.vcf.gz")
    s1_vcf = glob.glob(f"query_filtered_VCFs/S1_VCF/{base_name}*.vcf.gz")

    # Remove s2 from s1 list
    s1_vcf = [f for f in s1_vcf if "s2" not in f]

    return s2_vcf[0] if s2_vcf else None, s1_vcf[0] if s1_vcf else None


def extract_variant_data(vcf_path, chrom, pos, ref, alt):
    """
    Extract variant depth (DP) and variant allele frequency (VAF) from a VCF file.

    Parameters
    ----------
    vcf_path : str or None
        Path to the VCF file to query
    chrom : str
        Chromosome identifier (e.g., '1', 'chr1', 'X')
    pos : int or str
        Genomic position of the variant
    ref : str
        Reference allele sequence
    alt : str
        Alternative allele sequence

    Returns
    -------
    tuple
        A tuple containing:
        - dp : int
            Depth of coverage at the variant position (0 if not found)
        - vaf : float
            Variant allele frequency (0.0 if not found)
        - filter : str
            Variant FILTER field ('.' if not found)

    Raises
    ------
    FileNotFoundError
        If vcf_path is None or empty

    Notes
    -----
    Uses bcftools query to extract variant information. If the variant is not
    found or bcftools fails, returns (0, 0.0) as default values.
    """
    pos = int(pos)
    print(f"Extracting data for {chrom}:{pos} {ref} -> {alt} from VCF: {vcf_path}")
    # use bcftools with subprocess
    if not vcf_path:
        raise FileNotFoundError(f"VCF file not found for {chrom}:{pos} {ref} -> {alt}")
    try:
        result = subprocess.run(
            [
                "bcftools",
                "query",
                "-r",
                f"{chrom}:{pos}-{pos}",
                "-f",
                "'%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%DP\t[%VF]\n'",
                vcf_path
            ],
            capture_output=True,
            text=True,
            check=True,
        )
        for line in result.stdout.strip().split("\n"):
            if not line or line == "'":
                continue
            # Split the line into components
            c, p, r, a, filter_string, dp, af = line.strip().split("\t")
            #   chr, pos, ref, alt, dp, af for this line
            print(f"Extracted data: {c}:{p} {r} -> {a}: DP={dp}, VAF={af}, FILTER={filter_string}")
            # Check if the reference and alt match
            # Note: a can be multiple alleles separated by commas
            if r == ref and alt in a.split(","):
                return int(dp) if dp != "." else 0, float(af) if af != "." else 0.0, filter_string
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] bcftools query failed: {e.stderr}")
    return 0, 0.0, '.'


def write_output(out_path, records, s1_data, s2_data):
    """
    Write extracted variant data to a CSV file.

    Parameters
    ----------
    out_path : str
        Path to the output CSV file
    records : list of dict
        List of variant records from the original features file
    s1_data : dict
        Dictionary mapping variant keys to S1 data tuples (DP, VAF, FILTER)
        Keys are tuples of (CHROM, POS, REF, ALT)
    s2_data : dict
        Dictionary mapping variant keys to S2 data tuples (DP, VAF, FILTER)
        Keys are tuples of (CHROM, POS, REF, ALT)

    Returns
    -------
    None
    """
    if not records:
        print("No records to write.")
        return

    # Get fieldnames from the first record and add new columns
    fieldnames = list(records[0].keys()) + ["DP_s2", "VAF_s2", "FILTER_s2", "DP_s1", "VAF_s1", "FILTER_s1"]

    print(f"Writing output to {out_path} with {len(records)} records...")

    with open(out_path, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for rec in records:
            # Create key from record data - convert POS to string to match your key format
            key = (str(rec["CHROM"]), str(rec["POS"]), rec["REF"], rec["ALT"])

            # Get s2 data (DP, VAF, INFO) or defaults
            s2_dp, s2_vaf, s2_filter = s2_data.get(key, (0, 0.0, '.'))

            # Get s1 data (DP, VAF, INFO) or defaults
            s1_dp, s1_vaf, s1_filter = s1_data.get(key, (0, 0.0, '.'))

            # Create output row with all original sompy fields plus new data
            row = {
                **rec,  # All original fields from the record
                "DP_s2": s2_dp,
                "VAF_s2": s2_vaf,
                "FILTER_s2": s2_filter,
                "DP_s1": s1_dp,
                "VAF_s1": s1_vaf,
                "FILTER_s1": s1_filter
            }

            writer.writerow(row)

    print(f"Successfully wrote {len(records)} records to {out_path}")


def main(features_file, output_file):
    """
    Main function to extract VAF and DP data from VCF files based on features CSV.
    """
    sample_name, records = parse_features_file(features_file)
    s2_vcf, s1_vcf = find_vcf_files(sample_name)

    print(f"Using VCFs:\ns2: {s2_vcf}\ns1: {s1_vcf}")

    s2_data = {}
    s1_data = {}

    for rec in records:
        s2_data_row = extract_variant_data(
            s2_vcf, rec["CHROM"], rec["POS"], rec["REF"], rec["ALT"]
        )
        s1_data_row = extract_variant_data(
            s1_vcf, rec["CHROM"], rec["POS"], rec["REF"], rec["ALT"]
        )
        key = (rec["CHROM"], rec["POS"], rec["REF"], rec["ALT"])
        print(f"Processing {key}...")
        print(f"Data for {key}: s2={s2_data_row}, s1={s1_data_row}")
        s2_data[key] = s2_data_row
        s1_data[key] = s1_data_row
        print(
            f"Extracted data for {rec['CHROM']}:{rec['POS']} {rec['REF']} -> {rec['ALT']}: "
            f"DP_s2={s2_data_row[0]}, VAF_s2={s2_data_row[1]}, FILTER_s2={s2_data_row[2]}, "
            f"DP_s1={s1_data_row[0]}, VAF_s1={s1_data_row[1]}, FILTER_s2={s1_data_row[2]}"
        )
    print(f"Extracted data for {len(records)} records.")
    write_output(output_file, records, s1_data, s2_data)
    print(f"Output written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract VAF and DP from VCF files based on features CSV."
    )
    parser.add_argument("--features_file", help="Input features CSV file")
    parser.add_argument("--output_file", help="Output CSV file")
    args = parser.parse_args()
    main(args.features_file, args.output_file)
