#merging all csv data within all nested folders into one file in the directory\\Clingen\cg\Genetics\WGS\WGS_Excel\Sendout CSVs

import os
import glob
import pandas as pd

#Path
BASE_PATH = r"\\Clingen\cg\Genetics\WGS\WGS_Excel\Sendout CSVs"
#OUT_PATH = need to be determined

def merge_genetics_csv(base_path,output_path):
    all_dfs = []
    reference_headers = None
    
    # Search for all CSV files
    search_path = os.path.join(base_path, "**", "*.csv")
    csv_files = glob.glob(search_path, recursive=True)
    
    if not csv_files:
        print(f"No CSV files found in {base_path}.")
        return

    print(f"Found {len(csv_files)} CSV files. Validating headers...")

    for file_path in csv_files:
        file_name = os.path.basename(file_path)
        
        # Extract the 'file_batch' from the file name after 'eme_'
        if "_eme_" in file_name:
            file_batch = file_name.split("_eme_")[1].replace(".csv", "")
        else:
            print(f"Warning: Skipping {file_name} because it doesn't contain '_eme_'")
            continue
            
        try:
            # Efficiently read JUST the header row (nrows=0) to check names/length
            header_check_df = pd.read_csv(file_path, nrows=0)
            current_headers = list(header_check_df.columns)
            
            # Establish the baseline headers from the very first file
            if reference_headers is None:
                reference_headers = current_headers
            
            # Compare current headers to the reference headers
            if current_headers != reference_headers:
                print(f"\n[ERROR] Header mismatch detected in file: {file_name}")
                print(f"Expected ({len(reference_headers)} cols): {reference_headers}")
                print(f"Found ({len(current_headers)} cols): {current_headers}")
                print("Skipping this file to avoid misalignment.")
                continue  # Skips this file and moves to the next one
            
            # If valid, read the full data content
            df = pd.read_csv(file_path)
            
            # Append the file_batch data column
            df['file_batch'] = file_batch
            all_dfs.append(df)
            
        except Exception as e:
            print(f"Error processing file {file_name}: {e}")
            continue

    # Merge everything together
    if all_dfs:
        merged_df = pd.concat(all_dfs, ignore_index=True)
        
        # Save the master file
        output_filename = "Merged_WGS_CSV_data.csv"
        output_path = os.path.join(output_path, output_filename)
        merged_df.to_csv(output_path, index=False)
        
        print("\n--- Success! ---")
        print(f"Successfully processed and merged {len(all_dfs)} files.")
        print(f"Master file saved at: {output_path}")
    else:
        print("No valid data was merged.")

if __name__ == "__main__":
    merge_genetics_csv(BASE_PATH, OUT_PATH)
