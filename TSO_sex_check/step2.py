#!/usr/bin/env python

import sys
import csv
import math
import dxpy as dx
import pandas as pd


output = []
project=sys.argv[1]

# find output files from eggd_sex_check
files = list(dx.find_data_objects(
name="*_idxstat.tsv", name_mode="glob",
project=project,
describe=True))

for file in files:
    sample = file['describe']['name'].rstrip('_idxstat.tsv')

    # get chr1 and chrY read counts from eggd_sex_check output
    chr1 = chrY = 0
    with dx.open_dxfile(file['id']) as file:
        reader = csv.reader(file, delimiter="\t")
        for line in reader:
            if line[0] == "chr1":
                chr1 = int(line[2])
            elif line[0] == "chrY":
                chrY = int(line[2])

    # calculate inverse log of ratio between read counts
    norm_chr_y = chrY / chr1 if chr1 != 0 else 0
    epsilon = 1e-9  # small value to avoid log(0)
    ratio = -math.log(norm_chr_y + epsilon)

    output.append({
        'sample': sample,
        'chr1': chr1,
        'chrY': chrY,
        'ratio': ratio
        })

output_df = pd.DataFrame(output)
output_fn = f"{project}_grouped_outputs.tsv"
output_df.to_csv(output_fn, sep='\t', index=False)
