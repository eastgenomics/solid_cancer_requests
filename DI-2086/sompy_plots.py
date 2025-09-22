#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 31 15:22:14 2025

@author: arun
"""

import dxpy
import matplotlib.pyplot as plt
import pandas as pd
from io import StringIO

dxfile = dxpy.DXFile('file-J24FXZQ48fFz1gb8J1Jxg1Vq').read()
type(dxfile)
df = pd.read_csv(StringIO(dxfile))

# Select the project
project = "project-J1zbgv84kZ75y8B1vzXv6X2v"

# Pick up the list of files of interest, for now the stats.csv should do the trick
stats_files = list(dxpy.find_data_objects(
    classname="file",
    name="*.stats.csv",
    name_mode="glob",
    project=project,
    folder="/query_filtered_VCFs_removed_sompy_results"
))

# For each sample
shared_variants = []
unique_to_truth = []
unique_to_query = []
for file in stats_files:
    print(file)
    dxfile = dxpy.DXFile(file['id']).read()
    df = pd.read_csv(StringIO(dxfile))
    # Count the number of shared variants (i.e. tp)
    shared_variants.append(df[df['type'] == 'records']['tp'].iloc[0])
    # Count the number of variants found only in query (fp)
    unique_to_query.append(df[df['type'] == 'records']['fp'].iloc[0])
    # Count the munber of variants found only in truth (fn)
    unique_to_truth.append(df[df['type'] == 'records']['fn'].iloc[0])

# Make a database with the columns (sample, shared, truth only, query only)
stats_df = pd.DataFrame(
    [dxpy.describe(file['id'], fields={'name'}) for file in stats_files])
stats_df['sample_name'] = stats_df['name'].str.extract(r'^(\d+-\d+[SQK]\d+)')
stats_df['run'] = stats_df['name'].str.extract(r'(25TSOD\d{2})')
stats_df['shared'] = shared_variants
stats_df['truth only'] = unique_to_truth
stats_df['query only'] = unique_to_query


# Use the database to create a plot for each sequencing run
runs = stats_df['run'].unique()

for i, run in enumerate(runs):
    run_data = stats_df[stats_df['run'] == run].reset_index(drop=True)

    # Adjust figure size based on number of samples
    if len(run_data) > 25:
        fig_width = 16
    elif len(run_data) > 15:
        fig_width = 12
    else:
        fig_width = 8

    # Create a new figure for each run
    fig, ax = plt.subplots(figsize=(fig_width, 8))

    # Sort data by shared variants for better visualization
    run_data = run_data.sort_values('shared', ascending=False).reset_index(drop=True)

    # Create stacked bars
    p1 = ax.bar(range(len(run_data)), run_data['shared'],
                label='Shared', color='steelblue')
    p2 = ax.bar(range(len(run_data)), run_data['truth only'],
                bottom=run_data['shared'], label='truth only', color='salmon')
    p3 = ax.bar(range(len(run_data)), run_data['query only'],
                bottom=run_data['shared'] + run_data['truth only'],
                label='query only', color='darkseagreen')
    
    ax.set_ylabel('Number of variants', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    ax.set_title(f'Number of variants per sample from {run} (n={len(run_data)})')

    # Adjust x-axis labels based on number of samples
    ax.set_xticks(range(len(run_data)))
    if len(run_data) > 15:
        ax.set_xticklabels(run_data['sample_name'], rotation=90, fontsize=8)
    else:
        # Show all sample names for smaller runs
        ax.set_xticklabels(run_data['sample_name'], rotation=45, fontsize=8, ha='right')

    # Add grid for better readability
    ax.grid(axis='y', alpha=0.3, linestyle='--')

    # Show legend on each plot (plots are not stored)
    ax.legend(loc='upper right', framealpha=0.9)

    plt.tight_layout()
    plt.show()

stats_df.to_csv('sompy_plots_table.csv')
