#!/usr/bin/bash

project=$1
echo $project

# Generate a batch input file for each project
dx select "$project"
dx generate_batch_inputs \
	-i input_bam='(.*).bam$' \
	-i index_file='(.*).bam.bai' \
	-o "$project"
# select testing project and create per-sample output folders
# dx select project-J48v2jj4FBjyB9kZ0gk248j0  
dx mkdir -p "/output/sample_swap_investigation"

# run eggd_sex_check for each batch
dx run eggd_sex_check \
	--batch-tsv "${project}.0000.tsv" \
	-ifemale_threshold=6.44 \
	-imale_threshold=5.16 \
	--destination "/output/sample_swap_investigation/" \
	-y

