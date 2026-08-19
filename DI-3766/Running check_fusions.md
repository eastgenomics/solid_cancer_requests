## Running check_fusions.py ##

check_fusions.py will compare expected fusions to the fusions called for samples in a provided workbook.
The script will extract the StarAligner and Arriba.

If there are multiple fusions detected for the same gene, only the first one will be returned.

Inputs: 

- DNAnexus file ID corresponding to a fusion workbook
- A TSV containing Cases (SP) and fusions detected
 
outputs:

- A TSV file with the following columns: sample_id, fusion_name, Arriba_breakpoint1, Arriba_breakpoint2, Arriba_split_reads1, Arriba_split_reads2, Arriba_coverage1, Arriba_coverage2, Arriba_confidence, StarFusion_LeftBreakpoint, StarFusion_RightBreakpoint, StarFusion_JunctionReadCount

Note: if a fusion is not detected by one of the tools, the corresponding rows contain "."


## Set up ##

Identify the workbook of interest. This will be held on DNAnexus and will be the output of eggd_generate_fusion_workbook

Identify the expected fusions. The expected fusions were linked in an excel file in an email chain.
This is not being added here as one of the sheets includes PID. Use this file to create the expected fusions files by following the process:
- Identify which sheet corresponds to the run of interest (ie sheet 26-RESVal4 corresponds to run 260731_A01295_0820_BHNYC3DRX7_RES)
- Copy the 'Cases (SP)' and 'Fusion detected' columns into a local csv file
- remove 'SP-' from sample names
- reformat fusions to be Gene1::Gene2
- If there are multiple fusions for a sample, split the fusions onto multiple lines
- If there is text included (ie "missed fusion xyz") remove the additional text
- If there are single gene variants (ie rearrangements, ITD, ex skipping) remove addition additional text, any fusions where this gene is Gene1 or Gene2 will be returned 

## Example ## 
```
 python check_fusions.py \
 --workbook_file_id file-J9jKpfj4yXj1KZ3GK8bvj5Vx \
 --expected_fusions_tsv 260731_A01295_0820_BHNYC3DRX7_RES_expected.tsv 

 ```