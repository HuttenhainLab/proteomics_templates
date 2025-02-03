#!/usr/bin/bash
#SBATCH --job-name=lfq_pipeline
#SBATCH --output=%j.lfq_pipeline.out
#SBATCH --error=%j.lfq_pipeline.err
#SBATCH --time=4-00:00:00
#SBATCH --partition=normal
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=224G

# Required Arguments
	# These must be defined or the script will crash
r_script="./sherlock_ptm_pipeline.R"
project_name="TEST01"
spectronaut_report="./20250106_093125_MOR04_Report.tsv"
sample_annotations="./sample_annotations.csv"
protein_database="./HomoSapien_Canonical_2024_10_26.fasta"

# Optional Arguments
	# These arguments have assumed defaults in the script if not set
	# Define them if you need different behavior
is_case_control="TRUE"
normalization_mode="equalizeMedians"
sample_removal_substring="ctrl"

# Run the R script with all arguments
	# Remember double quotes if your file names have spaces
	# I've placed each argument on its own line for readability
Rscript "$r_script" \
	--project_name "$project_name" \
	--spectronaut_report "$spectronaut_report" \
	--sample_annotations "$sample_annotations" \
	--protein_database "$protein_database" \
	--is_case_control "$is_case_control" \
	--normalization_mode "$NORMALIZATION_ENABLED" \
	--sample_removal_substring "$REMOVE_THESE_SAMPLES"