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

# Define paths and parameters
ml ruse

R_SCRIPT="sherlock_lfq_pipeline.R"
PROJECT_NAME="MOR04"
INPUT_FILE="./20250106_093125_MOR04_Report.tsv"
IS_CASE_CONTROL="TRUE"
NORMALIZATION_ENABLED="TRUE"
REMOVE_THESE_SAMPLES="ctrl"

# Run the R script with four arguments
ruse Rscript $R_SCRIPT $PROJECT_NAME $INPUT_FILE $IS_CASE_CONTROL $NORMALIZATION_ENABLED $REMOVE_THESE_SAMPLES