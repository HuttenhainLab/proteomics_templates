#!/usr/bin/bash
#SBATCH --job-name=lfq_pipeline
#SBATCH --output=%j.lfq_pipeline.out
#SBATCH --error=%j.lfq_pipeline.err
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --mem=16G

# Load the ruse module
ml ruse

# Define paths and parameters

R_SCRIPT="sherlock_lfq_pipeline.R"
PROJECT_NAME="MOR04"
INPUT_FILE="./example_data.csv"
IS_CASE_CONTROL="TRUE"
NORMALIZATION_ENABLED="TRUE"
REMOVE_THESE_SAMPLES="ctrl"

# Run the R script with provided arguments
ruse Rscript $R_SCRIPT $PROJECT_NAME $INPUT_FILE $IS_CASE_CONTROL $NORMALIZATION_ENABLED $REMOVE_THESE_SAMPLES