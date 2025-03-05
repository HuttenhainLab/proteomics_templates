#!/usr/bin/bash
#SBATCH --job-name=lfq_pipeline       
#SBATCH --output=%j.lfq_pipeline.out
#SBATCH --error=%j.lfq_pipeline.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=normal
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem=12G

# Make sure R is loaded
ml R

# Define R arguments, file paths, and data processing parameters

R_SCRIPT="sherlock_lfq_pipeline.R"                # Name of R script
PROJECT_NAME="Example_dataset"                              # Output folder names
INPUT_FILE="./example_data_spectronaut.tsv"   # Path to Spec. Report

IS_CASE_CONTROL="TRUE"                            # Make replicate #s unique
                                                    # "TRUE" or "FALSE" only

NORMALIZATION_ENABLED="TRUE"                      # Normalize dataset
                                                    # "TRUE" - median normalization
                                                    # "FALSE" - no normalization

REMOVE_THESE_SAMPLES="ctrl"                       # Removes samples containing substring "ctrl" before normalization

# Run the R script with all specified arguments
Rscript $R_SCRIPT $PROJECT_NAME $INPUT_FILE $IS_CASE_CONTROL $NORMALIZATION_ENABLED $REMOVE_THESE_SAMPLES