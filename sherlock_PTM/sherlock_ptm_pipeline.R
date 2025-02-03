# title: "sherlock_lfq_pipeline"
# author: "Dain Brademan"
# date: "2025-01-24"


####################
###### README ######
####################
# 
# This script is a truncated version of the Huttenhain lab's phopshoproteomics pipeline (`phospho_pipeline_MsStats+GlobalProteome.rmd`) and is meant to run from top-to-bottom on Sherlock, Stanford's high-throughput computing cluster.
# Once you run this script, the Feature- and ProteinLevelData.csv output for the phospho can be fed into the end of the phospho pipeline for pairwise comparisons or other analyses.
#   
# Running this file will produce two directories:
#   • [NAME]_data:
#     • CleanedPreprocessedData.csv - Cleaned/preprocessed transition-level data for low-level QC
#     • FeatureLevelData.csv - Compiled transition-level data for pairwise comparisons
#     • ProteinLevelData.csv - Compiled protein-level data for pairwise comparisons
#   • [NAME]_figures:
#     • Various quality control figures
#     • Principal component analysis


###########################
###### Prerequisites ######
###########################
#
# • Install all required packages before attempting to run this file. 
#   • Refer to the README for a shortlist
# • Upload the Spectronaut export file to Oak
# • Configure an SBATCH file and run
# 

# Load Packages
library(devtools)
library (data.table)

library (magrittr)
library(dplyr)
library(tidyr)
library(stringr)

library(R.utils) # Actually a cran package

library(UniProt.ws)
library(purrr)

# Plotting packages
library (ggplot2)
library(ggrepel)
rotate.x.axis.text <- theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))

# MSstats Packages
library(MSstats)
library(MSstatsPTM)

# Krogan Lab Utility Functions- Load From GitHub
utils <- "https://raw.githubusercontent.com/HuttenhainLab/bp_utils/master/"
source(paste0(utils, "enrichmentTestFunctions.R"))
source(paste0(utils,"MSstats_V4_Functions.R"))

# Dain Utility Functions - Load From GitHub
#drb_utils <- "C:/Users/dainb/OneDrive/Documents/GitHub/drb_utils/"
drb_utils <- "https://raw.githubusercontent.com/HuttenhainLab/drb_utils/main/"
source(paste0(drb_utils,"data_management.R"))
source(paste0(drb_utils,"spectronaut_data_cleaning.R"))
source(paste0(drb_utils,"gene_ontology_enrichment.R"))

# Define functions to bind gene names to a data table based on UniProt ID column
map.gene.names <- function(protein.column, taxonomy.id = 9606) {
  unique.phosphosites <- unique(unlist(strsplit(as.character(protein.column), ";")))
  UniProt.IDs <- unique(sapply(unique.phosphosites, function(x) strsplit(x, "_")[[1]][1]))
  
  UniProt.mapping <- select(UniProt.ws(taxId = 9606), # 10090 = mouse
                            keys = UniProt.IDs,
                            columns = c("UniProtKB", "gene_primary"),
                            keytype = "UniProtKB")
  
  # Create a named vector for fast lookup
  gene.name.map <- setNames(UniProt.mapping$Gene.Names..primary., UniProt.mapping$From)
  
  # Function to replace UniProt IDs with gene names
  map_protein <- function(entry) {
    # Split the entry by semicolons
    
    phosphosites <- unlist(strsplit(entry, ";"))
    
    # Replace UniProt IDs with gene names
    mapped_sites <- sapply(phosphosites, function(site) {
      # Split the site by underscore to separate UniProt ID and the rest
      parts <- unlist(strsplit(site, "_"))
      UniProtID <- parts[1]
      
      # Check if the UniProt ID has a mapped gene name
      gene_name <- gene.name.map[UniProtID]
      
      # If a gene name is found, use it; otherwise, keep the UniProt ID
      if (!is.na(gene_name)) {
        # Recombine the gene name with the remaining part
        paste(c(gene_name, parts[-1]), collapse = "_")
      } else {
        # If no mapping exists, return the original site
        site
      }
    })
    
    # Rejoin the mapped sites with semicolons
    paste(mapped_sites, collapse = ";")
  }
  
  # Apply the mapping function to each row of the protein column
  mapped_column <- sapply(as.character(protein.column), map_protein)
  
  return(mapped_column)
}



################################
######  Argument Parsing  ######
################################

input.arguments <- commandArgs(trailingOnly = TRUE)

# Prepended onto output directories
project.name <- input.arguments[1]

# Relative or absolute path to Spectronaut report
path.to.spectronaut.report <- input.arguments[2] 
path.to.fasta <- input.arguments[3]
path.to.annotation.file <- input.arguments[4]


# If `is.case.control` = TRUE, makes replicate column unique
is.case.control <- as.logical(input.arguments[5])

# Normalize? 
normalization.enabled <- as.logical(input.arguments[6])

# If fifth argument specified, throw away samples matching the REGEX
if (length(input.arguments) == 7) {
  sample.removal.substring <- input.arguments[7]
}

# Nothing below this needs editing
#
#
#
######################################
###### Begin Automated Pipeline ######
######################################

# create subdirectories to store intermediate tables and figures
Create.Pipeline.Directories(project.name)

# Import Data
spectronaut.lfq.report <- as.data.table(read.csv(path.to.spectronaut.report, sep = "\t", stringsAsFactors = FALSE))

# Toss contaminants
spectronaut.lfq.report <- spectronaut.lfq.report %>%
  filter(!grepl("Cont_|NaN|TRYP_PIG", PG.ProteinGroups))


### Make Replicate column unique if set
if (is.case.control) {
  spectronaut.lfq.report$R.Replicate <- paste(spectronaut.lfq.report$R.Condition, spectronaut.lfq.report$R.Replicate, sep = ".")
}

### Reformat data for msStats Data Processing
# the optional arguments are recommended in the msStats user guide. Feel free to change if you know better.
global.prepared <- SpectronauttoMSstatsFormat(spectronaut.lfq.report,
                                              filter_with_Qvalue = TRUE, ## same as default
                                              qvalue_cutoff = 0.01, ## same as default
                                              removeProtein_with1Feature = FALSE,
                                              use_log_file = FALSE
)

rm(spectronaut.lfq.report)
Save.Csv.With.Timestamp(global.prepared, "CleanedPreprocessedData.csv", paste(project.name, "data", sep = "_"))

# try and save memory
gc()

### Remove any unwanted runs
# Tosses any rows where the 'Run' column matches a value in the 'runs.to.remove' vector listed at the top
if (exists(sample.removal.substring)) {
  global.prepared <- global.prepared[!grepl(sample.removal.substring, global.prepared$Condition), ]
}

# Quality Control Plots of the Raw Data
### View a distribution of the peptide ion transitions.
ggplot(global.prepared, aes(x=log10(Intensity))) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Peptide Transition Intensities", x = "Log10(Intensity)", y = "Frequences") +
  theme_bw()

ggsave(paste(paste(project.name, "figures", sep = "_"), "Histogram - Peptide Transition Intensities.pdf", sep = "/"))


### Count peptides per sample
unique_counts <- as.data.table(global.prepared) %>%
  filter(!is.na(Intensity)) %>%
  filter(Intensity != 0) %>%
  group_by(BioReplicate) %>%
  summarize(Unique_Peptide_Count = n_distinct(PeptideSequence))

ggplot(unique_counts, aes(x = BioReplicate, y = Unique_Peptide_Count)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_text(aes(label = Unique_Peptide_Count), vjust = -0.5, color = "black", size = 3) +
  labs(title = "Unique Peptide Sequences per Replicate",
       x = "Replicate",
       y = "Number of Peptide Sequences per Replicate") +
  theme_bw() +
  rotate.x.axis.text

ggsave(paste(paste(project.name, "figures", sep = "_"), "Histogram - Peptides per Sample.pdf", sep = "/"))


### View transition intensity distributions per sample
ggplot(data = global.prepared, aes (x = BioReplicate,
                                    y = log10(Intensity),
                                    fill = Condition)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "lightgray") +
  labs(title = "Peptide Transition Intensities per Sample",
       x = "Sample",
       y = "Log10 Transition Intensity") +
  rotate.x.axis.text

ggsave(paste(paste(project.name, "figures", sep = "_"), "Peptide Transition Intensities per Sample.pdf", sep = "/"))

# Use MSstats for protein summarization

if (normalization.enabled == TRUE) {
  normalization <- "equalizeMedians"
} else {
  normalization = "FALSE"
}

global.proteinSummarization = MSstats::dataProcess(global.prepared,
                                                   normalization = normalization,
                                                   logTrans = 2,				
                                                   featureSubset = 'highQuality',				
                                                   summaryMethod="TMP",
                                                   censoredInt='NA',				
                                                   MBimpute=FALSE,				
                                                   maxQuantileforCensored=0.999)

# Merge in protein names
global.proteinSummarization$ProteinLevelData$gene.name <- map.gene.names(global.proteinSummarization$ProteinLevelData$Protein)

# Write out data to CSVs
Save.Csv.With.Timestamp(global.proteinSummarization$FeatureLevelData, "FeatureLevelData.csv", paste(project.name, "data", sep = "_"))
Save.Csv.With.Timestamp(global.proteinSummarization$ProteinLevelData, "ProteinLevelData.csv", paste(project.name, "data", sep = "_"))


# Do Principal Component Analysis

# Columns are each sample's LogIntensity values for each protein
intensity.matrix <- as.matrix(dcast(as.data.table(global.proteinSummarization$ProteinLevelData), Protein ~ SUBJECT, value.var = "LogIntensities"),
                              rownames = "Protein")

complete.data.matrix <- intensity.matrix[complete.cases(intensity.matrix),]

# Do PCA
pcaOut <- prcomp(t(complete.data.matrix))

# Reshape from matrix to data table
pcaDT <- as.data.table(pcaOut$x, keep.rownames=TRUE)

# Add column named `mainGroup`. 
# This value is derived from the rowname by splitting on '_' and grabbing the first element
pcaDT[, Condition := tstrsplit(rn, "\\.")[[1]]] #transpose & split

# Add column named `batch`. 
# This value is derived from the rowname by splitting on '_' and grabbing the first element
pcaDT[, Batch := sapply(strsplit(rn, "\\."), function(x) tail(x, 1))] #transpose & split

pcaPercentVar <- round(100 * (pcaOut$sdev^2)/sum(pcaOut$sdev^2), 1)

# Visualize PC1 and PC2, color data points by batch value.
# You can also color data points by mainGroup or other factors that you scrape out above.
p <- ggplot (pcaDT, aes(x=PC1, y=PC2, color = Condition )) + 
  geom_point(alpha=1.0, size=4) + 
  ggrepel::geom_text_repel(aes(label=rn), show.legend = FALSE, size = 3) +
  theme_bw() + 
  xlab (sprintf ("PC1, %.1f%%", pcaPercentVar[1])) + 
  ylab (sprintf ("PC2, %.1f%%", pcaPercentVar[2])) + 
  ggtitle (sprintf ("PCA using %d Proteins (log intensity)", nrow(complete.data.matrix))) 
p

ggsave(paste(paste(project.name, "figures", sep = "_"), "Principal Component Analysis - Condition.pdf", sep = "/"))

# Plot protein intensities pre-normalization
ggplot (global.proteinSummarization$ProteinLevelData, 
        aes (x = interaction (SUBJECT), 
             y = LogIntensities, 
             fill = GROUP)
) + 
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "lightgray") +
  labs(title = "Protein Abundance Distributions",
       x = "Replicate",
       y = "Log2 Protein Abundance Intensity") +
  rotate.x.axis.text

ggsave(paste(paste(project.name, "figures", sep = "_"), "Protein Abundances per Sample.pdf", sep = "/"))

# Count proteins per condition
unique_counts <- global.proteinSummarization$ProteinLevelData %>%
  filter(LogIntensities != 0) %>%
  filter(!is.na(LogIntensities)) %>%
  filter(!is.infinite(LogIntensities)) %>%
  group_by(GROUP) %>%
  summarize(Unique_Protein_Count = n_distinct(Protein))

ggplot(unique_counts, aes(x = GROUP, y = Unique_Protein_Count)) +
  geom_bar(stat = "identity", fill = "skyblue", width = 0.5) +
  geom_text(aes(label = Unique_Protein_Count), vjust = -0.5, color = "black", size = 3) +
  labs(title = "Count of Protein Groups per Condition",
       x = "Condition",
       y = "Number of Protein Groups") +
  rotate.x.axis.text

ggsave(paste(paste(project.name, "figures", sep = "_"), "UniqueProteinGroupsPerCondition.pdf", sep = "/"))


# Calculate Protein CVs
proteomics_data <- global.proteinSummarization$ProteinLevelData %>%
  group_by(GROUP, Protein) %>%
  summarize(CV = sd(LogIntensities) / mean(LogIntensities) * 100)

# Plot the distribution of CV for each group
ggplot(proteomics_data, aes(x = GROUP, y = CV, fill = GROUP)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "lightgray") +
  labs(title = "Protein Coefficient of Variation Distribution by Group",
       x = "Please note CV's were calculated in Log2 space, not raw intensities",
       y = "Coefficient of Variation (%)") +
  rotate.x.axis.text

ggsave(paste(paste(project.name, "figures", sep = "_"), "Protein CVs per Condition.pdf", sep = "/"))