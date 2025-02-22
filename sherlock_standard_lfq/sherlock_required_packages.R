# title: "sherlock_lfq_pipeline"
# author: "Dain Brademan"
# date: "2025-01-07"
# output: html_document

####################
###### README ######
####################
# 
# This script is a truncated version of the Huttenhain lab's abundance proteomics pipeline (`lfq_pipeline.rmd`) and is meant to run from top-to-bottom on Sherlock, Stanford's high-throughput computing cluster.
# Once you run this script, the Feature- and ProteinLevelData.csv outputs can be fed into the normal lfq pipeline for pairwise comparisons or other analyses.
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

# • Install all required libraries through 
################################
###### Required Arguments ######
################################

# • Specify a project name

### Usage

This pipeline supports proteomics data searched using Spectronaut or MaxQuant. Use this pipeline only for abundance-based proteomics experiments. If you have PTM data, you will be better off using the PTM pipeline. You will need to following files in order to use this pipeline:
  
  **Spectronaut Search**
  
  -   Spectronaut Report (HuttenhainDIA_MsStats_FragmentIon)

-   This file is in .tsv format by default

**MaxQuant Search**
  
  -   Maxquant Evidence File (.txt)

-   Maxquant Protein Groups File  (.txt)

-   Maxquant Annotations File (manually built)

# Installing / Loading Required R Libraries

When I'm working on a script from the ground up, the first thing I like to do is install (if not already installed) & load all the libraries that I'm planning on using before I even start thinking about loading and manipulating any data. R libraries only need to be installed once via the *install.packages("Function")*
  
  ### Installing CRAN Packages
  
  Libraries only need to be installed once. After that, they're saved to your computer and you can view them under the **Packages** tab in the bottom right panel on your screen. How do you tell which libraries/packages you need to install? I usually just run the load statements below, and R will tell me if I don't have a package installed

The packages used in this code come from two different repositories, $\underline{C}omprehensive\space\underline{R}\space\underline{A}rchive\space\underline{N}etwork$ (CRAN) and BioConductor. They're installed differently.

CRAN packages can be installed several different ways. The first way is demonstrated in the code chunk below. Alternatively, you can click to the **Packages** tab over on the bottom right of the screen. You'll find an **Install** and **Update** button there that you can use to browse for packages

```{r}
# Don't need run this every time. 
# Only if you need to install a package.
# Just replace the Package in the parentheses with whatever you need installed
# If you get complaints about Rtools is required to build R packages, you can install it from the following URL: https://cran.r-project.org/bin/windows/Rtools/rtools42/rtools.html
# You don't really need to install it, but warnings are annoying

# install.packages("devtools")
# install.packages("data.table")
# install.packages("dplyr")
# install.packages("magrittr")
# install.packages("R.utils")
# install.packages("ggplot2")
# install.packages("ggrepel")
```

### Installing Bioconductor Packages

BioConductor has its own package manager. Before you can install anything from it, you first need to install the BioConductor manager from CRAN.

```{r}
# install.packages("BiocManager")
```

Then you can install all the BioConductor packages you want! Note, you may be prompted in the console window (bottom left) if you want to update packages. It's generally a good idea. Sometimes a package refuses to update. If that happens, close out of all your R tabs, clear your workspace, and try again. Sometimes incompatible packages can clash and cause really weird errors that are tough to troubleshoot, so I'd recommend keeping things as updated as possible.

```{r}
# Install a specific version of BioConductor core packages
# Don't use this unless you know what you're doing...<BiocManager::install(version = "3.18)>

# Install latest version of BioConductor core packages
# BiocManager::install()

```

Some BioConductor packages don't install with the core distribution. Do install these add-ons, you can specifically request these packages be installed

```{r}
# Install a set of specific non-core package. Format: c("Package1", "Package2", etc....)
# BiocManager::install(c("MSstats", "ComplexHeatmap"))
```

### Load Packages

Now that in theory everything is installed. Let's load our packages. If you get any errors regarding **Package not found**, just go back and install it! A general rule is the fully lowercase libraries are from CRAN and the CamelCase libraries are from BioConductor


# Load Packages
#library(devtools)
#library (data.table)
library (magrittr)
#library(dplyr)
#library(tidyr)
library(stringr)
#library(R.utils) # Actually a cran package

#library(UniProt.ws)
library(purrr)

# Plotting packages
#library (ggplot2)
#library(ggrepel)
rotate.x.axis.text <- theme(axis.text.x = element_text(angle = 90,vjust = 0.5, hjust=1))

# MSstats Packages
#library(MSstats)

# EnrichR Complex Heatmap Package
#library (ComplexHeatmap)

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


# Pipeline Parameters

Enter in values specific to your dataset in this block and run all blocks of code sequentially.

```{r}
# Data Type
# SP - Spectronaut
# MQ - MaxQuant
export.type <- "SP" 

# Appended to all exported data/figures
data.name <- "Example_LFQ"

# create subdirectories to store intermediate tables and figures
Create.Pipeline.Directories(data.name)

# Pipeline considers samples with the same "Replicate/BioReplicate" column to be biological replicates i.e. from the same patient
# This is not true for a typical APEX experiment where we instead have technical replicates
# If `is.case.control` = TRUE, pipeline will make replicate column unique by concatenating Condition & Replicate columns
# If FALSE, pipeline will leave replicate columns as-is. 
# This will impact the statistical model either way.
is.case.control = TRUE

# File paths to your mass spec data in either MaxQuant or Spectronaut format
if (export.type == "SP") {
  # Spectronaut specific parameters
  spectronaut.lfq.report <- "../example_data/example_data_spectronaut.tsv"
} else if (export.type == "MQ") {
  maxquant.evidence <- "../example_data/example_data_MQ_evidence.txt"
  maxquant.proteingroups <- "../example_data/example_data_MQ_proteinGroups.txt"
  maxquant.annotations <- "../example_data/example_data_MQ_annotations.csv"
} else {
  stop("Non-supported data type. Make sure data type is SP or MQ, otherwise pipeline will break.")
}

# Names of all possible contrasts will be generated from Condition names (formatted "ConditionName1-ConditionName2")
## Example: regexContrasts <- c("Fsk-DMSO", "H2O2_Fsk-Fsk", "H2O2-DMSO", "H2O2_Fsk-H2O2")
## This list should be strings that match your exact condition names that you entered in Spectronaut. 
## You can also use RegEex terms if you know what you're doing, but I find it more exact to use full condition names to avoid accidental     RegEx matches

global.regexContrasts <- c(
  # 0 minute APEX vs no biotinylation control
  "0m_40-ctrl_40", "0m_50-ctrl_50", "0m_60-ctrl_60", "0m_70-ctrl_70", "0m_80-ctrl_80",
  
  # 5 minute APEX vs no biotinylation control
  "5m_40-ctrl_40", "5m_50-ctrl_50", "5m_60-ctrl_60", "5m_70-ctrl_70", "5m_80-ctrl_80",
  
  # 5 minute APEX vs 0 minute APEX                          
  "5m_40-0m_40", "5m_50-0m_50", "5m_60-0m_60", "5m_70-0m_70", "5m_80-0m_80"
)


# Sometimes you want to remove certain runs because they're bad, or before data normalization. 
# You can also do this at the condition level for ease, but you'll need to modify the code in the post-normalization section to do this
# Specify the runs here. 
runs.to.remove <- c()
# Uncomment the below line to remove the "control" samples before median normalization
# runs.to.remove <- c("TT000243_AM04-01","TT000244_AM04-02","TT000245_AM04-03","TT000246_AM04-04","TT000247_AM04-05","TT000254_AM04-12","TT000255_AM04-13","TT000256_AM04-14","TT000257_AM04-15","TT000258_AM04-16")
```

# Import Data

This reads the proteomics search engine data verbatim and loads it into a data frame. You can see this data frame in the environment window. (default location: upper left)
```{r}

# Spectronaut specific data files
if (export.type == "SP") {
  spectronaut.lfq.report <- fread(spectronaut.lfq.report)
  
  # MaxQuant specific data files
} else if (export.type == "MQ") {
  maxquant.evidence <- fread(maxquant.evidence)
  maxquant.proteingroups <- fread(maxquant.proteingroups)
  maxquant.annotations <- fread(maxquant.annotations)
}
```

### Clean data

It's recommended to include commonly observed contaminant proteins when processing your proteomic datasets. However, Spectronaut does not remove these entries for you. Fortunately, all contaminant protein accessions are prepended with the text **Cont\_** so we can go through each row in our dataset and toss that data ourselves.

```{r}
if (export.type == "SP") {
  # remove all contaminant proteins.
  spectronaut.lfq.report <- spectronaut.lfq.report %>%
    filter(!grepl("Cont_|NaN|TRYP_PIG", PG.ProteinGroups))
}
```

### Make Replicate column unique if set

This prevents MSstats from treating experimental/technical replicates as bioreplicates. Bioreplicates are modeled differently within the statistical framework.

```{r}
if (is.case.control) {
  spectronaut.lfq.report$R.Replicate <- paste(spectronaut.lfq.report$R.Condition, spectronaut.lfq.report$R.Replicate, sep = ".")
}
```

### Reformat data for msStats Data Processing

Most proteomic pipelines by default do not export data in the format that msStat expects for protein grouping and normalization. Fortunately for us, msStats already has a built-in converter function that we can take advantage of with minimal effort. While it's not the worst thing to reformat this data ourselves, work smarter, not harder.

After this stage, we will no longer need to write special code to handle Spectronaut vs MaxQuant data

```{r}
if (export.type == "SP") {
  
  # the optional arguments are recommended in the msStats user guide. Feel free to change if you know better.
  global.prepared <- SpectronauttoMSstatsFormat(spectronaut.lfq.report,
                                                filter_with_Qvalue = TRUE, ## same as default
                                                qvalue_cutoff = 0.01, ## same as default
                                                removeProtein_with1Feature = FALSE,
                                                use_log_file = FALSE
  )
  
  rm(spectronaut.lfq.report)
  
  Save.Csv.With.Timestamp(global.prepared, "CleanedPreprocessedData.csv", paste(data.name, "data", sep = "_"))
  
} else if (export.type == "MQ") {
  
  global.prepared <- MaxQtoMSstatsFormat(maxquant.evidence, 
                                         maxquant.annotations, 
                                         maxquant.proteingroups)
  
  rm(maxquant.evidence, maxquant.annotations, maxquant.proteingroups)
}

```

### Remove any unwanted runs

This tosses any files listed in the *runs.to.remove* vector. I use this block of code to toss all my control samples before median normalization. 
```{r}
# Tosses any rows where the 'Run' column matches a value in the 'runs.to.remove' vector listed at the top
global.prepared <- global.prepared[!(global.prepared$Run %in% runs.to.remove), ]

```

# Quality Control Plots of the Raw Data

Now is a good time to take a look at the raw data.

### View a distribution of the peptide ion transitions.

This should be monomodal after the *SpectronauttoMSstatsFormat* filtering step

```{r}

ggplot(global.prepared, aes(x=log10(Intensity))) +
  geom_histogram(binwidth = 0.05, fill = "gray", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Peptide Transition Intensities", x = "Log10(Intensity)", y = "Frequences") +
  theme_bw()

ggsave(paste(paste(data.name, "figures", sep = "_"), "Histogram - Peptide Transition Intensities.pdf", sep = "/"))

```

### Count unique peptide sequences per experimental condition

For a general abundance proteomics APEX experiment, we would hope for 70,000+ peptide sequences per replicate. This particular dataset is lower than that due to various reasons, but what we are looking for are consistent numbers across related replicates and the APEX samples. APEX control samples (no construct or no biotin so no enrichment) hopefully will be lower.

```{r}
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

ggsave(paste(paste(data.name, "figures", sep = "_"), "Histogram - Peptides per Sample.pdf", sep = "/"))
```

### Get boxplot of quantitative values per sample

Just the same as the number of peptide sequences, we also expect related samples to have similar quantitative distributions. The APEX experimental workflow involved biotin labeling of proximal proteins and subsequent enrichment. As the control samples were not designed to undergo biotin labeling, the control samples ideally will have *lower* intensity distributions than the APEX samples.

```{r}

ggplot(data = global.prepared, aes (x = BioReplicate,
                                    y = log10(Intensity),
                                    fill = Condition)) +
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "lightgray") +
  labs(title = "Peptide Transition Intensities per Sample",
       x = "Sample",
       y = "Log10 Transition Intensity") +
  rotate.x.axis.text

ggsave(paste(paste(data.name, "figures", sep = "_"), "Peptide Transition Intensities per Sample.pdf", sep = "/"))

```

# Use MSstats for protein summarization

This is where we summarize peptides (or peptide transitions depending on your data export format) back into protein groups. Depending on your experiment, you're going to want to turn normalization on or off. For your typical APEX experiment, the control condition doesn't have Biotin Phenol meaning biotin enrichment does not occur. We expect an entirely different population of proteins to be enriched for in the BP samples, so its best to not do any sort of normalization at this point.

```{r}
global.proteinSummarization = MSstats::dataProcess(global.prepared,
                                                   normalization = 'FALSE',	# "equalizeMedians"	once you remove non-BP controls for APEX samples	
                                                   logTrans = 2,				
                                                   featureSubset = 'highQuality',				
                                                   summaryMethod="TMP",
                                                   censoredInt='NA',				
                                                   MBimpute=FALSE,				
                                                   maxQuantileforCensored=0.999)

# Merge in protein names
global.proteinSummarization$ProteinLevelData$gene.name <- map.gene.names(global.proteinSummarization$ProteinLevelData$Protein)

# Write out data to CSVs
Save.Csv.With.Timestamp(global.proteinSummarization$FeatureLevelData, "FeatureLevelData.csv", paste(data.name, "data", sep = "_"))
Save.Csv.With.Timestamp(global.proteinSummarization$ProteinLevelData, "ProteinLevelData.csv", paste(data.name, "data", sep = "_"))
```

# (Optional) Reload DataProcess results

If you ever want to go back & reanalyze data after the data process function, it can be frustrating to have to run through the entire pipeline again up to this point.

```{r}
# global.proteinSummarization <- list(ProteinLevelData = NULL, FeatureLevelData = NULL)

# global.proteinSummarization$ProteinLevelData <- fread(paste0(paste(data.name, "data/", sep = "_"), "20240326_ProteinLevelData.csv"))
# global.proteinSummarization$FeatureLevelData <- fread(paste0(paste(data.name, "data/", sep = "_"), "20240326_FeatureLevelData.csv"))

# global.proteinSummarization$ProteinLevelData$RUN <- factor(global.proteinSummarization$ProteinLevelData$RUN)
# global.proteinSummarization$ProteinLevelData$Protein <- factor(global.proteinSummarization$ProteinLevelData$Protein)
# global.proteinSummarization$ProteinLevelData$GROUP <- factor(global.proteinSummarization$ProteinLevelData$GROUP)

```

# Do Principal Component Analysis
look to see if there are serious batch effects, outlying samples, etc.
```{r}

# Go from long to wide format
# Row names are protein groups
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

# Optionally save the figure
ggsave(paste(paste(data.name, "figures", sep = "_"), "Principal Component Analysis - Condition.pdf", sep = "/"))

```

# Plot protein intensities pre-normalization

All of our peptide-level data has now been summarized into protein groups. Let's inspect the intensity distributions of this data. Note: the data process function already log2 transformed this data, so we don't have to ourselves.

```{r}

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

ggsave(paste(paste(data.name, "figures", sep = "_"), "Protein Abundances per Sample.pdf", sep = "/"))

```

```{r}

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

ggsave(paste(paste(data.name, "figures", sep = "_"), "UniqueProteinGroupsPerCondition.pdf", sep = "/"))
```

# Plot average protein CVs per condition

```{r}

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

ggsave(paste(paste(data.name, "figures", sep = "_"), "Protein CVs per Condition.pdf", sep = "/"))

```

# Native Carboxylase Enrichment

This is a quality control plot I borrowed from Ben Polacco.If you've normalized your data, it's useful. Without normalization the results are pretty difficult to interpret.

*APEX works by labeling neighboring/interacting proteins with biotin which is then used to purify labeled proteins. There are also proteins that are endogenously biotinylated which will co-purify with the APEX-labeled proteins. Here we look at a subset of these endogenous biotin proteins, and we inspect their post-normalization background levels, which are inversely related to the labeling-efficiency of APEX. More background after normalization implies there is less APEX-labeled signal.*
  
  If you see your control conditions at higher abundance compared to the rest of your samples, great, APEX labeling and enrichment is behaving as we expected!
  
  ```{r}

biotin.carboxylases.up <- c("O00763","P05165","P11498","Q13085","Q96RQ3")

native.biotin.data <- as.data.table(global.proteinSummarization$ProteinLevelData)
native.biotin.data <- native.biotin.data[grepl(paste(biotin.carboxylases.up, collapse = "|"), native.biotin.data$Protein), ]


# Plotting code
ggplot(native.biotin.data, aes(x = Protein, y = LogIntensities, fill = GROUP)) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "bar",
    position = position_dodge(width = 0.8),  # Adjust the width as needed
    width = 0.6  # Adjust the width of bars
  ) +
  stat_summary(
    fun.data = mean_sdl,
    fun.args = list(mult = 1),
    geom = "errorbar",
    position = position_dodge(width = 0.8),
    width = 0.4  # Adjust the width of error bars
  ) +
  labs(title = "Mean Log Intensities with Standard Deviation by GROUP",
       x = "GROUP", y = "Log Intensities") +
  theme_bw()

ggsave(paste(paste(data.name, "figures", sep = "_"), "NativeCarboxylases.pdf", sep = "/"))

rm(biotin.carboxylases.up, test)
```

# Do global pairwise comparisons to check for successful APEX enrichment

Another way you can check for successful APEX enrichment is to make volcano plots comparing protein abundances between non-BP vs BP conditions. What you should see is a mass up-regulation of most proteins. If you don't see that, your enrichment likely didn't work, or you have a ton of background binding. To make these comparisons, we can use the GroupComparison function from msStats. To use this function, you need the output from MSstats::DataProcess, and you also have to build a contrast matrix, essentially a numerical vector indicating what experimental conditions you want to compare.

I adapted Ben Polacco's old *makeContrast* function to help build the group comparison matrix without errors. Again, you can do this yourself if you're careful without using this function.

```{r}

# makes an MsStats-compatible contrast matrix using the condition comparisons in the parameter block
contrast.matrix <- Make.LFQ.Contrast.Matrix (input.data.frame = global.proteinSummarization, condition.vector = global.regexContrasts)
# Model-based comparison
# Do pairwise comparison via MSstats using the provide contrast matrix, then select the comparison results

global.pairwiseComparison <- MSstats::groupComparison(contrast.matrix, global.proteinSummarization)$ComparisonResult

global.pairwiseComparison$gene.name <- map.gene.names(global.pairwiseComparison$Protein)

Save.Csv.With.Timestamp(global.pairwiseComparison, "GroupComparisonsData.csv", paste(data.name, "data", sep = "_"))

```

# Volcano Plots

Do all the volcano plots for all experiment contrasts

```{r}

# This is essentially a foreach loop. 
# Takes each of the labels of the above pairwise comparison
# Subselects that specific data
# Finally, makes the volcano plot
lapply(global.regexContrasts, function(comparison) {
  
  pvalue.threshold <- 0.05
  log2FC.threshold <- 1 # Ruth likes log2(1.5)
  
  # define significant proteins
  thisPairwiseComparison <- as.data.table(global.pairwiseComparison %>% filter(Label == comparison))
  
  
  ## This chunks adds a new column to our pairwise comparison data table
  # For all proteins with |Log2FC| < 1 and pval > 0.05, value is "Not"
  # For proteins with |Log2FC| >= 1 and pval <= 0.05, value is "Up" or "Down" as appropriate
  thisPairwiseComparison[, Significance := "Not"]
  thisPairwiseComparison[pvalue < pvalue.threshold & abs(log2FC) > log2FC.threshold,
                         Significance := ifelse (log2FC > 0, "Up", "Down")]
  
  thisCondition1 <- strsplit(comparison, '-')[[1]][1]
  thisCondition2 <- strsplit(comparison, '-')[[1]][2]
  
  ## Render volcano plots.
  plot <- ggplot(thisPairwiseComparison, aes(x = log2FC, y = -log10(pvalue), color = Significance, label = gene.name)) +
    
    # circles representing proteins
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c(Not = "gray", Down = "#67a9cf", Up = "#ef8a62")) +
    ggrepel::geom_text_repel(data = thisPairwiseComparison[Significance != "Not"], size = 2, max.overlaps = 20) +
    labs(title = paste("Differential Protein Abundance Between", thisCondition1, "&", thisCondition2), x = paste("Protein Log2 Fold-change"), y = "-Log10(p-value)") +
    
    
    # significance labels and lines
    # vertical lines
    geom_vline(xintercept = c(-log2FC.threshold, log2FC.threshold), linetype = "dashed", color = "darkgray") +
    annotate("text", x = c(-1.1, 1.1), y = 0, label = c("log2FC = -1", "log2FC = 1"), vjust = 0, hjust = c(1, 0), color = "black") +
    
    # horizontal line
    geom_hline(yintercept = -log10(pvalue.threshold), linetype = "dashed", color = "darkgray") +
    annotate("text", x = Inf, y = -log10(log2FC.threshold + 0.02), label = "pval = 0.05", vjust = 1, hjust = 1, color = "black") +
    theme_bw()
  
  print(plot)
  
  ggsave(paste(paste(data.name, "figures", sep = "_"), paste("Protein - Volcano - ", comparison, ".pdf", sep = ""), sep = "/"))
})

```