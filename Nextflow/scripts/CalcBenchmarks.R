###### Script to calculate general benchmarks
## For more details about the metrics, see
## https://docs.google.com/spreadsheets/d/1tH08Sc78h3oHyRoWgjyqoL4r_943PaDvQwUbjldn2fo/edit#gid=728728122


library(stringr)
library(matrixStats)
library(jsonlite)

cat("\n#### Calculating benchmarks ###\n")

Benchmarks <- NULL

# numPeptides=0, numProteins=0, dynRangePep=0, propUniquePep=0, uniqueStrippedPep=0, percMissingPep=0,
# aucDiffRegPeptides=list(), tprPep0.01=list(), tprPep0.05=list(), tFDRPep0.01=list(), tFDRPep0.05=list(), 
# propMisCleavedPeps=list(),sumSquareDiffFCPep=0, sdWithinRepsPep=0, skewnessPeps=0, kurtosisPeps=0, sdPeps=0,
# # Protein level
# numQuantProtGroups=0, dynRangeProt=0, propUniqueProts=0, percMissingProt=0, meanPepPerProt=0, aucDiffRegProteins=list(), 
# tFDRProt0.01=list(), tFDRProt0.05=list(), tprProt0.01=list(), tprProt0.05=list(), sumSquareDiffFCProt=0, sdWithinRepsProt=0, propMisCleavedProts=0,
# propDiffRegWrongIDProt0.01=list(),propDiffRegWrongIDProt0.05=list(),skewnessProts=0, kurtosisProts=0, sdProts=0,
# # PTM level
# numProteoforms=0, numModPeptides=0, meanProteoformsPerProt=0, propModAndUnmodPep=0, aucDiffRegAdjModPep=list(),
# tFDRAdjModPep0.01=list(), tFDRAdjModPep0.05=list(), tprAdjModPep0.01=list(), tprAdjModPep0.05=list(),
# propDiffRegPepWrong0.01=list(),propDiffRegPepWrong0.05=list(), percOverlapModPepProt=0, sumSquareDiffFCModPep=0)

cat("\n### Reading files ###\n")
## TODO: change to statistics output
StatsPep <- read.csv("stand_pep_quant_merged.csv")
StatsProt <- read.csv("stand_prot_quant_merged.csv")
ExpDesign <- read.csv("exp_design.txt", sep="\t")
Params <- read_json("params.json")
  
#### Functionality
Functionality = list()
### Traceability
Traceability = list() 
## Spectra
Spectra = list() 

# How do we check that? Would require looking into intermediate results files
Spectra[["TraceableSpectra"]] <- F

# Set to FALSE for now as determining them would be tricky
Spectra[["UniversalSpectumIdentifiers"]] <- F

# Check for corresponding columns and missing value in them
id_cols <- grep("SpectraID|ScanNumber",colnames(StatsPep))
Spectra[["PeptideToSpectra"]] <- ifelse(length(id_cols > 0), ifelse(any(rowSums(!is.na(StatsPep[, id_cols, drop=F])) == 0), F, T), F)

id_cols <- grep("SpectraID|ScanNumber",colnames(StatsProt))
Spectra[["ProteinToSpectra"]] <- ifelse(length(id_cols > 0), ifelse(any(rowSums(!is.na(StatsProt[, id_cols, drop=F])) == 0), F, T), F)

Traceability[["Spectra"]] <- Spectra

## FileNames
FileNames = list()

# TODO: check for sdrf or experimental design file
FileNames[["ResultsToRawFiles"]] <- NA

# URL available anywhere? TODO check for sdrf
FileNames[["PublicRawFiles"]]  <- NA

Traceability[["FileNames"]] <- FileNames

## Parameters
Parameters=list()

# Check for settings in parameter file (yaml or sdrf): TODO: related to sdrf and/or yaml
Parameters[["Settings"]] > NA

# TODO Check for availability of non-default experimental design file (TODO standardize within workflows)
Parameters[["ExperimentalDesign"]] <- F

Traceability[["Parameters"]] <- Parameters
Functionality[["Traceability"]] <- Traceability

### Reproducibility
Reproducibility <- list()

## Files
Files <- list()

# Reproducible results. Has been run in a container? TODO check for call with docker or singularity
Files[["Identify"]]  <- NA
Reproducibility[["Files"]] <- Files

### Performance
Performance=list()

##  Identification
Identification <- list()
# How many PSMs? TODO: needs columns (or only one) with PSM numbers
Identification[["PSMNumber"]] <- NA

Identification[["PeptideNumber"]] <- nrow(StatsPep)

# Uniquely identified proteins. Assuming that protein accessions of protein groups are separated via ";"
Identification[["ProteinNumber"]] <- sum(!grepl(";",StatsProt$protein_group))

Identification[["ProteinGroupNumber"]] <- nrow(StatsProt)

# percentage of peptides identified in all samples
all_pep_samples <- grep("^abundance_", colnames(StatsPep))
Identification[["PeptideCoverage"]] <- sum(rowSums(is.na(StatsPep[,all_pep_samples,drop=F])) == 0)

# precentage of proteins identified in all samples
all_prot_samples <- grep("^abundance_", colnames(StatsProt))
Identification[["ProteinCoverage"]] <- sum(rowSums(is.na(StatsProt[,all_prot_samples,drop=F])) == 0)

# distribution of peptides per protein group (only 1-10)
tab <- table(unlist(StatsProt[,grep("^number_of_peptides_", colnames(StatsProt)),drop=F]))

Identification[["PeptidesPerProtein"]] <- as.data.frame(tab[1:10])
Performance[["Identification"]] <- Identification

## Quantification
Quantification=list()

# CV and correlation of peptides within replicates
tPep <- tProt <- tPep2 <- tProt2 <-tprotquant <- tpepquant <-  NULL
for (i in unique(ExpDesign$exp_condition)) {
  tquant <- as.matrix(StatsPep[,grep(paste0("^abundance_", i), colnames(StatsPep)), drop=F])
  tPep <- c(tPep, rowSds(tquant, na.rm=T) / rowMeans(tquant, na.rm=T))
  tPep2 <- c(tPep2, cor(log2(tquant), use="pairwise.complete.obs"))
  tquant <- 2^as.matrix(StatsProt[,grep(paste0("^abundance_", i), colnames(StatsProt)), drop=F])
  tProt <- c(tProt, rowSds(tquant, na.rm=T) / rowMeans(tquant, na.rm=T))
  tProt2 <- c(tProt2, cor((tquant), use="pairwise.complete.obs"))
  
}
tPep2[tPep2 == 1] <- NA
tProt2[tProt2 == 1] <- NA

# Calculate medians
Quantification[["CVPeptides"]] <- median(tPep, na.rm=T)
Quantification[["CVProteins"]] <- median(tProt, na.rm=T)

Quantification[["CorrelationPeptides"]] <- mean(tPep2, na.rm=T)
Quantification[["CorrelationProteins"]] <- mean(tProt2, na.rm=T)

# number of peptides with at least 50% coverage
tpepquant <- StatsPep[, all_pep_samples]
Quantification[["NumberOfPeptides"]] <- sum(rowSums(!is.na(tpepquant)) >= 0.5*ncol(tpepquant))

# 5% top vs. 5% bottom quantile
qs <- quantile(tpepquant, probs = seq(0,1,0.05), na.rm=T)
Quantification[["DynamicPeptideRange"]] <- qs[length(qs)-1] / qs[2]

# number of protein groups with at least 50% coverage
tprotquant <- StatsProt[, all_prot_samples]
Quantification[["NumberOfProteinGroups"]] <- sum(rowSums(!is.na(tprotquant)) >= 0.5*ncol(tprotquant))

# 5% top vs. 5% bottom quantile
qs <- quantile(tprotquant, probs = seq(0,1,0.05), na.rm=T)
Quantification[["DynamicProteinRange"]] <- 2^(qs[length(qs)-1] - qs[2])

Performance[["Quantification"]] <- Quantification

# GroundTruth=list()
#   FoldChangePrecision=NA,
#   AUCs=vector(),
#   FDRs5Perc=vector(),
#   FDRs1Perc=vector(),
#   ProteinLinearity=NA
# ),

## Statistics
Statistics=list()
# Calculate the average per comparison (columns)
tstat <- StatsPep[,grep("^differential_abundance_qvalue", colnames(StatsPep)), drop=F]
Statistics[["DifferentialRegulatedPeptides5Perc"]]  <- colSums(tstat < 0.05, na.rm=T) / ncol(tstat)
Statistics[["DifferentialRegulatedPeptides1Perc"]]  <- colSums(tstat < 0.01, na.rm=T) / ncol(tstat)
tstat <- StatsProt[,grep("^differential_abundance_qvalue", colnames(StatsProt)), drop=F]
Statistics[["DifferentialRegulatedProteins5Perc"]]  <- colSums(tstat < 0.05, na.rm=T) / ncol(tstat)
Statistics[["DifferentialRegulatedProteins1Perc"]]  <- colSums(tstat < 0.01, na.rm=T) / ncol(tstat)
# Percentages in full data
Statistics[["MissingPeptideValues"]] <- sum(is.na(tpepquant)) / length(as.matrix(tpepquant))
Statistics[["MissingProteinValues"]] <- sum(is.na(tprotquant)) / length(as.matrix(tprotquant))
Performance[["Statistics"]] <- Statistics

## Digestion (distribution of miscleavages)
Digestion=list()
Digestion[["Efficiency"]] <- as.data.frame(table(StatsPep$miscleavages))
Performance[["Digestion"]] <- Digestion

## PTMs
PTMs=list()
# percentages of different PTMs
mods <- str_extract_all( StatsPep$modified_peptide, "\\[[a-z,A-Z,0-9,\\.]*\\]|\\([a-z,A-Z]*\\)|\\<[a-z,A-Z]*\\>")
PTMs[["PTMDistribution"]] <- as.data.frame(table(unlist(mods)))
# How many PTMs per peptide
PTMs[["PTMOccupancy"]] <- as.data.frame(table(sapply(mods, length)))
Performance[["PTMs"]] <- PTMs
Functionality[["Performance"]] <- Performance

### Parameter
# Most of them anysways already in the Params file and sdrf!!!
Parameter=list()

## Identification
Identification<-list()
Identification[["DatabaseSize"]] <- as.numeric(system(paste0("grep '>' database.fasta | wc -l"), intern=T))
# Percentage of canonical sequences, assuming that the sequences are given as uniprot sequences
Identification[["CanonicalSequences"]] <- (nrow(StatsProt) - length(grep("\\-", StatsProt$protein_group))) / nrow(StatsProt)
# Whether there is localization in the workflow. TODO: relevant? Needs to be taken from Params
Identification[["PTMLocalization"]] <- NA
# Whether parsimony was used in the identification. TODO: relevant? Needs to be taken from Params
Identification[["Parsimony"]] <- Params$protein_inference
Parameter[["Identification"]] <- Identification

## Quantification
Quantification=list()
# TODO: add to all workflows
Quantification[["Alignment"]] <- Params$match_between_runs
# TODO: add as parameter?
Quantification[["Imputation"]] <- Params$imputation
Quantification[["Normalization"]] <- Params$normalization_method
Parameter[["Quantification"]] <- Quantification
Functionality[["Parameter"]] <- Parameter
Benchmarks[["Functionality"]] <- Functionality
cat("## Done\n")
write_json(Benchmarks, path="benchmarks.json")
