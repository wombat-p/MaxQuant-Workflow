# reading cmd arguments
args <- commandArgs(trailingOnly = TRUE)
normalyzerMethod <- strsplit(grep('--method', args, value = TRUE), split = '=')[[1]][[2]]
comps <- strsplit(grep('--comps', args, value = TRUE), split = '=')[[1]]
if (length(comps) > 1) {
   comps <- comps[[2]]
} else {
   comps <- ""
}

# merging experimental design file from sdrf parser with actual design
final_exp<-read.csv("Normalyzer_design.tsv",sep="\t")
# B<- read.csv("exp_file2.tsv", sep="\t")
# 
final_exp$Experiment <- paste0(make.names(final_exp$source_name), "_Tr_", final_exp$technical_replicate)
# C<-merge(A,B,by.x="file",by.y="Name")
# final_exp <- unique(C[,c("Experiment","group")])
# num_cond <- nrow(final_exp)
write.table(final_exp, "Normalyzer_design.tsv",sep="\t", row.names=F,quote=F)

# comparison set to everything versus first
if (comps == "") {
   print("setting statistical tests to compare all conditions versus the first")
   comps <- unique(final_exp[,"group"])
   comps <- paste0(comps[2:length(comps)],"-", comps[1])  
} else {
  comps <- unlist(strsplit(comps,","))
}

## changing peptides to modification string of "sequence_[ptm]_mod.id in peptide input for Normalyzer
# Create some modification string
peps <- read.csv("peptide_file.txt", sep="\t")
cmods <- grep("site\\.IDs$", colnames(peps), value=T)
mods <- unlist(sub("\\.site\\.IDs", "", cmods))
mods <- unlist(gsub("\\.", "", mods))
names(mods) <- cmods
for (c in cmods) {
  peps[,c] <- sapply(peps[,c], function(x) ifelse(is.na(x), "", paste0("_[",mods[c],"]_",unlist(strsplit(as.character(x),";")), collapse="")))
}
peps$Sequence <- paste0(peps$Sequence, peps[,cmods])
write.table(peps, "peptide_file.txt", row.names=F, sep="\t", quote=F)

## run Normalyzer
if (min(table(final_exp[,"group"])) > 1 & length(unique(final_exp[,"group"])) > 1) {
   NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.tsv", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.tsv", dataPath="peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   print("Now running differential expression analysis")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerProteins", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerProteins/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerPeptides", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerPeptides/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
} else {
  NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.tsv", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.tsv", dataPath="peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  print("No statistical testing as at least one sample group with only 1 replicate or only one sample group")
  write.csv(NA,"NormalyzerProteins/Normalyzer_stats.tsv")
  write.csv(NA,"NormalyzerPeptides/Normalyzer_stats.tsv")
}

## Preparing for standardized format
# Reading files
peptides <- read.csv("peptide_file.txt", sep="\t", row.names=1)
proteins <- read.csv("protein_file.txt", sep="\t", row.names=1)
norm_peptides <- read.csv(paste0("NormalyzerPeptides/", normalyzerMethod, "-normalized.txt"), sep="\t", row.names = 1)
norm_proteins <- read.csv(paste0("NormalyzerProteins/", normalyzerMethod, "-normalized.txt"), sep="\t", row.names = 1)
stats_peptides <- stats_proteins <- NULL
if (file.exists("NormalyzerPeptides/NormalyzerPeptides_stats.tsv")) {
  stats_peptides <- read.csv("NormalyzerPeptides/NormalyzerPeptides_stats.tsv", sep="\t", row.names=1)
  stats_proteins <- read.csv("NormalyzerProteins/NormalyzerProteins_stats.tsv", sep="\t", row.names=1)
}

# changing column names
peptides$missed_cleavages <- peptides$Missed.cleavages
peptides$charge <- peptides$Charges
peptides$protein_group <- peptides$Proteins
if (!is.null(stats_peptides)) {
  pval_cols <- colnames(stats_peptides)[grep("AdjPVal$", colnames(stats_peptides))]
  colnames(stats_peptides)[grep("AdjPVal$", colnames(stats_peptides))] <- paste0("differential_abundance_qvalue", sub("_AdjPVal","", pval_cols))
  pval_cols <- colnames(stats_proteins)[grep("AdjPVal$", colnames(stats_proteins))]
  colnames(stats_proteins)[grep("AdjPVal$", colnames(stats_proteins))] <- paste0("differential_abundance_qvalue", sub("_AdjPVal","", pval_cols))
  pval_cols <- colnames(stats_peptides)[grep("PValue$", colnames(stats_peptides))]
  colnames(stats_peptides)[grep("AdjPVal$", colnames(stats_peptides))] <- paste0("differential_abundance_pvalue", sub("_PValue","", pval_cols))
  pval_cols <- colnames(stats_proteins)[grep("PValue$", colnames(stats_proteins))]
  colnames(stats_proteins)[grep("PValue$", colnames(stats_proteins))] <- paste0("differential_abundance_pvalue", sub("_PValue","", pval_cols))
}

colnames(peptides) <- unlist(sub("^Experiment\\.", "number_of_psms_", colnames(peptides)))

# colnames(peptides) <- unlist(sub("^LFQ\\.intensity\\.","abundance_", colnames(peptides)))
# for (i in 1:nrow(final_exp)) {
#   colnames(peptides) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(peptides)))
#   colnames(proteins) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(proteins)))
# }

# changing generic names to experimental design
colnames(proteins) <- unlist(sub("^Razor\\.\\.\\.unique\\.peptides\\.", "number_of_peptides_", colnames(proteins)))
proteins$protein_group <- rownames(proteins)
for (s in 1:nrow(final_exp))  {
  substitute_from <- paste0(make.names(final_exp$source_name[s]),"_Tr_",final_exp$technical_replicate[s])
  substitute_with <- paste0(make.names(final_exp$group[s]),"_",final_exp$technical_replicate[s])
  # avoid setting X in front 
  if (grepl("^X", substitute_with)) substitute_with <- sub("^X","",substitute_with)
  colnames(proteins) <- sub(substitute_from, substitute_with, colnames(proteins))
  colnames(peptides) <- sub(substitute_from, substitute_with, colnames(peptides))
  colnames(norm_proteins) <- sub(paste0("^",substitute_from), paste0("abundance_", substitute_with), colnames(norm_proteins))
  colnames(norm_peptides) <- sub(paste0("^",substitute_from), paste0("abundance_", substitute_with), colnames(norm_peptides))
}

# getting relevant columns
norm_peptides[,grep("^abundance_", colnames(norm_peptides), value=T)] <- 2^(norm_peptides[,grep("^abundance_", colnames(norm_peptides), value=T)])
if (is.null(stats_peptides)) {
  proteins <- cbind(proteins[rownames(norm_proteins), c("protein_group", grep("^number_of_peptides_", colnames(proteins), value=T))], 
                  norm_proteins[,grep("^abundance_", colnames(norm_proteins), value=T)])
  peptides <- cbind(modified_peptide=rownames(norm_peptides), 
                  peptides[rownames(norm_peptides), c("protein_group", grep("^number_of_psms_", colnames(peptides), value=T))],
                  norm_peptides[,grep("^abundance_", colnames(norm_peptides), value=T)])
} else {
  proteins <- cbind(proteins[rownames(norm_proteins), c("protein_group", grep("^number_of_peptides_", colnames(proteins), value=T))], 
                  stats_proteins[rownames(norm_proteins), grep("^differential_abundance", colnames(stats_proteins), value=T)],
                  norm_proteins[,grep("^abundance_", colnames(norm_proteins), value=T)])
  peptides <- cbind(modified_peptide=rownames(norm_peptides), 
                  peptides[rownames(norm_peptides), c("protein_group", grep("^number_of_psms_", colnames(peptides), value=T))],
                  norm_peptides[,grep("^abundance_", colnames(norm_peptides), value=T)],
                  stats_peptides[rownames(norm_peptides), grep("^differential_abundance", colnames(stats_peptides), value=T)])
}



write.csv(proteins, "stand_prot_quant_merged.csv", row.names = F)
write.csv(peptides, "stand_pep_quant_merged.csv", row.names = F)
exp_design_out <- final_exp[, c("Run","group" )]
colnames(exp_design_out) <- c("raw_file","exp_condition")
write.table(exp_design_out, "exp_design_calcb.tsv", quote=F, sep="\t", row.names=F)

cat("Done\n")
