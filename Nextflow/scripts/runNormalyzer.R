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
A<-read.csv("exp_file.tsv",sep="\t")
B<- read.csv("exp_file2.tsv", sep="\t")

B[,1] <- tools::file_path_sans_ext(basename(B[,1]))
A$raw_file <- tools::file_path_sans_ext(basename(A$raw_file))
names(A) <- c("file", "group")
B$Experiment <- make.names(B$Experiment)
C<-merge(A,B,by.x="file",by.y="Name")
final_exp <- unique(C[,c("Experiment","group")])
num_cond <- nrow(final_exp)
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
if (min(table(final_exp[,"group"])) > 1) {
   NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.tsv", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.tsv", dataPath="peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   print("Now running differential expression analysis")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerProteins", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerProteins/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerPeptides", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerPeptides/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
} else {
  NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.tsv", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.tsv", dataPath="peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  print("No statistical testing as at least one sample group with only 1 replicate")
  write.csv(NA,"NormalyzerProteins/Normalyzer_stats.tsv")
  write.csv(NA,"NormalyzerPeptides/Normalyzer_stats.tsv")
}

## Preparing for standardized format
# Reading files
peptides <- read.csv("peptide_file.txt", sep="\t", row.names=1)
proteins <- read.csv("protein_file.txt", sep="\t", row.names=1)
norm_peptides <- read.csv(paste0("NormalyzerPeptides/", normalyzerMethod, "-normalized.txt"), sep="\t", row.names = 1)
norm_proteins <- read.csv(paste0("NormalyzerProteins/", normalyzerMethod, "-normalized.txt"), sep="\t", row.names = 1)

# changing column names
peptides$missed_cleavages <- peptides$Missed.cleavages
peptides$charge <- peptides$Charges
peptides$protein_group <- peptides$Proteins

colnames(peptides) <- unlist(sub("^Experiment\\.", "number_of_psms_", colnames(peptides)))

# colnames(peptides) <- unlist(sub("^LFQ\\.intensity\\.","abundance_", colnames(peptides)))
# for (i in 1:nrow(final_exp)) {
#   colnames(peptides) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(peptides)))
#   colnames(proteins) <- unlist(sub(final_exp[i,1], final_exp[i,2], colnames(proteins)))
# }
colnames(proteins) <- unlist(sub("^Razor\\.\\.\\.unique\\.peptides\\.", "number_of_peptides_", colnames(proteins)))
proteins$protein_group <- rownames(proteins)
for (s in final_exp$Experiment) {
  colnames(norm_proteins) <- sub(s, paste0("abundance_",s), colnames(norm_proteins))
  colnames(norm_peptides) <- sub(s, paste0("abundance_",s), colnames(norm_peptides))
}

# getting relevant columns
proteins <- cbind(proteins[rownames(norm_proteins), c("protein_group", grep("^number_of_peptides_", colnames(proteins), value=T))],
                  norm_proteins[,grep("^abundance_", colnames(norm_proteins), value=T)])
peptides <- cbind(modified_sequence=rownames(norm_peptides), peptides[rownames(norm_peptides), c("protein_group", grep("^number_of_psms_", colnames(peptides), value=T))],
                  norm_peptides[,grep("^abundance_", colnames(norm_peptides), value=T)])

write.csv(proteins, "stand_prot_quant_merged.csv", row.names = F)
write.csv(peptides, "stand_pep_quant_merged.csv", row.names = F)

cat("Done\n")