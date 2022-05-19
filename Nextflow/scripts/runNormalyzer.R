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

B[,1] <- sub(".raw","",B[,1])
names(B) <- c("file","group")
A$Experiment <- make.names(A$Experiment);
C<-merge(A,B,by.x="Name",by.y="file");
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

# run Normalyzer
if (min(table(final_exp[,"group"])) > 1) {
   NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.tsv", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.tsv", dataPath="peptide_file.txt", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F)
   print("Now running differential expression analysis")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerProteins", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerProteins/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
   NormalyzerDE::normalyzerDE(jobName="NormalyzerPeptides", comparisons=comps, designPath="Normalyzer_design.tsv", dataPath=paste0("./NormalyzerPeptides/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment", leastRepCount="0")
} else {
  NormalyzerDE::normalyzer(jobName="NormalyzerProteins", designPath="Normalyzer_design.txt", dataPath="protein_file.tsv", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  NormalyzerDE::normalyzer(jobName="NormalyzerPeptides", designPath="Normalyzer_design.txt", dataPath="peptide_file.tsv", zeroToNA = TRUE, inputFormat = "maxquantpep", outputDir="./",sampleColName="Experiment",requireReplicates=F,skipAnalysis=T)
  print("No statistical testing as at least one sample group with only 1 replicate")
  write.csv(NA,"NormalyzerProteins/Normalyzer_stats.tsv")
  write.csv(NA,"NormalyzerPeptides/Normalyzer_stats.tsv")
}
