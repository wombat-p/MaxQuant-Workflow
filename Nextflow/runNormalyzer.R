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
write.table(final_exp, "full_design.txt",sep="\t", row.names=F,quote=F)

# comparison set to everything versus first
if (comps == "") {
   print("setting statistical tests to compare all conditions versus the first")
   comps <- unique(final_exp[,"group"])
   comps <- paste0(comps[2:length(comps)],"-", comps[1])  
} else {
  comps <- unlist(strsplit(comps,","))
}

# run Normalyzer
NormalyzerDE::normalyzer(jobName="Normalyzer", designPath="full_design.txt", dataPath="protein_file.txt", zeroToNA = TRUE, inputFormat = "maxquantprot", outputDir="./",sampleColName="Experiment",requireReplicates=F)
print("Now running differential expression analysis")
NormalyzerDE::normalyzerDE(jobName="Normalyzer", comparisons=comps, designPath="full_design.txt", dataPath=paste0("./Normalyzer/",normalyzerMethod,"-normalized.txt"), outputDir="./", sampleCol="Experiment")
