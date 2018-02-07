problem.dir <- "labels/H3K4me3_TDH_immune/samples/tcell/McGill0004/problems/chr10:60000-17974675"
data.name <- gsub("[-/:]", "_", gsub("labels/|problems/|samples/", "", problem.dir))
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
out.file <- paste0(coverage.bedGraph, ".out")
system(paste("grep data_i=", out.file), intern=TRUE)
cmd.vec <- system(paste("grep ^PeakSegFPOP", out.file), intern=TRUE)
library(namedCapture)
pen.str <- str_match_named(cmd.vec[length(cmd.vec)], "penalty=(?<penalty>.*).db")
library(data.table)
coverage <- fread(coverage.bedGraph)
setnames(coverage, c("chrom","chromStart","chromEnd","count"))

library(coseg)
fit <- PeakSegFPOPchrom(coverage, as.numeric(pen.str))

labels.bed <- file.path(problem.dir, "labels.bed")
labels <- fread(labels.bed)
setnames(labels, c("chrom","chromStart","chromEnd","annotation"))

assign(data.name, list(
  coverage=data.frame(coverage),
  ##penalty=pen.str,
  labels=data.frame(labels)))
data.path <- file.path("~/lib64/R/library/cosegData/data", paste0(data.name, ".RData"))
save(list=data.name, file=data.path, compress="xz")
cat("scp thocking@guillimin.hpc.mcgill.ca:", data.path, " ~/R/cosegData/data", sep="")

error.file.vec <- system("grep code labels/*/*/problems/*/coverage.bedGraph.out|grep -v 0$|sed 's/:Exit.*//'", intern=TRUE)
recent.vec <- system("grep 'Begin PBS Epilogue Thu Sep 29' labels/*/*/problems/*/coverage.bedGraph.out", intern=TRUE)
recent.out.vec <- sub(":Begin.*", "", recent.vec)
system(paste("grep code", paste(recent.out.vec, collapse=" "), "|grep 1$|sed 's/out:Exit.*/err/'|xargs tail"))

sh.file.vec <- sub("out$", "sh", error.file.vec)
for(sh.file in sh.file.vec){
  cmd <- paste("qsub", sh.file)
  system(cmd)
}
