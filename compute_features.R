arg.vec <- "test/H3K36me3_AM_immune_McGill0106_chr16_46385801_88389383"
arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_features.R data_dir/sample_dir/problems/problem_dir")
}
problem.dir <- arg.vec[1]

library(data.table)

coverage <- fread(file.path(problem.dir, "coverage.bedGraph"))
setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))
bases <- with(coverage, chromEnd-chromStart)
long <- rep(coverage$count, bases)
diff.vec <- abs(diff(long))
feature.vec <- c(
  quartile=quantile(long),
  mean=mean(long),
  sd=sd(long),
  bases=sum(bases),
  data=nrow(coverage))
log.features <-
  c(feature.vec,
    `log+1`=log(feature.vec+1),
    log=log(feature.vec),
    log.log=log(log(feature.vec)))
feature.dt <- data.table(t(log.features))
write.table(
  feature.dt,
  file.path(problem.dir, "features.tsv"),
  quote=FALSE,
  row.names=FALSE,
  col.names=TRUE,
  sep="\t")
