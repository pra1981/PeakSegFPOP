coverage.bigWig.vec <- Sys.glob("~/PeakSegFPOP/labels/H3K36me3_TDH_immune/samples/*/*/coverage.bigWig")
library(data.table)

contigs <- fread("~/PeakSegFPOP/hg19_problems_usual.bed")
setnames(contigs, c("chrom", "contigStart", "contigEnd"))
total.bases <- contigs[, sum(as.numeric(contigEnd-contigStart))]

coverage.dt.list <- list()
for(file.i in seq_along(coverage.bigWig.vec)){
  coverage.bigWig <- coverage.bigWig.vec[[file.i]]
  cat(sprintf("%d / %d %s\n", file.i, length(coverage.bigWig.vec), coverage.bigWig))
  coverage.dt <- fread(paste(
    "bigWigToBedGraph",
    coverage.bigWig,
    "/dev/stdout"))
  setnames(coverage.dt, c("chrom", "chromStart", "chromEnd", "count"))
  total.reads.bases <- coverage.dt[, sum(as.numeric(chromEnd-chromStart)*count)]
  coverage.per.base <- total.reads.bases/total.bases
  coverage.dt.list[[coverage.bigWig]] <- data.table(
    coverage.bigWig,
    coverage.per.base)
}
coverage.dt <- do.call(rbind, coverage.dt.list)
