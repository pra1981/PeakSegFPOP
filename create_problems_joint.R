arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "chr1:17175658-29878082")

arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "24")

samples.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
problem.name <- arg.vec[2]

library(data.table)
library(PeakSegJoint)

peaks.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "problems", problem.name, "peaks.bed"))

peaks.list <- list()
coverage.list <- list()
for(sample.i in seq_along(peaks.bed.vec)){
  problem.dir <- dirname(peaks.bed)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.id <- basename(sample.dir)
  peaks.bed <- peaks.bed.vec[[sample.i]]
  cat(sprintf("%4d / %4d %s\n", sample.i, length(peaks.bed.vec), problem.dir))
  peaks.list[[peaks.bed]] <- tryCatch({
    sample.peaks <- fread(peaks.bed)
    setnames(
      sample.peaks,
      c("chrom", "chromStart", "chromEnd", "status", "mean"))
    data.table(sample.id, sample.peaks)
  }, error=function(e){
    ## do nothing
  })
  coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  coverage.list[[peaks.bed]] <- tryCatch({
    sample.cov <- fread(coverage.bedGraph)
    setnames(
      sample.cov,
      c("chrom", "chromStart", "chromEnd", "count"))
    data.table(sample.id, sample.cov)
  }, error=function(e){
    ## do nothing
  })
}
peaks <- do.call(rbind, peaks.list)
coverage <- do.call(rbind, coverage.list)

clustered <- clusterPeaks(peaks)
clusters <- data.table(clustered)[, list(
  peakStart=min(chromStart),
  peakEnd=max(chromEnd)
  ), by=cluster]
clusters[, bases := peakEnd - peakStart]
mid.between.clusters <- clusters[, (peakEnd[-.N]+peakStart[-1])/2]
clusters[, mid.before := c(NA, mid.between.clusters)]
clusters[, mid.after := c(mid.between.clusters, NA)]
clusters[, problemStart := as.integer(peakStart-bases)]
clusters[, problemEnd := as.integer(peakEnd+bases)]
clusters[problemStart < mid.before, problemStart := mid.before]
clusters[mid.after < problemEnd, problemEnd := mid.after]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample.id ~ ., scales="free")+
  geom_segment(aes(problemStart/1e3, 0,
                   xend=problemEnd/1e3, yend=0),
               color="red",
               size=3,
               data=clusters)+
  geom_step(aes(chromStart/1e3, count),
            data=coverage,
            color="grey50")+
  geom_segment(aes(chromStart/1e3, 0,
                   xend=chromEnd/1e3, yend=0),
               data=peaks,
               color="deepskyblue",
               size=2)

jointProblems <- file.path(samples.dir, "jointProblems")
##TODO make directories with coverage and sh files for each.
