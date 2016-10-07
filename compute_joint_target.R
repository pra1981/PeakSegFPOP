arg.vec <- "test/H3K4me3_TDH_other/jointProblems/chr4:88923952-88935469"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_coverage_target.R data_dir/sample_dir/problems/problem_dir")
}
jointProblem.dir <- normalizePath(arg.vec[1], mustWork=TRUE)

library(PeakSegJoint)
library(data.table)

problem.bed <- file.path(jointProblem.dir, "problem.bed")
problem <- fread(problem.bed)
setnames(problem, c("chrom",  "problemStart", "problemEnd", "problem.name"))
problem[, problemStart1 := problemStart + 1L]
setkey(problem, chrom, problemStart1, problemEnd)

jointProblems <- dirname(jointProblem.dir)
data.dir <- dirname(jointProblems)
samples.dir <- file.path(data.dir, "samples")
coverage.bedGraph.vec <- Sys.glob(file.path(
  samples.dir, "*", "problems", problem$problem.name, "coverage.bedGraph"))

coverage.list <- list()
for(coverage.i in seq_along(coverage.bedGraph.vec)){
  coverage.bedGraph <- coverage.bedGraph.vec[[coverage.i]]
  cat(sprintf("%4d / %4d %s\n", coverage.i, length(coverage.bedGraph.vec), coverage.bedGraph))
  sample.coverage <- fread(coverage.bedGraph)
  setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
  sample.coverage[, chromStart1 := chromStart + 1L]
  setkey(sample.coverage, chrom, chromStart1, chromEnd)
  problem.coverage <- foverlaps(sample.coverage, problem, nomatch=0L)
  problem.dir <- dirname(coverage.bedGraph)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.id <- basename(sample.dir)
  coverage.list[[coverage.bedGraph]] <- data.table(
    sample.id, problem.coverage)
}
coverage <- do.call(rbind, coverage.list)
setkey(coverage, sample.id, chrom, chromStart, chromEnd)

labels.bed <- file.path(jointProblem.dir, "labels.bed")
labels <- fread(labels.bed)
setnames(labels, c("sample.id", "chrom", "chromStart", "chromEnd", "annotation"))

fit <- PeakSegJointSeveral(coverage)
converted <- ConvertModelList(fit)

show.peaks <- 8
show.peaks.df <- subset(converted$peaks, peaks==show.peaks)
errors <- PeakErrorSamples(show.peaks.df, labels)

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample.id ~ ., scales="free")+
  scale_fill_manual(values=ann.colors)+
  geom_tallrect(aes(
    xmin=chromStart/1e3,
    xmax=chromEnd/1e3,
    fill=annotation),
                color="grey",
                alpha=0.5,
                data=labels)+
  geom_step(aes(chromStart/1e3, count),
            data=coverage,
            color="grey50")+
  ## geom_segment(aes(chromStart/1e3, 0,
  ##                  xend=chromEnd/1e3, yend=0),
  ##              data=subset(converted$peaks, peaks==max(peaks)),
  ##              color="deepskyblue",
  ##              size=2)+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               data=subset(converted$segments, peaks==max(peaks)),
               color="green")



