arg.vec <- "test/H3K4me3_TDH_other/jointProblems/chr4:88923952-88935469"

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 1){
  stop("usage: Rscript compute_coverage_target.R data_dir/sample_dir/problems/problem_dir")
}
jointProblem.dir <- normalizePath(arg.vec[1], mustWork=TRUE)

library(PeakSegJoint)
library(data.table)

converted <- problem.joint(jointProblem.dir)

labels.bed <- file.path(jointProblem.dir, "labels.tsv")
labels <- fread(labels.bed)
setnames(labels, c(
  "chrom", "chromStart", "chromEnd", "annotation",
  "sample.id", "sample.group"))

fit.error <- PeakSegJointError(converted, labels)

if(FALSE){
  show.peaks <- 8
  show.peaks.df <- subset(converted$peaks, peaks==show.peaks)
  show.errors <- fit.error$error.regions[[paste(show.peaks)]]
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
    geom_tallrect(aes(
      xmin=chromStart/1e3,
      xmax=chromEnd/1e3,
      linetype=status),
      color="black",
      data=show.errors)+
    scale_linetype_manual(
      "error type",
      limits=c("correct", 
               "false negative",
               "false positive"),
      values=c(correct=0,
               "false negative"=3,
               "false positive"=1))+
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
}

target.tsv <- file.path(jointProblem.dir, "target.tsv")
write(fit.error$target, target.tsv, sep="\t")
