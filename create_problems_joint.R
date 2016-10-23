## Edit the following definition to reflect your cluster
## configuration.
PBS.header <- "#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V"

arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "chr1:17175658-29878082")
arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "24")
arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "22")
arg.vec <- c(
  "test/H3K4me3_TDH_other/samples",
  "7")
arg.vec <- c(
  "/home/tdhock/PeakSegFPOP/test/input/samples",
  "chr10:18024675-38818835")
arg.vec <- c(
  "test/input/samples",
  "chr10:38868835-39154935")

arg.vec <- commandArgs(trailingOnly=TRUE)

samples.dir <- normalizePath(arg.vec[1], mustWork=TRUE)
problem.name <- arg.vec[2]

ann.colors <-
  c(noPeaks="#f6f4bf",
    peakStart="#ffafaf",
    peakEnd="#ff4c4c",
    peaks="#a445ee")
library(data.table)
library(PeakSegJoint)#for clusterPeaks

peaks.glob <- file.path(
  samples.dir, "*", "*", "problems", problem.name, "peaks.bed")
peaks.bed.vec <- Sys.glob(peaks.glob)
cat("Found", length(peaks.bed.vec), peaks.glob, "files.\n")
peaks.list <- list()
for(sample.i in seq_along(peaks.bed.vec)){
  peaks.bed <- peaks.bed.vec[[sample.i]]
  problem.dir <- dirname(peaks.bed)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  sample.group <- basename(group.dir)
  peaks.list[[peaks.bed]] <- tryCatch({
    sample.peaks <- fread(peaks.bed)
    setnames(
      sample.peaks,
      c("chrom", "chromStart", "chromEnd", "status", "mean"))
    data.table(sample.id, sample.group, sample.peaks)
  }, error=function(e){
    ## do nothing
  })
}
peaks <- do.call(rbind, peaks.list)
problems.list <- if(is.data.frame(peaks) && 0 < nrow(peaks)){
  clustered <- clusterPeaks(peaks)
  clusters <- data.table(clustered)[, list(
    clusterStart=min(chromStart),
    clusterEnd=max(chromEnd)
  ), by=cluster]
  clusters[, clusterStart1 := clusterStart + 1L]
  setkey(clusters, clusterStart1, clusterEnd)
  cat(nrow(peaks), "total peaks form",
      nrow(clusters), "overlapping peak clusters.\n")
  list(
    peaks=clusters[, .(clusterStart, clusterEnd)])
}else{
  cat("No predicted peaks.\n")
  list()
}

labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", problem.name, "labels.bed"))
labels.list <- list()
for(sample.i in seq_along(labels.bed.vec)){
  labels.bed <- labels.bed.vec[[sample.i]]
  problem.dir <- dirname(labels.bed)
  problems.dir <- dirname(problem.dir)
  sample.dir <- dirname(problems.dir)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  sample.group <- basename(group.dir)
  sample.labels <- fread(labels.bed)
  setnames(
    sample.labels,
    c("chrom", "labelStart", "labelEnd", "annotation"))
  labels.list[[labels.bed]] <- 
    data.table(sample.id, sample.group, sample.labels)
}
labels <- do.call(rbind, labels.list)
if(is.null(labels)){
  cat("No labels.\n")
}else{
  label.props <- labels[, list(
    prop.noPeaks=mean(annotation=="noPeaks")
    ), by=.(labelStart, labelEnd)]
  label.props[, labelStart1 := labelStart + 1L]
  setkey(label.props, labelStart1, labelEnd)
  over <- foverlaps(label.props, clusters, nomatch=NA)
  labels.with.no.peaks <- over[is.na(cluster),]
  labels.with.no.peaks[, bases := labelEnd - labelStart]
  labels.with.no.peaks[, reduce := as.integer(bases/3)]
  cat(
    "Found", nrow(labels.with.no.peaks),
    "labeled regions with no peaks out of",
    nrow(label.props), "total.\n")
  problems.list$labels <- labels.with.no.peaks[, data.table(
    clusterStart=labelStart+reduce,
    clusterEnd=labelEnd-reduce)]
}

problems <- do.call(rbind, problems.list)
if(is.data.table(problems) && 0 < nrow(problems)){
  setkey(problems, clusterStart, clusterEnd)
  problems[, bases := clusterEnd - clusterStart]
  mid.between.problems <- problems[, as.integer((clusterEnd[-.N]+clusterStart[-1])/2)]
  problems[, mid.before := c(NA_integer_, mid.between.problems)]
  problems[, mid.after := c(mid.between.problems, NA_integer_)]
  problems[, problemStart := as.integer(clusterStart-bases)]
  problems[, problemEnd := as.integer(clusterEnd+bases)]
  problems[problemStart < mid.before, problemStart := mid.before]
  problems[mid.after < problemEnd, problemEnd := mid.after]
  chrom <- peaks$chrom[1]
  problem.info <- problems[, data.table(
    problemStart,
    problemEnd,
    problem.name=sprintf("%s:%d-%d", chrom, problemStart, problemEnd))]
  problem.info[, problemStart1 := problemStart + 1L]
  setkey(problem.info, problemStart1, problemEnd)

  if(!is.null(labels)){
    labels[, labelStart1 := labelStart + 1L]
    setkey(labels, labelStart1, labelEnd)
    problems.with.labels <- foverlaps(problem.info, labels, nomatch=0L)
    setkey(problems.with.labels, problem.name)
  }

  coverage.bedGraph.vec <- Sys.glob(file.path(
    samples.dir, "*", "*", "problems", problem.name, "coverage.bedGraph"))
  if(FALSE){

    ## TODO: click on the problem diagram to show coverage and models
    ## in one problem?
    
    problem.diagram <- ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.group ~ ., scales="free", space="free")+
      geom_tallrect(aes(
        xmin=clusterStart/1e3, 
        xmax=clusterEnd/1e3), 
        alpha=0.5,
        data=problems)+
      scale_color_manual(values=ann.colors)+
      geom_segment(aes(
        labelStart/1e3, sample.id,
        xend=labelEnd/1e3, yend=sample.id,
        color=annotation),
        size=10,
        data=labels)+
      geom_segment(aes(
        chromStart/1e3, sample.id,
        xend=chromEnd/1e3, yend=sample.id),
        size=1,
        data=peaks)+
      geom_point(aes(
        chromStart/1e3, sample.id),
        data=peaks)

    coverage.list <- list()
    for(coverage.bedGraph in coverage.bedGraph.vec){
      problem.dir <- dirname(coverage.bedGraph)
      problems.dir <- dirname(problem.dir)
      sample.dir <- dirname(problems.dir)
      sample.id <- basename(sample.dir)
      group.dir <- dirname(sample.dir)
      sample.group <- basename(group.dir)
      sample.cov <- fread(coverage.bedGraph)
      setnames(
        sample.cov,
        c("chrom", "chromStart", "chromEnd", "count"))
      coverage.list[[coverage.bedGraph]] <- 
        data.table(sample.id, sample.group, sample.cov)
    }
    coverage <- do.call(rbind, coverage.list)

    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(sample.id ~ ., scales="free")+
      scale_fill_manual(values=ann.colors)+
      geom_tallrect(aes(
        xmin=labelStart/1e3,
        xmax=labelEnd/1e3,
        fill=annotation),
        color="grey",
        alpha=0.5,
        data=labels)+
      geom_vline(aes(xintercept=problemStart/1e3),
                 data=problems,
                 color="red")+
      geom_segment(aes(problemStart/1e3, 0,
                       xend=problemEnd/1e3, yend=0),
                   color="red",
                   size=3,
                   data=problems)+
      geom_step(aes(chromStart/1e3, count),
                data=coverage,
                color="grey50")+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=peaks,
                   color="deepskyblue",
                   size=2)
    
  }

  data.dir <- dirname(samples.dir)
  jointProblems <- file.path(
    data.dir, "problems", problem.name, "jointProblems")
  unlink(jointProblems, recursive=TRUE)
  joint.model.RData <- file.path(data.dir, "joint.model.RData")
  ## TODO make directories with sh files for each.
  makeProblem <- function(problem.i){
    problem <- problem.info[problem.i,]
    pname <- problem$problem.name
    ##cat(sprintf("%4d / %4d problems %s\n", problem.i, nrow(problem.info), pname))
    problem.dir <- file.path(jointProblems, pname)
    dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
    pout <- data.table(
      chrom,
      problem[, .(problemStart, problemEnd)],
      problem.name)
    write.table(
      pout,
      file.path(problem.dir, "problem.bed"),
      quote=FALSE,
      sep="\t",
      row.names=FALSE,
      col.names=FALSE)
    if(!is.null(labels) && pname %in% problems.with.labels$problem.name){
      problem.labels <- problems.with.labels[pname]
      write.table(
        problem.labels[, .(
          chrom, labelStart, labelEnd, annotation, sample.id, sample.group)],
        file.path(problem.dir, "labels.tsv"),
        quote=FALSE,
        sep="\t",
        row.names=FALSE,
        col.names=FALSE)
      ## Script for target.
      target.tsv <- file.path(problem.dir, "target.tsv")
      sh.file <- paste0(target.tsv, ".sh")
      script.txt <- paste0(PBS.header, "
#PBS -o ", target.tsv, ".out
#PBS -e ", target.tsv, ".err
#PBS -N JTarget", pname, "
", "Rscript ", normalizePath("compute_joint_target.R", mustWork=TRUE), " ",
problem.dir, " 
")
      writeLines(script.txt, sh.file)
    }
    ## Script for peaks.
    peaks.bed <- file.path(problem.dir, "peaks.bed")
    sh.file <- paste0(peaks.bed, ".sh")
    script.txt <- paste0(PBS.header, "
#PBS -o ", peaks.bed, ".out
#PBS -e ", peaks.bed, ".err
#PBS -N JPred", problem$problem.name, "
", "Rscript ", normalizePath("predict_problem_joint.R", mustWork=TRUE), " ",
joint.model.RData, " ", problem.dir, " 
")
    writeLines(script.txt, sh.file)
  }

  cat(
    "Creating ", nrow(problem.info),
    " joint segmentation problems for ", problem.name,
    "\n", sep="")
  library(parallel)
  nothing <- lapply(1:nrow(problem.info), makeProblem)

  write.table(
    problem.info[, .(chrom, problemStart, problemEnd)],
    paste0(jointProblems, ".bed"),
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
}else{
  cat("No joint problems.\n")
}
