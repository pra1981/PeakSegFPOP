## inputs: problem directory with problem.bed.

## outputs: if the problem directory does not have coverage.bedGraph,
## it will be created via intersectBed. If labels.bed is present, we
## also create target.bed.
arg.vec <- "labels/H3K36me3_AM_immune_folds2-4/McGill0002/problems/chr1:3995268-13052998"
arg.vec <- "labels/H3K36me3_TDH_immune/McGill0001/problems/chr11:96437584-134946516"
arg.vec <- "labels/small/McGill0106/problems/chr1:17175658-29878082/"
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 1){
  stop("usage: Rscript compute_coverage_target.R data_dir/sample_dir/problems/problem_dir")
}
problem.dir <- arg.vec[1]
problems.dir <- dirname(problem.dir)
sample.dir <- dirname(problems.dir)

library(coseg)
library(data.table)
library(PeakError)

problem.bed <- file.path(problem.dir, "problem.bed")
problem <- fread(problem.bed)
setnames(problem, c("chrom", "problemStart", "problemEnd"))

## First check if problem/coverage.bedGraph has been created.
prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
coverage.ok <- tryCatch({
  tail.cmd <- paste("tail -1", prob.cov.bedGraph)
  last.cov <- fread(tail.cmd)
  setnames(last.cov, c("chrom", "chromStart", "chromEnd", "coverage"))
  last.cov$chromEnd == problem$problemEnd
}, error=function(e){
  FALSE
})

## Create problem/coverage.bedGraph if it is not present.
if(!coverage.ok){
  ## Use intersectBed -sorted to avoid memory
  ## problems. coverage.bedGraph needs to be -a since that is
  ## reported in the output.
  coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
  intersect.cmd <- paste(
    "intersectBed -sorted",
    "-a", coverage.bedGraph,
    "-b", problem.bed,
    ">", prob.cov.bedGraph)
  cat(intersect.cmd, "\n")
  system(intersect.cmd)
}

## Check if problem/labels.bed exists.
is.labeled <- tryCatch({
  prob.lab.bed <- file.path(problem.dir, "labels.bed")
  problem.labels <- fread(prob.lab.bed)
  setnames(problem.labels, c("chrom", "chromStart", "chromEnd", "annotation"))
  0 < nrow(problem.labels)
}, error=function(e){
  FALSE
})

## If this problem is labeled, then compute the target interval.
if(is.labeled){
  error.list <- list()
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
    penalty_segments.bed <- paste0(pre, "_segments.bed")
    if(!file.exists(penalty_segments.bed)){
      penalty.db <- paste0(pre, ".db")
      fpop.cmd <- paste(
        "PeakSegFPOP", prob.cov.bedGraph, penalty.str, penalty.db)
      cat(fpop.cmd, "\n")
      seconds <- system.time({
        system(fpop.cmd)
      })[["elapsed"]]
      megabytes <- if(file.exists(penalty.db)){
        file.size(penalty.db)/1024/1024
      }else{
        0
      }
      penalty_timing.tsv <- paste0(pre, "_timing.tsv")
      timing <- data.table(
        penalty.str,
        megabytes,
        seconds)
      write.table(
        timing,
        penalty_timing.tsv,
        row.names=FALSE, col.names=FALSE,
        quote=FALSE, sep="\t")
      unlink(penalty.db)
    }
    penalty_loss.tsv <- paste0(pre, "_loss.tsv")
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "status"))
    penalty.segs <- fread(penalty_segments.bed, colClasses=list(NULL=5))
    setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status"))
    penalty.peaks <- penalty.segs[status=="peak",]
    penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    error.list[[penalty.str]] <<- with(penalty.error, data.table(
      penalty.loss,
      fn=sum(fn),
      fp=sum(fp)))
    do.call(rbind, error.list)[order(penalty),]
  }

  error.dt <- getError("0")
  error.dt <- getError("Inf")
  min.fp <- error.list[["Inf"]]$fp
  min.fn <- error.list[["0"]]$fn

  ## mx+b = lossInf => x = (lossInf-b)/m
  lossInf <- error.list[["Inf"]]$total.cost
  next.pen <- with(error.list[["0"]], (lossInf-total.cost)/peaks)
  while(!is.null(next.pen)){
    if(interactive()){
      gg <- ggplot()+
        geom_abline(aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        geom_vline(aes(xintercept=penalty),
                   data=data.table(penalty=next.pen))+
        geom_point(aes(penalty, mean.pen.cost*bases),
                   data=error.dt)
      print(gg)
    }
    next.str <- paste(next.pen)
    error.dt <- getError(next.str)
    print(error.dt[,.(penalty, peaks, fp, fn)])
    peaks.tab <- table(error.dt$peaks)
    ## one sufficient condition for having found the lower limit is
    ## having found one (p,p+1) pair with fp values (0,>0)
    fp.above.min <- error.dt[min.fp < fp,]
    fp.is.min <- error.dt[fp==min.fp,]
    last.above <- fp.above.min[.N,]
    first.min <- fp.is.min[1,]
    fp.above.is.next <- first.min$peaks == last.above$peaks-1
    fp.two.lambda <-
      any(1 < peaks.tab[paste(c(last.above$peaks, first.min$peaks))])
    fp.found <- fp.above.is.next || fp.two.lambda
    ## one sufficient condition for having found the upper limit is
    ## having found one (p,p+1) pair with fn values (>min.fn,min.fn)
    fn.above.min <- error.dt[min.fn < fn, ]
    fn.is.min <- error.dt[fn==min.fn, ]
    last.min <- fn.is.min[.N, ]
    first.above <- fn.above.min[1, ]
    last.is.next <- last.min$peaks == first.above$peaks+1
    fn.two.lambda <-
      any(1 < peaks.tab[paste(c(first.above$peaks, last.min$peaks))])
    fn.found <- last.is.next || fn.two.lambda
    ## Rather than searching for when fn becomes minimum, search upper
    ## penalty limit of the min error.
    error.dt[, errors := fn+fp]
    min.error <- error.dt[, min(errors)]
    error.is.min <- error.dt[errors==min.error, ]
    bigger.pen <- error.dt[max(error.is.min$penalty) < penalty,]
    last.min.err <- error.is.min[.N, ]
    first.bigger <- bigger.pen[1, ]
    min.is.next <- last.min.err$peaks == first.bigger$peaks+1
    err.two.lambda <-
      any(1 < peaks.tab[paste(c(first.bigger$peaks, last.min.err$peaks))])
    err.found <- min.is.next || err.two.lambda
    next.pen <- if(!fp.found){
      ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
      (last.above$total.cost-first.min$total.cost)/
        (first.min$peaks-last.above$peaks)
    }else if(!err.found){
      ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
      (first.bigger$total.cost-last.min.err$total.cost)/
        (last.min.err$peaks-first.bigger$peaks)
    }else{
      NULL
    }
  }#while(!is.null(pen))

  error.sorted <- error.dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]
  path <- error.sorted[, exactModelSelection(total.cost, peaks, peaks)]
  setkey(error.sorted, peaks)
  path$errors <- error.sorted[J(path$peaks), fp+fn]
  indices <- with(path, largestContinuousMinimum(
    errors, max.log.lambda-min.log.lambda))
  target <- with(path, data.table(
    min.log.lambda=min.log.lambda[indices$start],
    max.log.lambda=max.log.lambda[indices$end]))

  write.table(
    target,
    file.path(problem.dir, "target.tsv"),
    quote=FALSE,
    col.names=FALSE,
    row.names=FALSE)
}
