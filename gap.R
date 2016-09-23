library(coseg)
library(data.table)
library(PeakError)
data("hg19.gap", package="cosegData")
hg19.problems <- data.table(hg19.gap)[, data.table(
  problemStart=chromEnd[-.N],
  problemEnd=chromStart[-1]),
  by=chrom]
hg19.problems[, problem.name := paste0(chrom, ":", problemStart, "-", problemEnd)]
load("labels.RData")
sname <- "H3K36me3_TDH_immune"
one.set <- labels[sample.id=="McGill0001" & set.name==sname,]
setkey(one.set, chrom, chromStart, chromEnd)
setkey(hg19.problems, chrom, problemStart, problemEnd)
over.dt <- foverlaps(hg19.problems, one.set, nomatch=0L)
coverage.bedGraph <- "data/McGill0001/H3K36me3.bedGraph"

setkey(hg19.problems, problem.name)
labels.by.problem <- split(over.dt, over.dt$problem.name)
for(problem.i in seq_along(labels.by.problem)){
  pname <- names(labels.by.problem)[[problem.i]]
  problem <- hg19.problems[pname,]
  problem.dir <- file.path("labels", sname, pname)
  problem.labels <- labels.by.problem[[problem.i]]
  dir.create(problem.dir, showWarnings=FALSE, recursive=TRUE)
  problem.bed <- file.path(problem.dir, "problem.bed")
  write.table(
    problem, problem.bed,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
  prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  if(!file.exists(prob.cov.bedGraph)){
    ## Use intersectBed to avoid memory problems.
    cmd <- paste(
      "intersectBed -sorted",
      "-a", coverage.bedGraph,
      "-b", problem.bed,
      ">", prob.cov.bedGraph)
    system(cmd)
  }
  error.list <- list()
  getError <- function(penalty.str){
    stopifnot(is.character(penalty.str))
    pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
    penalty_segments.bed <- paste0(pre, "_segments.bed")
    if(!file.exists(penalty_segments.bed)){
      penalty.db <- paste0(pre, ".db")
      fpop.cmd <- paste(
        "./PeakSegFPOP", prob.cov.bedGraph, penalty.str, penalty.db)
      system(fpop.cmd)
      unlink(penalty.db)
    }
    penalty_loss.tsv <- paste0(pre, "_loss.tsv")
    penalty.loss <- fread(penalty_loss.tsv)
    setnames(penalty.loss, c(
      "penalty", "segments", "peaks", "bases",
      "mean.pen.cost", "total.cost", "status"))
    read.cmd <- paste("grep peak", penalty_segments.bed)
    penalty.peaks <- fread(read.cmd, colClasses=list(NULL=4:5))
    setnames(penalty.peaks, c("chrom","chromStart", "chromEnd"))
    penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
    error.list[[penalty.str]] <<- with(penalty.error, data.table(
      penalty.loss,
      fn=sum(fn),
      fp=sum(fp)))
  }
  getError("0")
  problem.coverage <- fread(prob.cov.bedGraph)
  setnames(problem.coverage, c("chrom", "chromStart", "chromEnd", "count"))
  problem.coverage[, weight := chromEnd-chromStart]
  problem.mean <- problem.coverage[, sum(weight*count)/sum(weight)]
  lossInf <- problem.coverage[, PoissonLoss(count, problem.mean, weight)]
  bases <- sum(problem.coverage$weight)
  errorInf <- PeakErrorChrom(Peaks(), problem.labels)
  error.list[["Inf"]] <- eInf <- data.table(
    penalty=Inf,
    segments=1,
    peaks=0,
    bases,
    mean.pen.cost=lossInf/bases,
    total.cost=lossInf,
    status="feasible",
    fn=sum(errorInf$fn),
    fp=sum(errorInf$fp))
  min.fp <- eInf$fp
  error.dt <- do.call(rbind, error.list)
  e0 <- error.list[["0"]]
  min.fn <- e0$fn
  done <- FALSE
  while(!done){
    error.dt <- do.call(rbind, error.list)[order(penalty),]
    peaks.tab <- table(error.dt$peaks)
    if(any(1 < peaks.tab)){
      print(peaks.tab)
      stop("found two lambda which yield the same number of peaks")
    }
    ## one sufficient condition for having found the lower limit is
    ## having found one (p,p+1) pair with fp values (0,>0)
    fp.above.min <- error.dt[min.fp < fp,]
    fp.is.min <- error.dt[fp==min.fp,]
    last.above <- fp.above.min[.N,]
    first.min <- fp.is.min[1,]
    fp.above.is.next <- first.min$peaks == last.above$peaks-1
    ## one sufficient condition for having found the upper limit is
    ## having found one (p,p+1) pair with fn values (>min.fn,min.fn)
    fn.above.min <- error.dt[min.fn < fn, ]
    fn.is.min <- error.dt[fn==min.fn, ]
    last.min <- fn.is.min[.N, ]
    first.above <- fn.above.min[1, ]
    last.is.next <- last.min$peaks == first.above$peaks+1
    if(!fp.above.is.next){
      ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
      next.pen <-
        (last.above$total.cost-first.min$total.cost)/
        (first.min$peaks-last.above$peaks)
      next.dt <- data.table(
        penalty=next.pen,
        cost=first.min[, total.cost + peaks*next.pen])
    }else if(!last.is.next){
      ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
      next.pen <-
        (first.above$total.cost-last.min$total.cost)/
        (last.min$peaks-first.above$peaks)
      next.dt <- data.table(
        penalty=next.pen,
        cost=first.above[, total.cost + peaks*next.pen])
    }else{
      done <- TRUE
    }
    if(!done){
      ggplot()+
        geom_abline(aes(slope=peaks, intercept=total.cost),
                    data=error.dt)+
        geom_vline(aes(xintercept=penalty),
                   data=next.dt)+
        geom_point(aes(penalty, mean.pen.cost*bases),
                   data=error.dt)
      next.str <- paste(next.pen)
      getError(next.str)
    }
  }#while(!done)
  error.dt[,.(penalty, peaks, fp, fn)]
  error.sorted <- error.dt[order(peaks),]
  path <- error.sorted[, exactModelSelection(total.cost, peaks, peaks)]
  setkey(error.sorted, peaks)
  path$errors <- error.sorted[J(path$peaks), fp+fn]
  indices <- with(path, largestContinuousMinimum(
    errors, max.log.lambda-min.log.lambda))
  target <- with(path, data.table(
    min.log.lambda=min.log.lambda[indices$start],
    max.log.lambda=max.log.lambda[indices$end]))
}

