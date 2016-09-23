## inputs: problem directory with coverage.bedGraph and labels.bed

## outputs: target.tsv
arg.vec <- "labels/H3K36me3_TDH_immune/McGill0001/problems/chr10:18024675-38818835"
arg.vec <- "labels/small/McGill0106/problems/chr1:17175658-29878082/"
arg.vec <- commandArgs(trailingOnly=TRUE)
if(length(arg.vec) != 1){
  stop("usage: Rscript computeTarget.R data_dir/sample_dir/problems/problem_dir")
}
problem.dir <- arg.vec[1]
library(coseg)
library(data.table)
library(PeakError)

prob.lab.bed <- file.path(problem.dir, "labels.bed")
prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
problem.labels <- fread(prob.lab.bed)
setnames(problem.labels, c("chrom", "chromStart", "chromEnd", "annotation"))

error.list <- list()
getError <- function(penalty.str){
  stopifnot(is.character(penalty.str))
  pre <- paste0(prob.cov.bedGraph, "_penalty=", penalty.str)
  penalty_segments.bed <- paste0(pre, "_segments.bed")
  if(!file.exists(penalty_segments.bed)){
    penalty.db <- paste0(pre, ".db")
    fpop.cmd <- paste(
      "./PeakSegFPOP", prob.cov.bedGraph, penalty.str, penalty.db)
    cat(fpop.cmd, "\n")
    system(fpop.cmd)
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
  next.pen <- if(!fp.found){
    ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
    (last.above$total.cost-first.min$total.cost)/
      (first.min$peaks-last.above$peaks)
  }else if(!fn.found){
    ## m2*x + b2 = m3*x + b3 => x = (b3-b2)/(m2-m3)
    (first.above$total.cost-last.min$total.cost)/
      (last.min$peaks-first.above$peaks)
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
  col.names=TRUE,
  row.names=FALSE)
