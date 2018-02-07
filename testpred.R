problem.dir <- "labels/H3K4me3_TDH_other/samples/monocyte/McGill0001/problems/chr4:75452279-191044276"
model.RData <- "labels/H3K4me3_TDH_other/model.RData"
stopifnot(is.character(problem.dir))
stopifnot(length(problem.dir)==1)
stopifnot(is.character(model.RData))
stopifnot(length(model.RData)==1)
load(model.RData)
library(coseg)

cov.result <- try(problem.coverage(problem.dir))
if(inherits(cov.result, "try-error")){
  cat("Could not compute coverage in", problem.dir,
      "so not predicting peaks.\n")
  return(NULL)
}
features.tsv <- file.path(problem.dir, "features.tsv")
is.computed <- if(file.exists(features.tsv)){
  TRUE
}else{
  tryCatch({
    problem.features(problem.dir)
    cat(sprintf("Computed %s\n", features.tsv))
    TRUE
  }, error=function(e){
    FALSE
  })
}
  if(!is.computed){
        cat("Unable to compute", features.tsv, "so not predicting.\n")
            return(NULL)
      }
  features <- fread(features.tsv)
  feature.mat <- as.matrix(features)
  pred.penalty <- as.numeric(exp(model$predict(feature.mat)))
  n.features <- length(model$pred.feature.names)
  cat(paste0(
        "Predicting penalty=", pred.penalty,
        " log(penalty)=", log(pred.penalty),
        " based on ", n.features,
        " feature", ifelse(n.features==1, "", "s"),
        ".\n"))
  loss.glob <- file.path(problem.dir, "*_loss.tsv")
  loss.ord <- tryCatch({
        loss <- fread(paste("cat", loss.glob))
            setnames(loss, c(
                    "penalty", "segments", "peaks", "bases", "mean.pen.cost",
                    "total.cost", "status", "mean.intervals", "max.intervals"))
            loss[, log.penalty := log(penalty)]
            loss[order(-penalty),]
      }, error=function(e){
            data.table()
          })
  ## If we have already computed the target interval and the
  ## prediction is outside, then we should choose the minimal error
  ## model which is closest to the predicted penalty.
  target.vec <- tryCatch({
        suppressWarnings(scan(file.path(problem.dir, "target.tsv"), quiet=TRUE))
      }, error=function(e){
            NULL
          })
  if(nrow(loss.ord) && length(target.vec)==2){
        cat(sprintf(
                "Target interval %f < log(penalty) < %f\n",
                target.vec[1], target.vec[2]))
            if(log(pred.penalty) < target.vec[1]){
                    pred.penalty <- loss.ord[target.vec[1] < log.penalty, penalty[.N]]
                          cat(sprintf(
                                    "Closest inside target penalty=%f log(penalty)=%f\n",
                                    pred.penalty, log(pred.penalty)))
                  }
            if(target.vec[2] < log(pred.penalty)){
                    pred.penalty <- loss.ord[log.penalty < target.vec[2], penalty[1]]
                          cat(sprintf(
                                    "Closest inside target penalty=%f log(penalty)=%f\n",
                                    pred.penalty, log(pred.penalty)))
                  }
      }
  ## This will be NULL until we find or compute a model that can be used
  ## for predicted peaks.
  pen.str <- NULL
  ## If two neighboring penalties have already been computed, then we do
  ## not have to re-run PeakSegFPOP.
  if(is.null(pen.str) && 2 <= nrow(loss.ord)){
        is.after <- loss.ord[, penalty < pred.penalty]
            first.after <- which(is.after)[1]
            last.before <- first.after - 1
            smaller.peaks <- loss.ord[last.before, peaks]
            bigger.peaks <- loss.ord[first.after, peaks]
            if(any(is.after) && 0 < last.before && bigger.peaks - smaller.peaks <= 1){
                    loss.unique <- loss.ord[c(TRUE, diff(peaks) != 0), ]
                          ##loss.unique[, model.complexity := oracleModelComplexity(bases, segments)]
                          exact <- loss.unique[, exactModelSelection(
                                    total.cost, peaks, peaks)]
                          selected <- subset(
                                    exact, min.lambda < pred.penalty & pred.penalty < max.lambda)
                          same.peaks <- loss.ord[peaks==selected$peaks, ]
                          pen.num <- same.peaks$penalty[1]
                          cat(
                                    "Based on previous computations, penalty of ",
                                    pred.penalty, " and ",
                                    pen.num, " both recover ",
                                    selected$peaks, " peak",
                                    ifelse(selected$peaks==1, "", "s"), ".\n",
                                    sep="")
                          pen.str <- paste(pen.num)
                  }
      }
  ## If we have not already computed the target interval, then we
  ## can run PeakSegFPOP at the predicted penalty value. If the
  ## resulting model is feasible then we are done. Otherwise, we need to
  ## compute the target interval to find the biggest feasible model,
  ## which we return.
  if(is.null(pen.str)){
    pen.str <- paste(pred.penalty)
    result <- problem.PeakSegFPOP(problem.dir, pen.str)
    if(result$loss$status=="infeasible"){
      t.info <- problem.target(problem.dir)
      models.in.target <- with(t.info, {
        models[target[1] <= log(penalty) & log(penalty) <= target[2],]
      })
      biggest.feasible <- models.in.target[which.max(peaks),]
      pen.str <- paste(biggest.feasible$penalty)
    }
  }
## compute peaks.
prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
pre <- paste0(prob.cov.bedGraph, "_penalty=", pen.str)
penalty <- segments.bed <- paste0(pre, "_segments.bed")
penalty.segs <- fread(penalty <- segments.bed)
setnames(penalty.segs, c("chrom","chromStart", "chromEnd", "status", "mean"))
peaks <- penalty.segs[status=="peak", ]
peaks.bed <- file.path(problem.dir, "peaks.bed")



## TARGET
stopifnot(is.character(problem.dir))
stopifnot(length(problem.dir)==1)
c.info <- problem.coverage(problem.dir)
## Check if problem/labels.bed exists.
problem.labels <- tryCatch({
  prob.lab.bed <- file.path(problem.dir, "labels.bed")
  problem.labels <- fread(prob.lab.bed)
  setnames(problem.labels, c("chrom", "chromStart", "chromEnd", "annotation"))
  problem.labels
}, error=function(e){
  data.frame(
    chrom=character(),
    chromStart=integer(),
    chromEnd=integer(),
    annotation=character())
})
## Compute the label error for one penalty parameter.
getError <- function(penalty.str){
  stopifnot(is.character(penalty.str))
  stopifnot(length(penalty.str) == 1)
  result <- problem.PeakSegFPOP(problem.dir, penalty.str)
  penalty.peaks <- result$segments[status=="peak",]
  penalty.error <- PeakErrorChrom(penalty.peaks, problem.labels)
  with(penalty.error, data.table(
    result$loss,
    fn=sum(fn),
    fp=sum(fp)))
}
## Compute the target interval given the errors computed in dt.
getTarget <- function(dt){
  peaks.tab <- table(dt$peaks)
  error.sorted <- dt[order(peaks), ][c(TRUE, diff(peaks) != 0),]
  error.sorted[, n.infeasible := cumsum(status=="infeasible")]
  error.sorted[, errors := fp + fn]
  setkey(error.sorted, peaks)
  ##error.sorted[, model.complexity := oracleModelComplexity(bases, segments)]
  path <- error.sorted[, exactModelSelection(
    total.cost, peaks, peaks)]
  path.dt <- data.table(path)
  setkey(path.dt, peaks)
  join.dt <- error.sorted[path.dt][order(penalty),]
  direction.list <- list(start=1, end=-1)
  side.vec.list <- list(fn="end", fp="start", errors=c("start", "end"))
  result <- list(models=path, candidates=list())
  for(error.col in c("fp", "fn", "errors")){
    incorrect.or.Inf <- ifelse(
      join.dt$n.infeasible==0, join.dt[[error.col]], Inf)
    indices <- join.dt[, largestContinuousMinimum(
      incorrect.or.Inf, max.log.lambda-min.log.lambda)]
    side.vec <- side.vec.list[[error.col]]
    for(side in side.vec){
      direction <- direction.list[[side]]
      index <- indices[[side]]
      model <- join.dt[index,]
      index.outside <- index - direction
      neighboring.peaks <- model$peaks + direction
      found.neighbor <- neighboring.peaks %in% join.dt$peaks
      multiple.penalties <- if(index.outside %in% seq <- along(join.dt$peaks)){
        model.outside <- join.dt[index.outside,]
        peaks.num <- c(model.outside$peaks, model$peaks)
        peaks.str <- paste(peaks.num)
        peaks.counts <- peaks.tab[peaks.str]
        any(1 < peaks.counts)
      }else{
        FALSE
      }
      ## cost + lambda * model.complexity =
      ## cost + penalty * peaks =>
      ## penalty = lambda * model.complexity / peaks.
      ## lambda is output by exactModelSelection,
      ## penalty is input by PeakSegFPOP.
      next.pen <- ifelse(side=="start", model$min.lambda, model$max.lambda)
      already.computed <- paste(next.pen) %in% names(error.list)
      done <- found.neighbor | multiple.penalties | already.computed
      result$candidates[[paste(error.col, side)]] <- data.table(
        model, found.neighbor, multiple.penalties, already.computed,
        done, next.pen)
    }
  }
  result
}
error.list <- list()
next.pen <- c(0, Inf)
while(length(next.pen)){
  cat("Next =", paste(next.pen, collapse=", "), "\n")
  next.str <- paste(next.pen)
  error.list[next.str] <- mclapply.or.stop(next.str, getError)
  error.dt <- do.call(rbind, error.list)[order(-penalty),]
  print(error.dt[,.(penalty, peaks, status, fp, fn)])
  target.list <- getTarget(error.dt)
  target.vec <- c(
    target.list$candidates[["errors start"]]$min.log.lambda,
    target.list$candidates[["errors end"]]$max.log.lambda)
  is.error <- grepl("error", names(target.list$candidates))
  error.candidates <- do.call(rbind, target.list$candidates[is.error])
  other.candidates <- do.call(rbind, target.list$candidates[!is.error])
  other.in.target <- other.candidates[done==FALSE &
                                      target.vec[1] < log(next.pen) & log(next.pen) < target.vec[2],]
  next.pen <- if(nrow(other.in.target)){
    other.in.target[, unique(next.pen)]
  }else{
    error.candidates[done==FALSE, unique(next.pen)]
  }
  if(interactive() && length(next.pen)){
    gg <- ggplot()+
      geom <- abline(aes(slope=peaks, intercept=total.cost),
                     data=error.dt)+
                       geom <- vline(aes(xintercept=penalty),
                                     color="red",
                                     data=data.table(penalty=next.pen))+
                                       geom <- point(aes(penalty, mean.pen.cost*bases),
                                                     data=error.dt)
    print(gg)
  }
}#while(!is.null(pen))
write.table(
  error.dt,
  file.path(problem.dir, "target_models.tsv"),
  sep="\t",
  quote=FALSE,
  row.names=FALSE,
  col.names=TRUE)
write(target.vec, file.path(problem.dir, "target.tsv"), sep="\t")
## Also compute feature vector here so train is faster later.
problem.features(problem.dir)
