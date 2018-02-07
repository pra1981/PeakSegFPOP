stopifnot(is.character(problem.dir))
stopifnot(length(problem.dir) == 1)
stopifnot(is.character(model.RData))
stopifnot(length(model.RData) == 1)
load(model.RData)
problem.coverage(problem.dir)
features.tsv <- file.path(problem.dir, "features.tsv")
is.computed <- if (file.exists(features.tsv)) {
  TRUE
}else {
  tryCatch({
    problem.features(problem.dir)
    cat(sprintf("Computed %s\n", features.tsv))
    TRUE
  }, error = function(e) {
    FALSE
  })
}
if (!is.computed) {
  cat("Unable to compute", features.tsv, "so not predicting.\n")
  return(NULL)
}
features <- fread(features.tsv)
feature.mat <- as.matrix(features)
pred.penalty <- as.numeric(exp(model$predict(feature.mat)))
n.features <- length(model$pred.feature.names)
cat(paste0("Predicting penalty=", pred.penalty, " log(penalty)=", 
           log(pred.penalty), " based on ", n.features, " feature", 
           ifelse(n.features == 1, "", "s"), ".\n"))
loss.glob <- file.path(problem.dir, "*_loss.tsv")
loss.file.vec <- Sys.glob(loss.glob)
loss.ord <- if (length(loss.file.vec)) {
  loss <- fread(paste("cat", loss.glob))
  setnames(loss, c("penalty", "segments", "peaks", "bases", 
                   "mean.pen.cost", "total.cost", "status", "mean.intervals", 
                   "max.intervals"))
  loss[, `:=`(log.penalty, log(penalty))]
  loss[order(-penalty), ]
}else {
  data.table()
}
target.vec <- tryCatch({
  suppressWarnings(scan(file.path(problem.dir, "target.tsv"), 
                        quiet = TRUE))
}, error = function(e) {
  NULL
})
if (nrow(loss.ord) && length(target.vec) == 2) {
  cat(sprintf("Target interval %f < log(penalty) < %f\n", 
              target.vec[1], target.vec[2]))
  if (log(pred.penalty) < target.vec[1]) {
    pred.penalty <- loss.ord[target.vec[1] < log.penalty, 
                             penalty[.N]]
    cat(sprintf("Closest inside target penalty=%f log(penalty)=%f\n", 
                pred.penalty, log(pred.penalty)))
  }
  if (target.vec[2] < log(pred.penalty)) {
    pred.penalty <- loss.ord[log.penalty < target.vec[2], 
                             penalty[1]]
    cat(sprintf("Closest inside target penalty=%f log(penalty)=%f\n", 
                pred.penalty, log(pred.penalty)))
  }
}
pen.str <- NULL
if (is.null(pen.str) && 2 <= nrow(loss.ord)) {
  is.after <- loss.ord[, penalty < pred.penalty]
  first.after <- which(is.after)[1]
  last.before <- first.after - 1
  smaller.peaks <- loss.ord[last.before, peaks]
  bigger.peaks <- loss.ord[first.after, peaks]
  if(any(is.after) && bigger.peaks - smaller.peaks <= 1) {
    loss.unique <- loss.ord[c(TRUE, diff(peaks) != 0), 
                            ]
    exact <- loss.unique[, exactModelSelection(total.cost, 
                                               peaks, peaks)]
    selected <- subset(exact, min.lambda < pred.penalty & 
                       pred.penalty < max.lambda)
    same.peaks <- loss.ord[peaks == selected$peaks, ]
    pen.num <- same.peaks$penalty[1]
    cat("Based on previous computations, penalty of ", 
        pred.penalty, " and ", pen.num, " both recover ", 
        selected$peaks, " peak", ifelse(selected$peaks == 
                                        1, "", "s"), ".\n", sep = "")
    pen.str <- paste(pen.num)
  }
}
if (is.null(pen.str)) {
  pen.str <- paste(pred.penalty)
  result <- problem.PeakSegFPOP(problem.dir, pen.str)
  if (result$loss$status == "infeasible") {
    t.info <- problem.target(problem.dir)
    models.in.target <- with(t.info, {
      models[target[1] <= log(penalty) & log(penalty) <= 
             target[2], ]
    })
    biggest.feasible <- models.in.target[which.max(peaks), 
                                         ]
    pen.str <- paste(biggest.feasible$penalty)
  }
}
prob.cov.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
pre <- paste0(prob.cov.bedGraph, "_penalty=", pen.str)
penalty_segments.bed <- paste0(pre, "_segments.bed")
penalty.segs <- fread(penalty_segments.bed)
setnames(penalty.segs, c("chrom", "chromStart", "chromEnd", 
                         "status", "mean"))
peaks <- penalty.segs[status == "peak", ]
peaks.bed <- file.path(problem.dir, "peaks.bed")
cat("Writing ", peaks.bed, " with ", nrow(peaks), " peak", 
    ifelse(nrow(peaks) == 1, "", "s"), " based on ", penalty_segments.bed, 
    ".\n", sep = "")
write.table(peaks, peaks.bed, quote = FALSE, sep = "\t", 
            col.names = FALSE, row.names = FALSE)
peaks
