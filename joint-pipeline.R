problem.joint.predict.many <- function
### Compute all joint peak predictions for one separate problem.
(prob.dir
### project/problems/problemID
){
  joint.dir.vec <- Sys.glob(file.path(
    prob.dir, "jointProblems", "*"))
  problems.dir <- dirname(prob.dir)
  set.dir <- dirname(problems.dir)
  joint.model.RData <- file.path(set.dir, "joint.model.RData")
  mclapply.or.stop(joint.dir.vec, function(joint.dir){
    PeakSegJoint::problem.joint.predict(joint.model.RData, joint.dir)
  })
### Nothing.
}

problem.joint.targets.train <- function
### Compute all target intervals then learn a penalty function.
(data.dir
### project directory.
){
  labels.tsv.vec <- Sys.glob(file.path(
    data.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
  mclapply.or.stop(labels.tsv.vec, function(labels.tsv){
    problem.joint.target(dirname(labels.tsv))
  })
  problem.joint.train(data.dir)
### Nothing.
}

problem.joint.train <- function
### Learn a penalty function for joint peak prediction.
(data.dir
### project directory.
){
  joint.model.RData <- file.path(data.dir, "joint.model.RData")
  target.tsv.vec <- Sys.glob(file.path(
    data.dir, "problems", "*", "jointProblems", "*", "target.tsv"))
  cat("Found", length(target.tsv.vec), "target.tsv files for training.\n")
  problems.list <- list()
  for(target.tsv.i in seq_along(target.tsv.vec)){
    target.tsv <- target.tsv.vec[[target.tsv.i]]
    target.vec <- scan(target.tsv, quiet=TRUE)
    problem.dir <- dirname(target.tsv)
    segmentations.RData <- file.path(problem.dir, "segmentations.RData")
    load(segmentations.RData)
    if(any(is.finite(target.vec))){
      problems.list[[problem.dir]] <- list(
        features=segmentations$features,
        target=target.vec)
    }
  }
  cat("Training using", length(problems.list), "finite targets.\n")
  set.seed(1)
  n.observations <- length(problems.list)
  n.folds <- ifelse(n.observations < 10, 3, 5)
  fold.vec <- sample(rep(1:n.folds, l=n.observations))
  squared.hinge <- function(x){
    ifelse(x<1,(x-1)^2,0)
  }
  error.loss.list <- list()
  for(validation.fold in unique(fold.vec)){
    ##cat(sprintf("%d"))
    is.validation <- fold.vec == validation.fold
    is.train <- !is.validation
    train.problems <- problems.list[is.train]
    fit <- IntervalRegressionProblems(
      train.problems, max.iterations=1e4, verbose=0)
    for(problem.i in seq_along(problems.list)){
      problem <- problems.list[[problem.i]]
      log.penalty.vec <- fit$predict(problem$features)
      too.lo <- log.penalty.vec < problem$target[1]
      too.hi <- problem$target[2] < log.penalty.vec
      left.term <- squared.hinge(log.penalty.vec-problem$target[1])
      right.term <- squared.hinge(problem$target[2]-log.penalty.vec)
      error.loss.list[[paste(validation.fold, problem.i)]] <- data.table(
        validation.fold,
        surrogate.loss=left.term + right.term,
        incorrect.targets=too.lo | too.hi,
        set.name=ifelse(is.train[problem.i], "train", "validation"),
        regularization=fit$regularization.vec,
        problem.i)
    }
  }
  error.loss <- do.call(rbind, error.loss.list)
  set.error.loss <- error.loss[, list(
    mean.surrogate.loss=mean(surrogate.loss),
    mean.incorrect.targets=mean(incorrect.targets)
  ), by=.(validation.fold, set.name, regularization)]
  validation.min <- set.error.loss[set.name=="validation", {
    .SD[mean.surrogate.loss==min(mean.surrogate.loss),]
  }, by=validation.fold]
  validation.min.simplest <- validation.min[, {
    .SD[regularization==max(regularization),]
  }, by=validation.fold]
  if(FALSE){
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(validation.fold ~ .)+
      geom_point(aes(-log(regularization), mean.surrogate.loss,
                     color=set.name),
                 data=validation.min.simplest)+
      geom_line(aes(-log(regularization), mean.surrogate.loss,
                    group=paste(validation.fold, set.name),
                    color=set.name),
                data=set.error.loss)
    ggplot()+
      theme_bw()+
      theme(panel.margin=grid::unit(0, "lines"))+
      facet_grid(validation.fold ~ .)+
      geom_line(aes(-log(regularization), mean.incorrect.targets,
                    group=paste(validation.fold, set.name),
                    color=set.name),
                data=set.error.loss)

  }
  mean.reg <- mean(validation.min.simplest$regularization)
  joint.model <- IntervalRegressionProblems(
    problems.list,
    initial.regularization=mean.reg,
    factor.regularization=NULL,
    verbose=0)
  cat("Learned regularization parameter and weights:\n")
  is.selected <- joint.model$param.mat != 0
  print(joint.model$param.mat[is.selected, , drop=FALSE])
  pred.log.penalty <- sapply(problems.list, function(problem){
    joint.model$predict(problem$features)
  })
  lower.limit <- sapply(problems.list, with, target[1])
  upper.limit <- sapply(problems.list, with, target[2])
  pred.dt <- data.table(
    too.lo=as.logical(pred.log.penalty < lower.limit),
    lower.limit,
    pred.log.penalty,
    upper.limit,
    too.hi=as.logical(upper.limit < pred.log.penalty))
  pred.dt[, status := ifelse(
    too.lo, "low",
    ifelse(too.hi, "high", "correct"))]
  cat("Train errors:\n")
  pred.dt[, list(targets=.N), by=status]
  save(joint.model, problems.list, file=joint.model.RData)
### Nothing.
}

problem.joint <- function
### Fit a joint model.
(jointProblem.dir
### path/to/jointProblem
){
  segmentations.RData <- file.path(jointProblem.dir, "segmentations.RData")
  problem.bed <- file.path(jointProblem.dir, "problem.bed")
  problem <- fread(problem.bed)
  setnames(problem, c("chrom",  "problemStart", "problemEnd", "problem.name"))
  problem[, problemStart1 := problemStart + 1L]
  setkey(problem, chrom, problemStart1, problemEnd)
  jointProblems <- dirname(jointProblem.dir)
  prob.dir <- dirname(jointProblems)
  probs.dir <- dirname(prob.dir)
  data.dir <- dirname(probs.dir)
  samples.dir <- file.path(data.dir, "samples")
  coverage.bedGraph.vec <- Sys.glob(file.path(
    samples.dir, "*", "*", "problems",
    problem$problem.name, "coverage.bedGraph"))
  cat("Found", length(coverage.bedGraph.vec), "samples to jointly segment.\n")
  coverage.list <- list()
  for(coverage.i in seq_along(coverage.bedGraph.vec)){
    coverage.bedGraph <- coverage.bedGraph.vec[[coverage.i]]
    ## cat(sprintf(
    ##   "%4d / %4d %s\n",
    ##   coverage.i, length(coverage.bedGraph.vec), coverage.bedGraph))
    sample.coverage <- fread(coverage.bedGraph)
    setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
    sample.coverage[, chromStart1 := chromStart + 1L]
    setkey(sample.coverage, chrom, chromStart1, chromEnd)
    problem.coverage <- foverlaps(sample.coverage, problem, nomatch=0L)
    problem.dir <- dirname(coverage.bedGraph)
    problems.dir <- dirname(problem.dir)
    sample.dir <- dirname(problems.dir)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    sample.path <- paste0(sample.group, "/", sample.id)
    coverage.list[[sample.path]] <- data.table(
      sample.id, sample.group, problem.coverage)
  }
  coverage <- do.call(rbind, coverage.list)
  setkey(coverage, sample.id, chrom, chromStart, chromEnd)
  profile.list <- ProfileList(coverage)
  fit <- PeakSegJointSeveral(coverage)
  segmentations <- ConvertModelList(fit)
  segmentations$features <- featureMatrix(profile.list)
  cat("Writing segmentation and features to", segmentations.RData, "\n")
  save(segmentations, file=segmentations.RData)
  segmentations$coverage <- coverage
  segmentations
### Model from ConvertModelList.
}

problem.joint.predict <- function
### Compute peak predictions for a joint problem.
(joint.model.RData,
### path/to/joint.model.RData file (which contains an object called joint.model)
 jointProblem.dir
### problem to predict peaks.
){
  load(joint.model.RData)
  converted <- problem.joint(jointProblem.dir)
  log.penalty <- joint.model$predict(converted$features)
  stopifnot(length(log.penalty)==1)
  selected <- subset(
    converted$modelSelection,
    min.log.lambda < log.penalty & log.penalty < max.log.lambda)
  loss.tsv <- file.path(jointProblem.dir, "loss.tsv")
  pred.dt <- if(selected$peaks == 0){
    unlink(loss.tsv)
    data.table()
  }else{
    selected.loss <- converted$loss[paste(selected$peaks), "loss"]
    flat.loss <- converted$loss["0", "loss"]
    loss.dt <- data.table(
      loss.diff=flat.loss-selected.loss)
    write.table(
      loss.dt,
      loss.tsv,
      quote=FALSE,
      sep="\t",
      col.names=FALSE,
      row.names=FALSE)
    pred.df <- subset(converted$peaks, peaks==selected$peaks)
    chrom <- paste(converted$coverage$chrom[1])
    with(pred.df, data.table(
      chrom,
      chromStart,
      chromEnd,
      name=paste0(sample.group, "/", sample.id),
      mean))
  }
  peaks.bed <- file.path(jointProblem.dir, "peaks.bed")
  cat("Writing ",
      nrow(pred.dt), " peaks to ",
      peaks.bed,
      "\n", sep="")
  write.table(
    pred.dt, peaks.bed,
    quote=FALSE,
    sep="\t",
    col.names=FALSE,
    row.names=FALSE)
### Nothing.
}

problem.joint.target <- function
### Compute target interval for a joint problem.
(jointProblem.dir
### Joint problem directory.
){
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
      facet_grid(sample.group + sample.id ~ ., scales="free")+
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
        size=1,
        fill=NA,
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
                data=converted$coverage,
                color="grey50")+
      geom_segment(aes(chromStart/1e3, 0,
                       xend=chromEnd/1e3, yend=0),
                   data=show.peaks.df,
                   color="deepskyblue",
                   size=2)
    ## geom_segment(aes(chromStart/1e3, mean,
    ##                  xend=chromEnd/1e3, yend=mean),
    ##              data=
    ##              color="green")
  }
  cat("Train error:\n")
  print(data.table(fit.error$modelSelection)[, .(
    min.log.lambda, max.log.lambda, peaks, errors)])
  target.tsv <- file.path(jointProblem.dir, "target.tsv")
  cat(
    "Writing target interval (",
    paste(fit.error$target, collapse=", "),
    ") to ", 
    target.tsv,
    "\n", sep="")
  write(fit.error$target, target.tsv, sep="\t")
### Nothing.
}

problem.joint.plot <- function
### Plot one chunk.
(chunk.dir
### project/problems/problemID/chunks/chunkID
){
  if(!require(ggplot2)){
    stop("please install ggplot2")
  }
  chunks.dir <- dirname(chunk.dir)
  prob.dir <- dirname(chunks.dir)
  prob.name <- basename(prob.dir)
  probs.dir <- dirname(prob.dir)
  proj.dir <- dirname(probs.dir)
  problem.dir.vec <- Sys.glob(file.path(
    proj.dir, "samples", "*", "*", "problems", prob.name))
  chunk <- fread(file.path(chunk.dir, "chunk.bed"))
  setnames(chunk, c("chrom", "chunkStart", "chunkEnd"))
  chunk[, chunkStart1 := chunkStart + 1L]
  setkey(chunk, chunkStart1, chunkEnd)
  labels <- fread(file.path(chunk.dir, "labels.tsv"))
  cat("Read",
      nrow(labels),
      "labels.\n")
  jointProblems <- fread(file.path(prob.dir, "jointProblems.bed"))
  setnames(jointProblems, c("chrom", "problemStart", "problemEnd"))
  jointProblems[, problemStart1 := problemStart + 1L]
  jointProblems[, problem.name := sprintf(
    "%s:%d-%d", chrom, problemStart, problemEnd)]
  setkey(jointProblems, problemStart1, problemEnd)
  probs.in.chunk <- foverlaps(jointProblems, chunk, nomatch=0L)
  probs.in.chunk$sample.group <- "problems"
  probs.in.chunk$sample.id <- "joint"
  cat("Read",
      nrow(jointProblems),
      "joint problems, plotting",
      nrow(probs.in.chunk),
      "in chunk.\n")
  coverage.list <- list()
  separate.peaks.list <- list()
  for(sample.i in seq_along(problem.dir.vec)){
    problem.dir <- problem.dir.vec[[sample.i]]
    problems.dir <- dirname(problem.dir)
    sample.dir <- dirname(problems.dir)
    sample.id <- basename(sample.dir)
    group.dir <- dirname(sample.dir)
    sample.group <- basename(group.dir)
    sample.coverage <- fread(file.path(problem.dir, "coverage.bedGraph"))
    setnames(sample.coverage, c("chrom", "chromStart", "chromEnd", "count"))
    sample.peaks <- fread(file.path(problem.dir, "peaks.bed"))
    setnames(sample.peaks, c("chrom", "peakStart", "peakEnd", "status", "mean"))
    sample.peaks[, peakStart1 := peakStart + 1L]
    sample.coverage[, chromStart1 := chromStart + 1L]
    setkey(sample.coverage, chromStart1, chromEnd)
    setkey(sample.peaks, peakStart1, peakEnd)
    chunk.cov <- foverlaps(sample.coverage, chunk, nomatch=0L)
    chunk.peaks <- foverlaps(sample.peaks, chunk, nomatch=0L)
    coverage.list[[problem.dir]] <- data.table(
      sample.id, sample.group, chunk.cov)
    if(nrow(chunk.peaks)){
      separate.peaks.list[[problem.dir]] <- data.table(
        sample.id, sample.group, chunk.peaks)
    }
  }
  coverage <- do.call(rbind, coverage.list)
  separate.peaks <- do.call(rbind, separate.peaks.list)
  separate.peaks$peak.type <- "separate"
  cat("Read",
      length(coverage.list),
      "samples of coverage.\n")
  cat("Read",
      length(separate.peaks.list),
      "samples of separate peak predictions.\n")
  joint.peaks.list <- list()
  for(joint.i in 1:nrow(probs.in.chunk)){
    prob <- probs.in.chunk[joint.i,]
    tryCatch({
      peaks <- fread(file.path(
        prob.dir, "jointProblems", prob$problem.name, "peaks.bed"))
      setnames(peaks, c("chrom", "peakStart", "peakEnd", "sample.path", "mean"))
      peaks[, sample.id := sub(".*/", "", sample.path)]
      peaks[, sample.group := sub("/.*", "", sample.path)]
      joint.peaks.list[[prob$problem.name]] <- peaks
    }, error=function(e){
      NULL
    })
  }
  cat("Read",
      length(joint.peaks.list),
      "joint peak predictions.\n")
  joint.peaks <- do.call(rbind, joint.peaks.list)
  ann.colors <-
    c(noPeaks="#f6f4bf",
      peakStart="#ffafaf",
      peakEnd="#ff4c4c",
      peaks="#a445ee")
  gg <- ggplot()+
    theme_bw()+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(sample.group + sample.id ~ ., scales="free")+
    scale_y_continuous(
      "aligned read coverage",
      breaks=function(limits){
        lim <- floor(limits[2])
        if(lim==0){
          Inf
        }else{
          lim
        }
      })+
    scale_x_continuous(paste(
      "position on",
      coverage$chrom[1],
      "(kb = kilo bases)"))+
    ## geom_tallrect(aes(
    ##   xmin=problemStart/1e3,
    ##   xmax=problemEnd/1e3),
    ##   alpha=0.5,
    ##   size=3,
    ##   color="black",
    ##   fill=NA,
    ##   data=probs.in.chunk)+
    geom_segment(aes(
      problemStart/1e3, 0,
      xend=problemEnd/1e3, yend=0),
      size=1,
      color="blue",
      data=probs.in.chunk)+
    geom_point(aes(
      problemStart/1e3, 0),
      color="blue",
      data=probs.in.chunk)+
    geom_tallrect(aes(
      xmin=chromStart/1e3, 
      xmax=chromEnd/1e3,
      fill=annotation), 
      alpha=0.5,
      data=labels)+
    scale_fill_manual("label", values=ann.colors)+
    scale_color_manual(values=c(separate="black", joint="deepskyblue"))+
    scale_size_manual(values=c(separate=2, joint=3))+
    geom_step(aes(
      chromStart/1e3, count),
      data=coverage,
      color="grey50")
  if(length(joint.peaks)){
    joint.peaks$peak.type <- "joint"
    gg <- gg+
      geom_point(aes(
        peakStart/1e3, 0,
        color=peak.type,
        size=peak.type),
        data=joint.peaks)+
      geom_segment(aes(
        peakStart/1e3, 0,
        xend=peakEnd/1e3, yend=0,
        color=peak.type,
        size=peak.type),
        data=joint.peaks)
  }
  gg <- gg+
    geom_segment(aes(
      peakStart/1e3, 0,
      xend=peakEnd/1e3, yend=0,
      color=peak.type,
      size=peak.type),
      data=separate.peaks)+
    geom_point(aes(
      peakStart/1e3, 0,
      color=peak.type,
      size=peak.type),
      data=separate.peaks)
  n.rows <- length(coverage.list) + 2
  mypng <- function(base, g){
    f <- file.path(chunk.dir, base)
    cat("Writing ",
        f,
        "\n", sep="")
    png(f, res=100, width=1000, height=60*n.rows)
    print(g)
    dev.off()
    thumb.png <- sub(".png$", "-thumb.png", f)
    cmd <- sprintf("convert %s -resize 230 %s", f, thumb.png)
    system(cmd)
  }
  mypng("figure-predictions-zoomout.png", gg)
  gg.zoom <- gg+
    coord_cartesian(
      xlim=chunk[, c(chunkStart, chunkEnd)/1e3],
      expand=FALSE)
  mypng("figure-predictions.png", gg.zoom)
### Nothing
}  
