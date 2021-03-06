library(data.table)
labels.bed.vec <- Sys.glob("labels/*/*/problems/*/labels.bed")

timing.data.list <- list()
problem.size.list <- list()
for(labels.bed in labels.bed.vec){
  problem.dir <- dirname(labels.bed)
  coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
  n.data <- wc(coverage.bedGraph)
  problem.size.list[[problem.dir]] <- data.table(
    problem.dir, n.data)
  tryCatch({
    timing.cmd <- paste("cat", file.path(problem.dir, "*timing.tsv"))
    timing.dt <- fread(timing.cmd)
    setnames(timing.dt, c("penalty", "megabytes", "seconds"))
    loss.cmd <- paste("cat", file.path(problem.dir, "*loss.tsv"))
    loss.dt <- fread(loss.cmd)
    setnames(loss.dt, c("penalty", "segments", "peaks", "bases", "mean.pen.cost", "total.cost", "status", "mean.intervals", "max.intervals"))
    setkey(loss.dt, penalty)
    setkey(timing.dt, penalty)
    problem.dt <- timing.dt[loss.dt]
    timing.data.list[[problem.dir]] <- data.table(
      problem.dir,
      n.data,
      problem.dt)
  }, error=function(e){
    ##do nothing.
  })
}
timing.data <- do.call(rbind, timing.data.list)
problem.size <- do.call(rbind, problem.size.list)

thresh <- 1e4
some.data <- timing.data[penalty < Inf & thresh < n.data,]
some.size <- problem.size[thresh < n.data, ]
ggplot()+
  geom_point(aes(log10(n.data), -Inf),
             data=some.size)+
  geom_point(aes(log10(n.data), log10(seconds)),
             data=some.data)

setkey(timing.data, problem.dir)
target.bed.vec <- Sys.glob("labels/*/*/problems/*/target.tsv")
problems.with.target <- dirname(target.bed.vec)
timing.data[problems.with.target][, list(
  total.minutes=sum(seconds/60)
  ), by=problem.dir]
