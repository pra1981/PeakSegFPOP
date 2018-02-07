models.vec <- Sys.glob(file.path("labels/*/samples/*/*/problems/*/target_models.tsv"))

target.models.list <- list()
for(models.i in seq_along(models.vec)){
  models.tsv <- models.vec[[models.i]]
  models <- fread(models.tsv)
  problem.dir <- dirname(models.tsv)
  cat(sprintf("%4d / %4d %s\n", models.i, length(models.vec), problem.dir))
  timings <- fread(paste("cat", file.path(problem.dir, "*timing.tsv")))
  setnames(timings, c("penalty", "megabytes", "seconds"))
  setkey(timings, penalty)
  setkey(models, penalty)
  target.models.list[[problem.dir]] <- data.table(
    problem.dir, timings[models])
}

