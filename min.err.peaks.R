library(data.table)
library(namedCapture)
pattern <- paste0(
  "walltime=",
  "(?<hours>[0-9]+)",
  ":",
  "(?<minutes>[0-9]+)",
  ":",
  "(?<seconds>[0-9]+)")
min.err.peaks.list <- list()
for(interval.type in c("problems","problems_infeasibleInf")){
  target.models.vec <- Sys.glob(file.path(
    "labels/*_TDH_immune/samples/*/*",
    interval.type, "*/target_models.tsv"))
  for(target.models.i in seq_along(target.models.vec)){
    target.models.file <- target.models.vec[[target.models.i]]
    cat(sprintf(
      "%4d / %4d %s\n",
      target.models.i,
      length(target.models.vec),
      target.models.file))
    target.models <- fread(target.models.file)
    ##target.models[, errors := ifelse(status=="infeasible", Inf, fp+fn)] #old way
    target.models[, errors := fp+fn] 
    min.incorrect.labels <- target.models[, min(errors)]
    min.err.models <- target.models[errors == min.incorrect.labels, ]
    prob.dir <- dirname(target.models.file)
    out.file <- file.path(prob.dir, "target.tsv.out")
    cmd <- paste("grep Resources", out.file)
    Resources <- system(cmd, intern=TRUE)
    stopifnot(length(Resources)==1)
    match.row <- str_match_named(Resources, pattern, list(
      hours=as.integer, minutes=as.integer, seconds=as.integer))
    timing.glob <- file.path(prob.dir, "*timing.tsv")
    timing <- fread(paste("cat", timing.glob))
    setnames(timing, c("penalty", "megabytes", "seconds"))
    problems.dir <- dirname(prob.dir)
    sample.dir <- dirname(problems.dir)
    group.dir <- dirname(sample.dir)
    samples.dir <- dirname(group.dir)
    set.dir <- dirname(samples.dir)
    min.err.peaks.list[[target.models.file]] <- min.err.models[, data.table(
      interval.type,
      cpu.hours=sum(timing$seconds)/60/60,
      walltime.hours=with(match.row, hours + minutes/60 + seconds/60/60),
      min.incorrect.labels,
      penalties=nrow(target.models),
      min.peaks=peaks[1],
      max.peaks=peaks[.N],
      bases=bases[1],
      n.feasible=sum(status=="feasible"),
      n.infeasible=sum(status=="infeasible"),
      experiment=gsub("_.*", "", basename(set.dir)),
      problem=basename(prob.dir),
      sample=basename(sample.dir),
      group=basename(group.dir))]
  }
}
min.err.peaks <- do.call(rbind, min.err.peaks.list)

min.err.peaks

save(min.err.peaks, file="min.err.peaks.RData")
