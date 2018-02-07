library(data.table)
library(namedCapture)
code.pattern <- paste0(
  "(?<out>.*.out)",
  ".*?",
  "(?<code>-?[0-9]+)$")
glob.R <- "labels/*/samples/*/*/targets.R"
stopifnot(is.character(glob.R))
stopifnot(length(glob.R)==1)
R.vec <- Sys.glob(glob.R)
glob.out <- sub("R$", "out", glob.R)
stopifnot(glob.R != glob.out)
grep.cmd <- paste("grep code", glob.out)
grep.lines <- system(grep.cmd, intern=TRUE)
results <- data.table(code=NA_character_, R=R.vec)
setkey(results, R)
if(length(grep.lines)){
  code.mat <- str_match_named(grep.lines, code.pattern)
  results[sub("out$", "R", code.mat[, "out"]), code := code.mat[, "code"] ]
}
print(table(results$code, useNA="always"))
targets.R.vec <- results[code!=0 | is.na(code), R]
sample.dir.vec <- results[code!=0, dirname(R)]

for(targets.R in targets.R.vec){
  system(paste("qsub", targets.R))
}

for(sample.dir in sample.dir.vec){
  loss.dt <- fread(paste0("cat ", sample.dir, "/*_loss.tsv"))
}


sample.dir.vec <- normalizePath(Sys.glob("labels/*/samples/*/*"), mustWork=TRUE)
if(FALSE){
  bigwig.vec <- file.path(sample.dir.vec, "coverage.bigWig")
  file.rename(bigwig.vec, file.path(sample.dir.vec, "norm.bigWig"))
}

i.vec <- 1001:1851
i.vec <- 1:1000
for(sample.dir.i in i.vec){
  sample.dir <- normalizePath(sample.dir.vec[[sample.dir.i]], mustWork=TRUE)
  cat(sprintf("%4d / %4d %s\n", sample.dir.i, length(sample.dir.vec), sample.dir))
  script <- paste0('#!/home/thocking/bin/Rscript
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V
#PBS -o ', sample.dir, '/targets.out
#PBS -e ', sample.dir, '/targets.err
#PBS -N ', sub(".*samples/", "", sample.dir), '
sample.dir <- "', sample.dir, '"
coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
coverage.bigWig <- file.path(sample.dir, "coverage.bigWig")
norm.bigWig <- file.path(sample.dir, "norm.bigWig")
if(!file.exists(coverage.bigWig)){
  print(coverage.bigWig)
  if(file.exists(coverage.bedGraph)){
    cmd <- paste(
      "bedGraphToBigWig",
      coverage.bedGraph,
      "/home/thocking/PeakSegFPOP/hg19_chromInfo.txt",
      coverage.bigWig)
    system(cmd)
  }else if(file.exists(norm.bigWig)){
    PeakSegPipeline::denormalizeBigWig(norm.bigWig, coverage.bigWig)
  }else{
    stop("no coverage data")
  }
}
sample.models.list <- list()
sample.iterations.list <- list()
sample.features.list <- list()
target.tsv.vec <- Sys.glob(file.path(sample.dir, "problems", "*", "target.tsv"))
library(data.table)
library(PeakSegPipeline)
for(target.tsv in target.tsv.vec){
  prob.dir <- dirname(target.tsv)
  print(prob.dir)
  tryCatch({
    target.list <- problem.target(prob.dir)
  }, error=function(e){
    print("ERROR COMPUTING TARGET, DELETING LOSS FILES AND RECOMPUTING")
    unlink(file.path(prob.dir, "*_loss.tsv"))
    target.list <- problem.target(prob.dir)
  })
  coverage.bedGraph <- file.path(prob.dir, "coverage.bedGraph")
  bedGraph.lines <- wc(coverage.bedGraph)
  ##unlink(coverage.bedGraph)
  system(paste("gzip", coverage.bedGraph))
  timing.dt <- fread(paste0("cat ", prob.dir, "/*_timing.tsv"))
  setnames(timing.dt, c("penalty", "megabytes", "seconds"))
  setkey(target.list$models, penalty)
  setkey(timing.dt, penalty)
  models.dt <- timing.dt[target.list$models]
  features.dt <- fread(file.path(prob.dir, "features.tsv"))
  sample.features.list[[prob.dir]] <- data.table(prob.dir, features.dt)
  sample.models.list[[prob.dir]] <- data.table(prob.dir, bedGraph.lines, models.dt)
  sample.iterations.list[[prob.dir]] <- data.table(prob.dir, target.list$target.iterations)
}
sample.features <- do.call(rbind, sample.features.list)
sample.models <- do.call(rbind, sample.models.list)
sample.iterations <- do.call(rbind, sample.iterations.list)
targets.RData <- file.path(sample.dir, "targets.RData")
save(sample.models, sample.iterations, sample.features, file=targets.RData)
')
  targets.R <- file.path(sample.dir, "targets.R")
  cat(script, file=targets.R)
  Sys.chmod(targets.R, "0755")
  system(paste(
    "qsub",
    targets.R))
}

targets.RData.vec <- Sys.glob("labels/*/samples/*/*/targets.RData")

models.dt.list <- list()
iterations.dt.list <- list()
features.dt.list <- list()
for(targets.RData.i in seq_along(targets.RData.vec)){
  targets.RData <- targets.RData.vec[[targets.RData.i]]
  objs <- load(targets.RData)
  if(!"sample.features" %in% objs){
    cat(sprintf("%4d / %4d %s\n", targets.RData.i, length(targets.RData.vec), targets.RData))
    targets.R <- sub("Data$", "", targets.RData)
    system(targets.R)
    objs <- load(targets.RData)
    if(!"sample.features" %in% objs){
      stop("old script")
    }
  }
  models.dt.list[[targets.RData]] <- sample.models
  features.dt.list[[targets.RData]] <- sample.features
  iterations.dt.list[[targets.RData]] <- sample.iterations
}
models.dt <- do.call(rbind, models.dt.list)
features.dt <- do.call(rbind, features.dt.list)
iterations.dt <- do.call(rbind, iterations.dt.list)

fwrite(models.dt, "target.intervals.models.csv")
fwrite(features.dt, "target.intervals.features.csv")
fwrite(iterations.dt, "target.intervals.iterations.csv")


