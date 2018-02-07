#!/home/thocking/bin/Rscript
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V
#PBS -o /home/thocking/PeakSegFPOP/target.intervals.out
#PBS -e /home/thocking/PeakSegFPOP/target.intervals.err
#PBS -N Tint
library(data.table)
setwd("/home/thocking/PeakSegFPOP")
targets.R.vec <- Sys.glob("labels/*/samples/*/*/targets.R")

pattern <- paste0(
  "/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/",
  "(?<set_name>[^/]+)",
  "/samples/",
  "(?<sample_group>[^/]+)",
  "/",
  "(?<sample_id>[^/]+)",
  "/problems/",
  "(?<problem>",
  "(?<chrom>chr[^:]+)",
  ":",
  "(?<chromStart>[0-9]+)",
  "-",
  "(?<chromEnd>[0-9]+)",
  ")")
prob.dir <- "/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/H3K9me3_TDH_BP/samples/tcell/ERS358697/problems/chr4:68270000-75427379"
namedCapture::str_match_named(prob.dir, pattern, list(chromStart=as.integer, chromEnd=as.integer))
dt <- data.table(prob.dir)
change <- function(DT){
  ##df <- namedCapture::str_match_named(DT$prob.dir, pattern, list(chromStart=as.integer, chromEnd=as.integer))
  DT[, prob.dir := sub("/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/", "", prob.dir)]
}

models.dt.list <- list()
iterations.dt.list <- list()
features.dt.list <- list()
for(targets.i in seq_along(targets.R.vec)){
  targets.R <- targets.R.vec[[targets.i]]
  targets.RData <- sub("R$", "RData", targets.R)
  if(!file.exists(targets.RData)){
    cat(sprintf("%4d / %4d %s\n", targets.i, length(targets.R.vec), targets.R))
    system(targets.R)
  }
  tryCatch({
    objs <- load(targets.RData)
  }, error=function(e){
    cat(sprintf("%4d / %4d %s ERROR\n", targets.i, length(targets.R.vec), targets.R))
  })
  if(is.null(sample.models)){
    labels.bed <- sub("targets.R$", "labels.bed", targets.R)
    if(file.exists(labels.bed)){
      print(labels.bed)
    }
  }else{
    change(sample.models)
    change(sample.features)
    change(sample.iterations)
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


