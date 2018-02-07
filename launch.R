library(data.table)
problem.dir.vec <- Sys.glob("labels/H3K27ac-H3K4me3_TDHAM_BP/samples/*/*/problems/*")
prob.dt <- data.table(problem.dir=problem.dir.vec)
prob.dt[, target.tsv := file.path(problem.dir, "target.tsv")]
prob.dt[, labels.bed := file.path(problem.dir, "labels.bed")]
prob.dt[, target.tsv.exists := file.exists(target.tsv)]
prob.dt[, labels.bed.exists := file.exists(labels.bed)]
prob.dt[, table(labels.bed.exists, target.tsv.exists)]
sh.vec <- prob.dt[labels.bed.exists & !target.tsv.exists, paste0(target.tsv, ".sh")]

for(sh in sh.vec){
  cmd <- paste("qsub", sh)
  system(cmd)
}

set.name <- "ATAC_JV_adipose"
set.name <- "H3K27ac-H3K4me3_TDHAM_BP"
set.name <- "H3K4me3_TDH_immune"

for(set.name in dir("labels")){
  set.dir <- file.path("labels", set.name)

  ## This block creates a vector of shell script paths for every labeled
  ## problem in set.name
  labels.bed.vec <- Sys.glob(file.path(set.dir, "samples/*/*/problems/*/labels.bed"))
  sh.vec <- sub("labels.bed", "coverage.bedGraph.sh", labels.bed.vec)

  sample.dir.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*"))
  for(sample.dir in sample.dir.vec){
    labels <- fread(file.path(sample.dir, "labels.bed"))
    cmd <- paste("Rscript create_problems.R hg19_problems.bed", sample.dir)
    cat(cmd, "\n")
    system(cmd)
  }
}

##download label files.
for(set.name in dir("labels")){
  f <- file.path("labels", set.name, "labels", "original.txt")
  u <- sprintf(
    "http://cbio.mines-paristech.fr/~thocking/chip-seq-chunk-db/annotations/%s.txt",
    set.name)
  dir.create(dirname(f), showWarnings=FALSE,recursive=TRUE)
  download.file(u, f)
}

## This block creates a vector of shell script paths for every labeled
## problem which does not have a corresponding target_models.tsv
labels.bed.vec <- Sys.glob(file.path("labels", set.name, "samples/*/*/problems/*/labels.bed"))
models.vec <- sub("labels.bed", "target_models.tsv", labels.bed.vec)
has.models <- file.exists(models.vec)
sh.vec <- sub("labels.bed", "target.tsv.sh", labels.bed.vec[!has.models])

for(sh in sh.vec){
  cmd <- paste("qsub", sh)
  system(cmd)
}

## Train a model for each data set.
for(set.name in dir("labels")){
  set.dir <- file.path("labels", set.name)
  samples.dir <- file.path(set.dir, "samples")
  model.RData <- file.path(samples.dir, "model.RData")
  cmd <- paste("Rscript train_model.R", samples.dir, model.RData)
  system(cmd)
}

## This code block finds all .out files with non-zero exit codes.
library(namedCapture)
pattern <- paste0(
  "(?<bedGraph>.*.out)",
  ".*?",
  "(?<code>-?[0-9]+)$")
"grep -nH -e 'Begin PBS Prologue Mon' labels/H3K36me3_AM_immune/*/problems/*/coverage.bedGraph.out"
code.vec <- system("grep code labels/*/problems/*/jointProblems.bed.out", intern=TRUE)
meanings <- c(
  "0"="success",
  "1"="failure",
  ## both of the following are possible memory errors.
  "-11"="walltime exceeded limit",
  "265"="SIGKILL (kill -9)",
  "271"="SIGTERM (canceljob)")
code.mat <- str_match_named(code.vec, pattern)
code.int <- as.integer(code.mat[, "code"])
table(code.int)
failed <- code.mat[code.int!=0, "bedGraph"]
sh.vec <- sub(".out", ".sh", failed)
err.vec <- sub(".out", ".err", failed)

system(paste("grep 'Begin PBS Pro'", paste(failed, collapse=" "), "|grep Wed"))

cmd <- paste("tail -n 50", paste(paste0(code.mat[code.int==1, "bedGraph"], ".err"), collapse=" "))
system(cmd)

## This block finds all jointProblems.bed.sh files which have no
## corresponding out file.
notDone <- function(glob.sh){
  library(data.table)
  library(namedCapture)
  code.pattern <- paste0(
    "(?<out>.*.out)",
    ".*?",
    "(?<code>-?[0-9]+)$")
  stopifnot(is.character(glob.sh))
  stopifnot(length(glob.sh)==1)
  sh.vec <- Sys.glob(glob.sh)
  glob.out <- sub("sh$", "out", glob.sh)
  stopifnot(glob.sh != glob.out)
  grep.cmd <- paste("grep code", glob.out)
  grep.lines <- system(grep.cmd, intern=TRUE)
  results <- data.table(code=NA_character_, sh=sh.vec)
  setkey(results, sh)
  if(length(grep.lines)){
    code.mat <- str_match_named(grep.lines, code.pattern)
    results[sub("out$", "sh", code.mat[, "out"]), code := code.mat[, "code"] ]
  }
  print(table(results$code, useNA="always"))
  results[order(code),]
}

set.name <- "H3K27me3_RL_cancer"

bigwig.vec <- notDone(file.path(
  "labels", set.name, "samples/*/*/coverage.bigWig.sh"))

target.vec <- notDone(file.path(
  "labels", set.name, "samples/*/*/problems/*/target.tsv.sh"))
target.vec[, labeled := file.exists(sub("target.tsv.sh", "labels.bed", sh))]
target.vec[labeled==TRUE, table(code, useNA="always")]
target.vec[labeled==TRUE & is.na(code),]
glob.out <- sub("sh$", "out", target.vec[labeled==TRUE, sh])

## each jointProblems.bed.sh file does peak prediction for one
## problem, for all samples.
glob.sh <- file.path(
  "labels", set.name, "problems/*/jointProblems.bed.sh")
glob.out <- sub("sh$", "out", glob.sh)
jointProblems <- notDone(glob.sh)

pos.pattern <- paste0(
  "(?<chrom>chr.*?)",
  ":",
  "(?<chromStart>[0-9]+)",
  "-",
  "(?<chromEnd>[0-9]+)")
pos.mat <- jointProblems[, str_match_named(sh, pos.pattern)]
jointProblems[, chromStart := as.integer(pos.mat[, "chromStart"])]
jointProblems[, chromEnd := as.integer(pos.mat[, "chromEnd"])]
jointProblems[, chromBases := chromEnd-chromStart]

walltimes <- function(glob.out){
  cmd <- paste("grep ^Resources", glob.out)
  lines.dt <- fread(cmd, sep="=", header=FALSE)
  file.out.vec <- 
  pattern <- paste0(
    "(?<hours>[0-9]+)",
    ":",
    "(?<minutes>[0-9]+)",
    ":",
    "(?<seconds>[0-9]+)")
  match.df <- str_match_named(lines.dt$V6, pattern, list(
    hours=as.integer,
    minutes=as.integer,
    seconds=as.integer))
  data.table(
    file.out=sub("out:.*", "out", lines.dt$V1),
    minutes=with(match.df, hours*60 + minutes + seconds/60))
}
walltime.dt <- walltimes(glob.out)

## sizes
jointProblems.all <- jointProblems[, {
  bed <- sub(".sh$", "", sh)
  if(file.exists(bed))fread(bed)
}, by=.(sh, chromBases)]
setnames(jointProblems.all, c("sh", "chromBases", "chrom", "problemStart", "problemEnd"))
jointProblems.all[, experiment := basename(sub("_.*", "", sh))]
jointProblems.all[, bases := problemEnd-problemStart]
jointProblems.all[, list(
  problems=.N,
  median.bases=median(bases),
  mean.bases=mean(bases),
  max.bases=max(bases)
  ), by=experiment]
jointProblems.all[order(-bases),]

jointProblems.all[, quantile(bases, probs=seq(0, 1, by=0.1)), by=experiment]
library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_histogram(aes(log10(bases)), data=jointProblems.all)

(joint.models <- notDone("labels/*/joint.model.RData.sh"))

## This block finds all jobPeaks.sh files which have no
jpeaks <- notDone(file.path(
  "labels",
  set.name,
  "jobs/*/jobPeaks.sh"))

## Which have no corresponding peaks.bed file?
sh.vec <- jpeaks[!file.exists(sub(".sh$", "", sh)), sh]

sh.vec <- jpeaks[is.na(code), sh]
sh.vec <- jpeaks[code!=0, sh]
sh.vec <- jpeaks[code==-11, sh]
sh.vec <- jpeaks[code==1, sh]

for(sh in sh.vec){
  cmd <- paste("qsub", sh)
  cat(cmd, "\n")
  system(cmd)
}

## Edit the following definition to reflect your cluster
## configuration.
PBS.header <- "#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V"
jprobs.not.done <- data.table(
  jprobs.bed=sub("peaks.bed.sh$", "jointProblems.bed", sh.vec))[, fread(jprobs.bed), by=jprobs.bed]
setnames(jprobs.not.done, c("jprobs.bed", "chrom", "problemStart", "problemEnd"))
n.jobs <- 300
jprobs.not.done[, job := 1:n.jobs]
jprobs.not.done[, table(job)]
data.dir <- "labels/H3K27ac-H3K4me3_TDHAM_BP"
coverage.bigWig.vec <- Sys.glob(file.path(
  data.dir, "samples", "*", "*", "coverage.bigWig"))
jobs.dir <- file.path(data.dir, "jobs")

Rscript <- function(...){
  code <- sprintf(...)
  if(any(grepl("'", code))){
    print(code)
    stop("there can not be any ' in code")
  }
  sprintf("Rscript -e '%s'", code)
}
jprobs.not.done[, {
  job.dir <- normalizePath(file.path(jobs.dir, job), mustWork=TRUE)
  dir.create(job.dir, showWarnings=FALSE, recursive=TRUE)
  fwrite(
    .SD[,.(chrom, problemStart, problemEnd, basename(dirname(jprobs.bed)))],
    file.path(job.dir, "jobProblems.bed"),
    sep="\t", col.names=FALSE)
  script.txt <- paste0(
    PBS.header, "
#PBS -o ", file.path(job.dir, "jobPeaks"), ".out
#PBS -e ", file.path(job.dir, "jobPeaks"), ".err
#PBS -N Job", job, "
", Rscript('PeakSegJoint::problem.joint.predict.job("%s")', job.dir))
  writeLines(script.txt, file.path(job.dir, "jobPeaks.sh"))
}, by=job]

jobstat <- notDone(file.path(
  "labels", set.name, "jobs/*/jobPeaks.sh"))



## Look for overlapping joint problems.
library(data.table)
dt <- fread("cat labels/H3K36me3_TDH_immune/problems/*/jointProblems.bed")
setnames(dt, c("chrom", "problemStart", "problemEnd"))
setkey(dt, chrom, problemStart, problemEnd)
dt[, next.start := c(problemStart[-1], NA), by=chrom]
dt[, dist.to.next := next.start - problemEnd]
dt[dist.to.next < 0,][order(chrom, problemStart),]


bg.out.vec <- Sys.glob("labels/*/samples/*/*/problems/*/coverage.bedGraph.out")

for(bg.out in bg.out.vec){
  all.loss <- fread(paste("grep '^[ 0-9]*:'", bg.out, "|sed 's/:/ /'"))
  if(length(all.loss)==5){
    setnames(all.loss, c("row", "penalty", "peaks", "fp", "fn"))
  }else{
    stop("unknown columns")
  }
  all.loss[, cum.zeros := cumsum(penalty==0)]
  all.loss[cum.zeros==max(cum.zeros),]
}

for(bg.i in seq_along(bg.out.vec)){
  bg.out <- bg.out.vec[[bg.i]]
  bg.sh <- sub("out$", "sh", bg.out)
  cat(sprintf(
    "%4d / %4d %s\n", bg.i, length(bg.out.vec),
    bg.sh))
  cmd <- paste("bash", bg.sh)
  system(cmd)
}
