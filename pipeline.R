arg.vec <- "test/input"

arg.vec <- commandArgs(trailingOnly=TRUE)

set.dir <- normalizePath(arg.vec, mustWork=TRUE)

system.or.stop <- function(cmd){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
}

## First convert labels.
convert.cmd <- paste("Rscript convert_labels.R", set.dir)
system.or.stop(convert.cmd)

## Create problems for each sample.
create.cmd <- paste("Rscript create_problems_all.R", set.dir)
system.or.stop(create.cmd)

## Compute target interval for each problem.
samples.dir <- file.path(set.dir, "samples")
labels.bed.vec <- Sys.glob(file.path(
  samples.dir, "*", "*", "problems", "*", "labels.bed"))
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  target.cmd <- paste("Rscript compute_coverage_target.R", sample.dir)
  system.or.stop(target.cmd)
}

## Train single-sample model.
train.cmd <- paste("Rscript train_model.R", set.dir)
system.or.stop(train.cmd)

## Single-sample prediction and peak clustering, one job for each
## problem.
sh.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems.bed.sh"))
for(sh in sh.vec){
  predict.cmd <- paste("bash", sh)
  system.or.stop(predict.cmd)
}

## Compute target intervals for multi-sample problems... does not take
## much time, TODO combine with train_model_joint.R step
labels.tsv.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "labels.tsv"))
for(labels.tsv in labels.tsv.vec){
  target.cmd <- paste("Rscript compute_joint_target.R", dirname(labels.tsv))
  system.or.stop(target.cmd)
}
## Train joint model.
train.cmd <- paste("Rscript train_model_joint.R", set.dir)
system.or.stop(train.cmd)

## Joint prediction.
joint.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*"))
joint.model.RData <- file.path(set.dir, "joint.model.RData")
for(joint.dir in joint.dir.vec){
  predict.cmd <- paste(
    "Rscript predict_problem_joint.R",
    joint.model.RData, joint.dir)
  system.or.stop(predict.cmd)
}

## Plots.
chunk.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*"))
for(chunk.dir in chunk.dir.vec){
  plot.cmd <- paste(
    "Rscript plot_chunk.R",
    chunk.dir)
  system.or.stop(plot.cmd)
}

## TODO: make this a separate script, summarize other info like
## predicted peaks, specific peaks, etc.
loss.tsv.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "jointProblems", "*", "loss.tsv"))
joint.peaks.dt.list <- list()
input.pred.list <- list()
for(loss.tsv in loss.tsv.vec){
  loss.diff <- scan(loss.tsv, quiet=TRUE)
  joint.prob.dir <- dirname(loss.tsv)
  jointProblems.dir <- dirname(joint.prob.dir)
  prob.dir <- dirname(jointProblems.dir)
  prob.peaks <- fread(file.path(joint.prob.dir, "peaks.bed"))
  setnames(prob.peaks, c(
    "chrom", "peakStart", "peakEnd", "sample.path", "mean"))
  prob.peaks[, sample.id := sub(".*/", "", sample.path)]
  prob.peaks[, sample.group := sub("/.*", "", sample.path)]
  prob.peaks[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]
  group.dt <- prob.peaks[, list(samples=.N), by=sample.group]
  group.vec <- group.dt[order(sample.group), paste0(
    sample.group, ":", samples)]
  n.Input <- sum(prob.peaks$sample.group=="Input")
  input.pred.list[[loss.tsv]] <- data.table(
    prob.peaks[1, .(chrom, peakStart, peakEnd, peak.name)],
    n.Input,
    n.samples=nrow(prob.peaks),
    n.groups=nrow(group.dt),
    loss.diff,
    sample.counts=paste(group.vec, collapse="/"),
    separate.problem=basename(prob.dir),
    joint.problem=basename(joint.prob.dir))
  joint.peaks.dt.list[[loss.tsv]] <- prob.peaks
}
joint.peaks.dt <- do.call(rbind, joint.peaks.dt.list)
input.pred <- do.call(rbind, input.pred.list)
sample.dir.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*"))
sample.path.vec <- sub(".*samples/", "", sample.dir.vec)
peak.mat <- matrix(
  FALSE,
  length(sample.path.vec),
  nrow(input.pred),
  dimnames=list(
    sample=sample.path.vec,
    peak=input.pred$peak.name))
i.mat <- joint.peaks.dt[, cbind(sample.path, peak.name)]
peak.mat[i.mat] <- TRUE
d.mat <- dist(peak.mat, "manhattan")
tree <- hclust(d.mat, method="average")
plot(tree, hang=-1)

labels.bed.vec <- Sys.glob(file.path(
  set.dir, "samples", "*", "*", "labels.bed"))
all.labels.list <- list()
for(labels.bed in labels.bed.vec){
  sample.dir <- dirname(labels.bed)
  sample.id <- basename(sample.dir)
  group.dir <- dirname(sample.dir)
  sample.group <- basename(group.dir)
  sample.labels <- fread(labels.bed)
  setnames(sample.labels, c("chrom", "labelStart", "labelEnd", "annotation"))
  all.labels.list[[labels.bed]] <- data.table(
    sample.id, sample.group, sample.labels)
}
all.labels <- do.call(rbind, all.labels.list)
input.labels <- all.labels[sample.group=="Input", list(
  prop.noPeaks=mean(annotation=="noPeaks")
), by=.(chrom, labelStart, labelEnd)]
setkey(input.labels, chrom, labelStart, labelEnd)

setkey(input.pred, chrom, peakStart, peakEnd)
labeled.input <- foverlaps(input.pred, input.labels, nomatch=0L)
library(WeightedROC)
thresh.dt <- labeled.input[, data.table(WeightedROC(
  n.Input, ifelse(prop.noPeaks==0, 1, -1)))]
thresh.best <- thresh.dt[which.min(FP+FN),]
## threshold is smallest n.Input that is classified as non-specific.
input.pred[, specificity := ifelse(
  n.Input >= thresh.best$threshold, "non-specific", "specific")]
library(ggdendro)
gg.tree <- dendro_data(tree)
gg.tree$labels$sample.path <- gg.tree$labels$label
gg.tree$labels$sample.group <- sub("/.*", "", gg.tree$labels$label)
gg.tree$labels$sample.id <- sub(".*/", "", gg.tree$labels$label)

ggplot()+
  geom_segment(aes(
    y, x, yend=xend, xend=yend),
    data=gg.tree$segments)+
  geom_text(aes(
    y, x, label=sample.id, color=sample.group),
    hjust=1,
    data=gg.tree$labels)+
  scale_x_continuous(
    "number of peaks",
    limits=max.dist*c(-0.2, 1))+
  scale_y_continuous("sample", breaks=NULL)

max.dist <- max(gg.tree$segments$y)
ggplot()+
  scale_fill_continuous(low="white", high="black")+
  geom_point(aes(
    peakStart/1e3, chrom, color=specificity, fill=log10(loss.diff)),
    shape=21,
    data=input.pred)

viz <- list(
  tree=ggplot()+
    geom_segment(aes(
      y, x, yend=xend, xend=yend),
      data=gg.tree$segments)+
    geom_text(aes(
      y, x, label=sample.id, color=sample.group,
      clickSelects=sample.path),
      hjust=1,
      data=gg.tree$labels)+
    guides(color="none")+
    scale_x_continuous(
      "number of peaks",
      limits=max.dist*c(-0.2, 1))+
    scale_y_continuous("sample", breaks=NULL),
  genome=ggplot()+
    scale_color_manual(values=c(specific=NA, "non-specific"="black"))+
    geom_point(aes(
      peakStart/1e3, chrom,
      color=specificity,
      fill=sample.group,
      showSelected=sample.path),
      shape=21,
      alpha=0.2,
      data=joint.peaks.dt),
  selector.types=list(sample.path="multiple"))
library(animint)
animint2dir(viz, "test/input/figure-genome")

figure.png.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*", "figure-predictions-zoomout.png"))
relative.vec <- sub("/", "", sub(set.dir, "", figure.png.vec))
a.vec <- sprintf('
<a href="%s">
  <img src="%s" />
</a>
', relative.vec, sub("-zoomout", "", relative.vec))
writeLines(a.vec, file.path(set.dir, "index.html"))
