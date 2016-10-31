arg.vec <- "test/demo"

arg.vec <- commandArgs(trailingOnly=TRUE)

library(xtable)
library(ggplot2)
library(animint)
library(WeightedROC)
library(ggdendro)
library(namedCapture)
library(data.table)
library(coseg)

if(length(arg.vec) != 1){
  stop("usage: Rscript plot_all.R project_dir")
}
set.dir <- normalizePath(arg.vec, mustWork=TRUE)

orderChrom <- function(chrom.vec, ...){
  stopifnot(is.character(chrom.vec))
  chr.pattern <- paste0(
    "chr",
    "(?<before>[^_]+)",
    "(?<after>_.*)?")
  value.vec <- unique(chrom.vec)
  chr.mat <- str_match_named(value.vec, chr.pattern)
  did.not.match <- is.na(chr.mat[, 1])
  if(any(did.not.match)){
    print(value.vec[did.not.match])
    stop("chroms did not match ", chr.pattern)
  }
  rank.vec <- order(
    suppressWarnings(as.numeric(chr.mat[, "before"])),
    chr.mat[, "before"],
    chr.mat[, "after"])
  names(rank.vec) <- value.vec
  order(rank.vec[chrom.vec], ...)
} 

## Plot each labeled chunk.
chunk.dir.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*"))
mclapply.or.stop(chunk.dir.vec, function(chunk.dir){
  PeakSegJoint::problem.joint.plot(chunk.dir)
})

unsorted.problems <- fread(file.path(set.dir, "problems.bed"))
setnames(unsorted.problems, c("chrom", "problemStart", "problemEnd"))
chr.pattern <- paste0(
  "chr",
  "(?<before>[^_]+)",
  "(?<after>_.*)?")
chr.mat <- str_match_named(unsorted.problems$chrom, chr.pattern)
problems <- unsorted.problems[order(
  suppressWarnings(as.numeric(chr.mat[, "before"])),
  chr.mat[, "before"],
  chr.mat[, "after"],
  problemStart),]
problems[, problem.name := sprintf(
  "%s:%d-%d", chrom, problemStart, problemEnd)]
problems[, separate.problem := factor(problem.name, problem.name)]

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
  separate.problem <- factor(basename(prob.dir), problems$problem.name)
  input.pred.list[[loss.tsv]] <- data.table(
    prob.peaks[1, .(chrom, peakStart, peakEnd, peak.name)],
    n.Input,
    n.samples=nrow(prob.peaks),
    n.groups=nrow(group.dt),
    loss.diff,
    sample.counts=paste(group.vec, collapse="/"),
    separate.problem,
    joint.problem=basename(joint.prob.dir))
  joint.peaks.dt.list[[loss.tsv]] <- data.table(
    separate.problem, prob.peaks)
}
joint.peaks.dt <- do.call(rbind, joint.peaks.dt.list)
input.pred <- do.call(rbind, input.pred.list)[orderChrom(chrom, peakStart),]
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
##plot(tree, hang=-1)

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

if(nrow(input.labels)){
  setkey(input.pred, chrom, peakStart, peakEnd)
  labeled.input <- foverlaps(input.pred, input.labels, nomatch=0L)
  thresh.dt <- labeled.input[, data.table(WeightedROC(
    n.Input, ifelse(prop.noPeaks==0, 1, -1)))]
  thresh.best <- thresh.dt[which.min(FP+FN),]
  ## threshold is smallest n.Input that is classified as non-specific.
  input.pred[, specificity := ifelse(
    n.Input >= thresh.best$threshold, "non-specific", "specific")]
}else{
  input.pred[, specificity := "unknown"]
}

gg.tree <- dendro_data(tree)
gg.tree$labels$sample.path <- gg.tree$labels$label
gg.tree$labels$sample.group <- sub("/.*", "", gg.tree$labels$label)
gg.tree$labels$sample.id <- sub(".*/", "", gg.tree$labels$label)
rownames(gg.tree$labels) <- gg.tree$labels$label
max.dist <- max(gg.tree$segments$y)
max.dist.neg <- -max.dist*3
problems[, dist := (.N:1)*max.dist.neg/.N]
problem.width <- diff(problems$dist[1:2])/2
setkey(problems, separate.problem)
sample.problem.peaks <- joint.peaks.dt[, list(
  peaks=.N
), by=.(sample.path, separate.problem)]
setkey(sample.problem.peaks, separate.problem)
problem.peak.show <- problems[sample.problem.peaks]
problem.peak.show[, sample.pos := gg.tree$labels[paste(sample.path), "x"]]
group.means <- data.table(gg.tree$labels)[, list(
  mean.x=mean(x)
), by=sample.group]
peak.samples.counts <- sample.problem.peaks[, list(
  samples=.N,
  peak.samples=sum(peaks)
  ), by=separate.problem]
peaks.counts <- input.pred[, list(peaks=.N), by=separate.problem]
problems[, peaks := 0L]
problems[, samples := 0L]
problems[, peak.samples := 0L]
problems[peaks.counts, peaks := peaks.counts$peaks]
problems[peak.samples.counts, peak.samples := peak.samples.counts$peak.samples]
problems[peak.samples.counts, samples := peak.samples.counts$samples]

ggplot()+
  theme_grey()+
  geom_tile(aes(
    dist, sample.pos, fill=peaks),
    data=problem.peak.show)+
  scale_color_discrete(
    breaks=group.means[order(-mean.x), sample.group])+
  scale_fill_continuous(
    low="white", high="red")+
  geom_segment(aes(
    y, x, yend=xend, xend=yend),
    data=gg.tree$segments)+
  geom_text(aes(
    y, x, label=sample.id, color=sample.group),
    hjust=0,
    vjust=-0.5, 
    data=gg.tree$labels)+
  scale_x_continuous(
    "separate problem / number of peaks")+
  scale_y_continuous("sample", breaks=NULL)

ggplot()+
  scale_fill_continuous(low="white", high="black")+
  geom_point(aes(
    peakStart/1e3, chrom, color=specificity, fill=log10(loss.diff)),
    shape=21,
    data=input.pred)

n.samples <- nrow(gg.tree$labels)
h.pixels <- (n.samples+5)*15
chrom.limits <- problems[, list(
  min.dist=min(dist),
  max.dist=max(dist)
), by=chrom]
vline.dt <- chrom.limits[, data.table(
  dist=unique(c(min.dist-problem.width, max.dist+problem.width)))]
input.pred[, peakBases := peakEnd - peakStart]
joint.peaks.dt[, peakBases := peakEnd - peakStart]
viz <- list(
  title="PeakSegFPOP + PeakSegJoint predictions",
  tree=ggplot()+
    geom_text(aes(
      min.dist-problem.width, n.samples+0.5,
      label=sub("chr", "", chrom)),
      size=11,
      hjust=0,
      data=chrom.limits)+
    coord_cartesian(
      xlim=c(max.dist.neg-problem.width, max.dist),
      ylim=c(0.5, n.samples+1.5),
      expand=FALSE)+
    geom_vline(aes(
      xintercept=dist),
      size=0.25,
      data=vline.dt)+
    theme_grey()+
    theme_animint(width=1000, height=h.pixels)+
    geom_tallrect(aes(
      xmin=dist-problem.width,
      xmax=dist+problem.width,
      clickSelects=separate.problem),
      alpha=0.5,
      size=0.5,
      data=problems)+
    geom_tile(aes(
      dist, sample.pos, fill=peaks,
      clickSelects=separate.problem,
      tooltip=paste0(
        peaks, " peak",
        ifelse(peaks==1, "", "s"),
        " for ", sample.path,
        " in ", separate.problem
      )),
      size=0.5,
      data=problem.peak.show)+
    scale_color_discrete(
      breaks=group.means[order(-mean.x), sample.group])+
    scale_fill_continuous(
      low="white", high="red",
      breaks=sort(
        unique(as.integer(quantile(problem.peak.show$peaks))),
        decreasing=TRUE))+
    geom_segment(aes(
      y, x, yend=xend, xend=yend),
      size=1,
      data=gg.tree$segments)+
    geom_text(aes(
      y, x,
      key=sample.path,
      label=sample.id, color=sample.group),
      hjust=0,
      size=11,
      data=gg.tree$labels)+
    geom_point(aes(
      y, x,
      key=sample.path,
      color=sample.group),
      data=gg.tree$labels)+
    scale_x_continuous(
      paste(
        "Left heatmap: select genomic region (separate segmentation problem),",
        "Right tree: distance = number of peaks different"),
      breaks=unique(as.integer(seq(0, max.dist, l=6)))
    )+
    scale_y_continuous("sample", breaks=NULL),
  genome=ggplot()+
    theme_grey()+
    theme_animint(
      width=1000,
      height=h.pixels,
      update_axes=c("x"))+
    geom_vline(aes(
      showSelected=separate.problem,
      xintercept=peakStart/1e3,
      color=loss.diff,
      linetype=specificity),
      size=3,
      data=input.pred)+
    geom_tallrect(aes(
      xmin=problemStart/1e3,
      xmax=problemEnd/1e3,
      showSelected=separate.problem,
      tooltip=paste(
        "click to zoom to problem",
        separate.problem),
      href=paste0(
        "http://genome.ucsc.edu/cgi-bin/hgTracks?position=",
        separate.problem)),
      fill="white",
      color=NA,
      alpha=0.5,
      data=problems)+
    ## geom_vline(aes(
    ##   showSelected=separate.problem,
    ##   xintercept=peakStart/1e3,
    ##   size=loss.diff,
    ##   tooltip=paste(
    ##     "click to zoom to peak",
    ##     peak.name,
    ##     sample.counts,
    ##     "Poisson loss difference over model with no peaks = ",
    ##     loss.diff),
    ##   href=sprintf(
    ##     "http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s:%d-%d",
    ##     chrom, peakStart-peakBases, peakEnd+peakBases),
    ##   linetype=specificity),
    ##   color="grey",
    ##   data=input.pred)+
    ## scale_size_continuous(range=c(1, 9))+
    scale_linetype_manual(values=c(
      specific="solid",
      "non-specific"="dotted",
      unknown="dashed"
    ))+
    guides(fill="none")+
    scale_color_continuous(
      low="white", high="red")+
    scale_x_continuous("position on chromosome (kb = kilo bases)")+
    scale_y_discrete("sample"),
    ## geom_segment(aes(
    ##   peakStart/1e3, sample.path,
    ##   xend=peakEnd/1e3, yend=sample.path,
    ##   color=sample.group,
    ##   showSelected2=sample.group,
    ##   showSelected=separate.problem),
    ##   data=joint.peaks.dt)+
  selector.types=list())

figure.png.vec <- Sys.glob(file.path(
  set.dir, "problems", "*", "chunks", "*", "figure-predictions-zoomout.png"))
if(0 == length(figure.png.vec)){
  chunk.info <- data.table()
  chunks.html <- ""
}else{
  relative.vec <- sub("/", "", sub(set.dir, "", figure.png.vec))
  zoomin.png.vec <- sub("-zoomout", "", relative.vec)
  chunk.dir.vec <- dirname(zoomin.png.vec)
  chunks.dir.vec <- dirname(chunk.dir.vec)
  separate.prob.dir.vec <- dirname(chunks.dir.vec)
  g.pos.pattern <- paste0(
    "(?<chrom>chr.+?)",
    ":",
    "(?<chromStart>[0-9 ,]+)",
    "-",
    "(?<chromEnd>[0-9 ,]+)")
  pos2df <- function(path.vec){
    problem <- basename(path.vec)
    df <- str_match_named(
      problem,
      g.pos.pattern,
      list(
        chromStart=as.integer,
        chromEnd=as.integer))
    data.frame(problem, df)
  }
  chunk.info <- data.table(
    separate=pos2df(separate.prob.dir.vec),
    chunk=pos2df(chunk.dir.vec),
    zoomin.png=zoomin.png.vec)
  chunk.info[, image := sprintf('
<a href="%s">
  <img src="%s" />
</a>
', zoomin.png.vec, sub(".png$", "-thumb.png", zoomin.png.vec))]
  chunk.info[, chunk := sprintf({
    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
  }, chunk.problem, chunk.problem)]
  chunk.counts <- chunk.info[, list(chunks=.N), by=separate.problem]
  problems[, labeled.chunks := 0L]
  setkey(chunk.counts, separate.problem)
  problems[chunk.counts, labeled.chunks := chunk.counts$chunks]
  viz$genome <- viz$genome+
    geom_tallrect(aes(
      xmin=chunk.chromStart/1e3,
      xmax=chunk.chromEnd/1e3,
      tooltip=paste(
        "click to plot chunk",
        chunk.problem),
      href=file.path("..", zoomin.png),
      showSelected=separate.problem),
      color=NA,
      alpha=0.2,
      fill="yellow",
      data=chunk.info)
  chunks.xt <- xtable(chunk.info[, .(chunk, image)])
  chunks.html <- print(chunks.xt, type="html", sanitize.text.function=identity)
}

viz$genome <- viz$genome+
  geom_point(aes(
    peakStart/1e3, sample.path,
    fill=sample.group,
    key=paste(chrom, peakStart, sample.path),
    tooltip=paste(
      "click to zoom to peak",
      peak.name),
    href=sprintf(
      "http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s:%d-%d",
      chrom, peakStart-peakBases, peakEnd+peakBases),
    showSelected2=sample.group,
    showSelected=separate.problem),
    color="black",
    size=3.5,
    data=joint.peaks.dt)
animint2dir(viz, file.path(set.dir, "figure-genome"))

label.counts <- all.labels[, list(
  samples=.N
  ), by=.(chrom, labelStart, labelEnd)]
label.counts[, labelStart1 := labelStart + 1L]
problems[, problemStart1 := problemStart + 1L]
setkey(label.counts, chrom, labelStart1, labelEnd)
setkey(problems, chrom, problemStart1, problemEnd)
labels.in.problems <- foverlaps(label.counts, problems, nomatch=0L)
label.problem.counts <- labels.in.problems[, list(
  labels=sum(samples),
  labeled.regions=.N
  ), by=separate.problem]
problems[, labels := 0L]
setkey(label.problem.counts, separate.problem)
setkey(problems, separate.problem)
problems[label.problem.counts, labels := label.problem.counts$labels]
problems[, labeled.regions := 0L]
problems[label.problem.counts, labeled.regions := label.problem.counts$labeled.regions]
problems[, problem := sprintf({
  '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s">%s</a>'
}, separate.problem, separate.problem)]

problems.xt <- xtable(problems[, .(problem, samples, peaks, peak.samples, labeled.chunks, labeled.regions, labels)])
problems.html <- print(problems.xt, type="html", sanitize.text.function=identity)
html.vec <- c(
  '<title>PeakSegFPOP + PeakSegJoint predictions</title>',
  '<h2>Output: predicted peaks</h2>',
  '<ul>',
  sprintf('<li>
%d genomic regions were found to have
at least one sample with a peak
(<a href="peaks_summary.bigBed">peaks_summary.bigBed</a>,
 <a href="peaks_summary.tsv">peaks_summary.tsv</a>).
</li>', nrow(input.pred)),
  sprintf('<li>
%d peaks were detected overall, across %d samples
(<a href="peaks_matrix.tsv">peaks_matrix.tsv</a>).
</li>',
    nrow(joint.peaks.dt),
    length(unique(sample.problem.peaks$sample.path))),
  '<li>
<a href="figure-genome/index.html">Interactive heatmap, cluster tree, peak predictions</a>.
</li>',
  '</ul>',
  '<h2>Input: labeled genomic windows</h2>',
  sprintf(
    '%d labeled genomic windows (chunks) were defined in <a href="labels">labels/*.txt</a>',
    nrow(chunk.info)),
  chunks.html,
  '<h2>Input: genomic segmentation problems</h2>',
  sprintf(
    '%d problems were defined in <a href="problems.bed">problems.bed</a>',
    nrow(problems)),
  problems.html
  )
writeLines(html.vec, file.path(set.dir, "index.html"))

## Write peaks_summary.bigBed 
peak.mat.dt <- data.table(
  peak=colnames(peak.mat),
  ifelse(t(peak.mat), 1, 0))
fwrite(
  peak.mat.dt,
  file.path(set.dir, "peaks_matrix.tsv"),
  sep="\t")
fwrite(
  input.pred,
  file.path(set.dir, "peaks_summary.tsv"),
  sep="\t")
peaks.bed <- file.path(set.dir, "peaks_summary.bed")
bed.dt <- input.pred[, .(chrom, peakStart, peakEnd, sample.counts)]
fwrite(
  bed.dt,
  peaks.bed,
  sep="\t",
  col.names=FALSE)
sizes.dt <- problems[, list(chromEnd=max(problemEnd)), by=chrom]
sizes.tsv <- file.path(set.dir, "chrom_sizes.tsv")
fwrite(
  sizes.dt,
  sizes.tsv,
  sep="\t",
  col.names=FALSE)
cmd <- paste(
  "bedToBigBed",
  peaks.bed, sizes.tsv,
  file.path(set.dir, "peaks_summary.bigBed"))
status <- system(cmd)
if(status != 0){
  stop("error code ", status, " for command\n", cmd)
}
## TODO: labels.bigBed and track hub.

