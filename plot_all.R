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

## The UCSC links in the HTML tables will be zoomed out from the peak
## this number of times.
zoom.out.times <- 10
zoom.factor <- (zoom.out.times-1)/2

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
LAPPLY <- lapply
LAPPLY <- mclapply.or.stop
LAPPLY(chunk.dir.vec, function(chunk.dir){
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

joint.glob <- file.path(
  set.dir, "jobs", "*")
jobPeaks.RData.dt <- data.table(
  jobPeaks.RData=Sys.glob(file.path(joint.glob, "jobPeaks.RData")))
if(nrow(jobPeaks.RData.dt)==0){
  stop(
    "no predicted joint peaks found; to do joint peak prediction run ",
    file.path(joint.glob, "jobPeaks.sh"))
}
cat(
  "Reading predicted peaks in",
  nrow(jobPeaks.RData.dt),
  "jobPeaks.RData files.\n",
  sep=" ")
##load all jobPeaks files.
jobPeaks <- jobPeaks.RData.dt[, {
    load(jobPeaks.RData)
    jobPeaks
}, by=jobPeaks.RData]
jobPeaks[, peak.name := sprintf("%s:%d-%d", chrom, peakStart, peakEnd)]

## Count total samples using directories.
sample.dir.vec <- Sys.glob(file.path(set.dir, "samples", "*", "*"))
sample.path.vec <- sub(".*samples/", "", sample.dir.vec)
sample.group.tab <- table(sub("/.*", "", sample.path.vec))
sample.group.totals <- data.table(
  sample.group=names(sample.group.tab),
  samples.total=as.integer(sample.group.tab))
setkey(sample.group.totals, sample.group)

if(FALSE){ #old slow code
  joint.glob <- file.path(
    set.dir, "problems", "*")
  loss.tsv.vec <- Sys.glob(file.path(
    joint.glob, "jointProblems", "*", "loss.tsv"))
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
    group.dt <- prob.peaks[, list(
      samples.with.peaks=.N
    ), by=sample.group]
    setkey(group.dt, sample.group)
    all.group.counts <- group.dt[sample.group.totals]
    all.group.counts[is.na(samples.with.peaks), samples.with.peaks := 0L]
    all.group.counts[, samples.without.peaks := samples.total - samples.with.peaks]
    all.group.counts[, samples.prop := samples.with.peaks / samples.total]
    count.mat <- all.group.counts[, rbind(samples.with.peaks, samples.without.peaks)]
    exact <- fisher.test(count.mat)
    ##group.ord <- all.group.counts[order(samples.prop),]
    group.ord <- all.group.counts[order(sample.group), list(
      groups=paste(sample.group, collapse=",")
    ), by=samples.prop][order(samples.prop),]
    group.vec <- group.dt[order(sample.group), paste0(
      sample.group, ":", samples.with.peaks)]
    n.Input <- sum(prob.peaks$sample.group=="Input")
    separate.problem <- factor(basename(prob.dir), problems$problem.name)
    most1 <- group.ord[.N,]
    least1 <- group.ord[1,]
    if(nrow(group.ord)==1){
      most2 <- least2 <- data.table(samples.prop=NA, groups=NA)
    }else{
      most2 <- group.ord[.N-1,]
      least2 <- group.ord[2,]
    }
    input.pred.list[[loss.tsv]] <- data.table(
      prob.peaks[1, .(chrom, peakStart, peakEnd, peak.name)],
      n.Input,
      n.samples=nrow(prob.peaks),
      n.groups=nrow(group.dt),
      most.freq.group=most1$groups,
      most.freq.prop=most1$samples.prop,
      most.next.group=most2$groups,
      most.next.prop=most2$samples.prop,
      least.freq.group=least1$groups,
      least.freq.prop=least1$samples.prop,
      least.next.group=least2$groups,
      least.next.prop=least2$samples.prop,
      exact.pvalue=exact$p.value,
      loss.diff,
      sample.counts=paste(group.vec, collapse="/"),
      separate.problem,
      joint.problem=basename(joint.prob.dir))
    joint.peaks.dt.list[[loss.tsv]] <- data.table(
      separate.problem, prob.peaks)
  }
  joint.peaks.dt <- do.call(rbind, joint.peaks.dt.list)
  input.pred <- do.call(rbind, input.pred.list)[orderChrom(chrom, peakStart),]
}

joint.peaks.dt <- jobPeaks[, {
  mean.vec <- means[[1]]
  list(
    separate.problem=problem.name,
    sample.path=names(mean.vec),
    mean=as.double(mean.vec),
    sample.id=sub(".*/", "", names(mean.vec)),
    sample.group=sub("/.*", "", names(mean.vec)),
    peakBases=peakEnd-peakStart)
}, by=list(chrom, peakStart, peakEnd, peak.name)]

group.counts <- joint.peaks.dt[, {
  tab <- sort(table(sample.group))
  nonzero <- tab[tab!=0]
  list(
    n.Input=sum(sample.group=="Input"),
    n.samples=.N,
    n.groups=length(tab),
    sample.counts=paste(paste0(
      names(nonzero), ":", nonzero), collapse=","))
  }, by=list(chrom, peakStart, peakEnd, peak.name, peakBases)]

group.counts.wide <- dcast(
  joint.peaks.dt, chrom + peakStart + peakEnd + peak.name ~ sample.group, length,
  value.var="peakBases")#to avoid message.
group.counts.mat <- as.matrix(
  group.counts.wide[, sample.group.totals$sample.group, with=FALSE])
rownames(group.counts.mat) <- group.counts.wide$peak.name

group.prop.mat <- group.counts.mat / matrix(
  sample.group.totals$samples.total,
  nrow(group.counts.mat),
  ncol(group.counts.mat),
  byrow=TRUE)

group.prop.tall <- melt(data.table(
  group.counts.wide[,list(chrom, peakStart, peakEnd)],
  group.prop.mat),
  id.vars=c("chrom", "peakStart", "peakEnd"),
  variable.name="sample.group",
  value.name="samples.prop")
setkey(group.prop.tall, chrom, peakStart, peakEnd, samples.prop)
group.prop.groups <- group.prop.tall[, list(
  groups=paste(sample.group, collapse=",")
), by=list(chrom, peakStart, peakEnd, samples.prop)]

most.least <- group.prop.groups[, data.table(
  most.freq.group=groups[.N],
  most.freq.prop=samples.prop[.N],
  most.next.group=if(.N==1)NA_character_ else groups[.N-1],
  most.next.prop=if(.N==1)NA_real_ else samples.prop[.N-1],
  least.freq.group=groups[1],
  least.freq.prop=samples.prop[1],
  least.next.group=if(.N==1)NA_character_ else groups[2],
  least.next.prop=if(.N==1)NA_real_ else samples.prop[2]
), by=list(chrom, peakStart, peakEnd)]

input.pred <- most.least[group.counts, on=list(chrom, peakStart, peakEnd)]
input.pred[, exact.pvalue := apply(
  group.counts.mat[peak.name,], 1, function(samples.with.peaks){
    count.mat <- rbind(
      samples.with.peaks,
      sample.group.totals$samples.total-samples.with.peaks)
    exact <- fisher.test(count.mat)
    exact$p.value
  })]
setkey(jobPeaks, peak.name)
input.pred[, loss.diff := jobPeaks[input.pred$peak.name, loss.diff] ]
input.pred[, separate.problem := jobPeaks[input.pred$peak.name, problem.name] ]
input.pred[, peakBases := peakEnd - peakStart]
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

## Save samples/groupID/sampleID/joint_peaks.bedGraph files.
setkey(joint.peaks.dt, sample.path)
out.path.vec <- unique(joint.peaks.dt$sample.path)
joint.peaks.dt[, mean.str := sprintf("%.2f", mean)]
for(out.path in out.path.vec){
  sample.peaks <- joint.peaks.dt[out.path]
  out.file <- file.path(set.dir, "samples", out.path, "joint_peaks.bedGraph")
  fwrite(
    sample.peaks[, .(chrom, peakStart, peakEnd, mean.str)],
    out.file,
    sep="\t",
    col.names=FALSE,
    quote=FALSE)
}

specific.html.list <- list()
for(sg in sample.group.totals$sample.group){
  specific.html.list[[sg]] <- sprintf(
    '<h3>%s</h3>', sg)
  group.peak.list <- list(
    most=input.pred[
      most.freq.group==sg &
      most.freq.prop==1 &
      most.next.prop==0,][order(-loss.diff), .(
        peak.name, loss.diff, samples=sample.counts,
        chrom, peakStart, peakEnd, peakBases)],
    least=input.pred[
      least.freq.group==sg &
      least.freq.prop==0 &
      least.next.prop==1,][order(-loss.diff), .(
        peak.name, loss.diff, samples=least.next.group,
        chrom, peakStart, peakEnd, peakBases)])
  for(most.or.least in names(group.peak.list)){
    group.peaks <- group.peak.list[[most.or.least]]
    pre.msg <- paste0(
      "<p>",
      nrow(group.peaks),
      " genomic region",
      ifelse(nrow(group.peaks)==1, " was", "s were"))
    msg <- if(most.or.least=="most"){
      paste0(
        pre.msg,
        " predicted to have a peak in each ",
        sg, " sample, and no peaks in other samples.</p>")
    }else{
      paste0(
        pre.msg,
        " predicted to have no peaks in any ",
        sg, " samples, and at least one other",
        " group with peaks in all samples.</p>")
    }
    tab <- if(nrow(group.peaks)){
      group.peaks[, peak := sprintf('
<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?position=%s:%d-%d">%s</a>
', chrom,
as.integer(peakStart-peakBases*zoom.factor),
as.integer(peakEnd+peakBases*zoom.factor),
peak.name)]
      xt <- xtable(group.peaks[, .(peak, peakBases, loss.diff, samples)])
      print(
        xt, type="html", sanitize.text.function=identity)
    }else ""
    specific.html.list[[paste(sg, most.or.least)]] <- c(msg, tab)   
  }
}
specific.html.vec <- do.call(c, specific.html.list)

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
(<a href="hub.txt">track hub</a>,
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
  problems.html,
  '<h2>Output: predicted group-specific peaks</h2>',
  '<p>These are great candidates for re-labeling.</p>',
  specific.html.vec
  )
writeLines(html.vec, file.path(set.dir, "index.html"))

## Write peaks_summary.bed
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
bed.dt <- input.pred[specificity != "non-specific",]
max.samples <- max(bed.dt$n.samples)
bed.dt[, score := as.integer((n.samples/max.samples)*1000) ]
fwrite(
  bed.dt[, .(chrom, peakStart, peakEnd, sample.counts, score)],
  peaks.bed,
  sep="\t",
  col.names=FALSE)

