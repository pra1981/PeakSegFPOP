arg.vec <- c("test/input", "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-", "hg19", "email@domain.com")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 4){
  stop("usage: Rscript create_track_hub.R data_dir http://url_prefix hg19 email@domain.com")
}
pre <- "http://hgdownload.soe.ucsc.edu/goldenPath/"

system.or.stop <- function(cmd){
  cat(cmd, "\n")
  code <- system(cmd)
  if(code != 0){
    stop("non-zero exit code ", code)
  }
}
options(warn=2)

library(data.table)
data.dir <- arg.vec[1]
url.prefix <- arg.vec[2]
genome <- arg.vec[3]
email <- arg.vec[4]

## First make sure we have the chromInfo file for this genome.
chromInfo.txt <- paste0(genome, "_chromInfo.txt")
if(!file.exists(chromInfo.txt)){
  chromInfo.url <- paste0(pre, genome, "/database/chromInfo.txt.gz")
  chromInfo.gz <- paste(chromInfo.txt, ".gz")
  download.file(chromInfo.url, chromInfo.gz)
  system.or.stop(paste("zcat", gz, ">", chromInfo.txt))
}

## Then create bedGraph files if necessary.
bedGraph.file.vec <- Sys.glob(file.path(
  data.dir, "samples", "*", "*", "coverage.bedGraph"))
for(bedGraph.file in bedGraph.file.vec){
  bigWig <- sub("bedGraph$", "bigWig", bedGraph.file)
  if(!file.exists(bigWig)){
    cmd <- paste("bedGraphToBigWig", bedGraph.file, chromInfo.txt, bigWig)
    system.or.stop(cmd)
  }
}
bigWig.glob <- file.path(data.dir, "samples", "*", "*", "coverage.bigWig")
bigWig.file.vec <- Sys.glob(bigWig.glob)
if(length(bigWig.file.vec)==0){
  stop("no ", bigWig.glob, " files")
}
url.vec <- paste0(url.prefix, bigWig.file.vec)
sample.path.vec <- dirname(bigWig.file.vec)
sample.id.vec <- basename(sample.path.vec)
group.path.vec <- dirname(sample.path.vec)
group.id.vec <- basename(group.path.vec)
group.names <- unique(group.id.vec)
##dput(RColorBrewer::brewer.pal(Inf, "Set3"))
maybe.short <- c(
  "#8DD3C7",
  ##"#FFFFB3",#yellow
  "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", 
  "#B3DE69", "#FCCDE5",
  "#D9D9D9",#grey
  "#BC80BD", "#CCEBC5", "#FFED6F"
)
group.colors <- rep(maybe.short, l=length(group.names))
names(group.colors) <- group.names
data.name <- basename(data.dir)

joint_peaks.bedGraph.vec <- sub(
  "coverage.bigWig$", "joint_peaks.bedGraph", bigWig.file.vec)
for(joint_peaks.bedGraph in joint_peaks.bedGraph.vec){
  if(file.exists(joint_peaks.bedGraph)){
    joint_peaks.bigWig <- sub("bedGraph$", "bigWig", joint_peaks.bedGraph)
    cmd <- paste(
      "bedGraphToBigWig", joint_peaks.bedGraph,
      chromInfo.txt, joint_peaks.bigWig)
    system.or.stop(cmd)
  }
}

## Write genomes.txt
writeLines(paste0("
genome ", genome, "
trackDb trackDb.txt
"), file.path(data.dir, "genomes.txt"))

## Write hub.txt
writeLines(paste0("
hub ", data.name, "
shortLabel ", data.name, "
longLabel ", data.name, "
genomesFile genomes.txt
email ", email), file.path(data.dir, "hub.txt"))

## create jointProblems.bigBed
jproblems.glob <- file.path(data.dir, "problems", "*", "jointProblems.bed")
jprobs <- fread(paste("cat", jproblems.glob))
setnames(jprobs, c("chrom", "problemStart", "problemEnd"))
sizes.dt <- fread(chromInfo.txt)
names(sizes.dt)[1:2] <- c("chrom", "chromEnd")
join.dt <- jprobs[sizes.dt, on=list(chrom)]
join.dt[, problemStart := ifelse(problemStart < 0, 0, problemStart)]
join.dt[, problemEnd := ifelse(problemEnd < chromEnd, problemEnd, chromEnd)]
setkey(join.dt, chrom, problemStart, problemEnd)
write.table(
  join.dt[, .(chrom, problemStart, problemEnd)],
  file.path(data.dir, "jointProblems.bed"),
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE)

bedToBigBed <- function(bed, opt=""){
  bigBed <- sub("bed$", "bigBed", bed)
  cmd <- paste(
    "bedToBigBed",
    opt,
    bed, chromInfo.txt,
    bigBed)
  system.or.stop(cmd)
  bigBed
}
bed.num.vec <- c(
  all_labels=9,
  problems=3,
  jointProblems=3,
  peaks_summary=5)
bigBed.list <- list()
for(bed.name in names(bed.num.vec)){
  bed.file <- file.path(data.dir, paste0(bed.name, ".bed"))
  if(file.exists(bed.file)){
    bigBed.list[[bed.name]] <- bedToBigBed(bed.file)
  }
}

bed.track.vec <- paste0("
 track ", names(bigBed.list), "
 type bigBed ", bed.num.vec[names(bigBed.list)], "
 parent ", data.name, " on
 shortLabel _model_", names(bigBed.list), "
 longLabel _model_", names(bigBed.list), "
 visibility pack
 itemRgb ", ifelse(names(bigBed.list)=="all_labels", "on", "off"), "
 spectrum ", ifelse(names(bigBed.list)=="peaks_summary", "on", "off"), "
 bigDataUrl ", paste0(url.prefix, unlist(bigBed.list)))

track.id.vec <- paste0(group.id.vec, "_", sample.id.vec)
track.vec <- paste0("
 track ", track.id.vec, "
 type bigWig
 container multiWig
 shortLabel ", track.id.vec, "
 longLabel ", group.id.vec, " | ", sample.id.vec, "
 visibility full
 graphType poins
 aggregate transparentOverlay
 showSubtrackColorOnUi on
 parent ", data.name, " on
 maxHeightPixels 25:20:8
 autoScale on

  track ", track.id.vec, "Counts
  bigDataUrl ", url.vec, "
  shortLabel ", track.id.vec, "Counts
  longLabel ", group.id.vec, " | ", sample.id.vec, " | Counts
  parent ", track.id.vec, " on
  type bigWig
  color ", apply(col2rgb(group.colors[group.id.vec]), 2, paste, collapse=","), "

  track ", track.id.vec, "Peaks
  bigDataUrl ", sub("coverage.bigWig$", "joint_peaks.bigWig", url.vec), "
  shortLabel ", track.id.vec, "Peaks
  longLabel ", group.id.vec, " | ", sample.id.vec, " | Peaks
  parent ", track.id.vec, " on
  type bigWig
  color 0,0,0
")

u.group.vec <- unique(group.id.vec)
equals.vec <- paste0(u.group.vec, "=", u.group.vec)
track.content <- paste0("
track ", data.name, "
superTrack on show
shortLabel ", data.name, "
longLabel ", data.name, "

", paste(bed.track.vec, collapse="\n"), "

", paste(track.vec, collapse="\n"), "
")

writeLines(track.content, file.path(data.dir, "trackDb.txt"))
