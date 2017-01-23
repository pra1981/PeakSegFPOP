arg.vec <- c("test/input", "http://hubs.hpc.mcgill.ca/~thocking/PeakSegFPOP-", "hg19", "email@domain.com")

arg.vec <- commandArgs(trailingOnly=TRUE)

if(length(arg.vec) != 4){
  stop("usage: Rscript create_track_hub.R data_dir http://url_prefix hg19 email@domain.com")
}

library(data.table)
data.dir <- arg.vec[1]
url.prefix <- arg.vec[2]
genome <- arg.vec[3]
email <- arg.vec[4]
bigWig.file.vec <- Sys.glob(file.path(data.dir, "samples", "*", "*", "coverage.bigWig"))
url.vec <- paste0(url.prefix, bigWig.file.vec)
sample.path.vec <- dirname(bigWig.file.vec)
sample.id.vec <- basename(sample.path.vec)
group.path.vec <- dirname(sample.path.vec)
group.id.vec <- basename(group.path.vec)
group.names <- unique(group.id.vec)
maybe.short <- RColorBrewer::brewer.pal(length(group.names), "Set3")
group.colors <- rep(maybe.short, l=length(group.names))
names(group.colors) <- group.names
data.name <- basename(data.dir)

writeLines("
genome hg19
trackDb trackDb.txt
", file.path(data.dir, "genomes.txt"))
writeLines(paste0("
hub ", data.name, "
shortLabel ", data.name, "
longLabel ", data.name, "
genomesFile genomes.txt
email ", email), file.path(data.dir, "hub.txt"))

track.id.vec <- paste0(group.id.vec, "_", sample.id.vec)
track.vec <- paste0("
    track ", track.id.vec, "
    type bigWig
    parent ", data.name, " on
    shortLabel ", track.id.vec, "
    longLabel ", group.id.vec, " | ", sample.id.vec, "
    bigDataUrl ", url.vec, "
    maxHeightPixels 25:25:8
    color ", apply(col2rgb(group.colors[group.id.vec]), 2, paste, collapse=","), "
    autoScale on")


track.content <- paste0(
  "track ", data.name, "
compositeTrack on
shortLabel ", data.name, "
longLabel ", data.name, "
dragAndDrop subTracks
priority 1
type bed 5
visibility pack
", paste(track.vec, collapse="\n"))

writeLines(track.content, file.path(data.dir, "trackDb.txt"))
