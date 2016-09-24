library(data.table)
regions.RData.vec <- Sys.glob("chunks/H*/*/regions.RData")
labels.list <- list()
for(regions.RData in regions.RData.vec){
  load(regions.RData)
  chunk.dir <- dirname(regions.RData)
  chunk.id <- basename(chunk.dir)
  set.dir <- dirname(chunk.dir)
  set.name <- basename(set.dir)
  labels.list[[regions.RData]] <- data.table(set.name, chunk.id, regions)
}
labels <- do.call(rbind, labels.list)
save(labels, file="labels.RData")

labels[, chunk.name := paste0(set.name, "/", chunk.id)]
test.chunk.vec <-
  c("H3K36me3_AM_immune/11", "H3K36me3_AM_immune/12", "H3K36me3_AM_immune/15", 
    "H3K36me3_AM_immune/19", "H3K36me3_AM_immune/20", "H3K36me3_AM_immune/9")
sname <- "H3K36me3_AM_immune"
train.labels <- labels[set.name==sname & (!chunk.name %in% test.chunk.vec),]
set.dir <- file.path("labels", paste0(sname, "_folds2-4"))
train.by.sample <- split(train.labels, train.labels$sample.id, drop=TRUE)
for(sample.id in names(train.by.sample)){
  sample.dir <- file.path(set.dir, sample.id)
  dir.create(sample.dir, showWarnings=FALSE, recursive=TRUE)
  sample.labels <- train.by.sample[[sample.id]]
  write.table(
    sample.labels[, .(chrom, chromStart, chromEnd, annotation)],
    file.path(sample.dir, "labels.bed"),
    sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
  coverage.bedGraph <- file.path(sample.dir, "coverage.bedGraph")
  unlink(coverage.bedGraph)
  file.symlink(
    file.path("~/PeakSegFPOP/data", sample.id, "H3K36me3.bedGraph"),
    coverage.bedGraph)
}
