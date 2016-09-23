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
