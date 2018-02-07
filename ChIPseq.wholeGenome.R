model.RData.vec <- Sys.glob("labels/*/model.RData")
data.sets.list <- list()
for(model.RData in model.RData.vec){
  objs <- load(model.RData)
  if(all(c("features", "targets") %in% objs)){
    data.sets.list[[model.RData]] <- list(
      feature.mat=features,
      target.mat=targets)
  }
}

joint.model.RData.vec <- Sys.glob("labels/*/joint.model.RData")
for(joint.model.RData in joint.model.RData.vec){
  objs <- load(joint.model.RData)
  list.of.features <- lapply(problems.list, "[[", "features")
  feature.mat <- t(rbind(
    ##sum=sapply(list.of.features, colSums),
    mean=sapply(list.of.features, colMeans)))
  stopifnot(nrow(feature.mat)==length(problems.list))
  data.sets.list[[joint.model.RData]] <- list(
    feature.mat=feature.mat,
    target.mat=t(sapply(problems.list, "[[", "target")))
}

saveRDS(data.sets.list, "ChIPseq.wholeGenome.rds")
