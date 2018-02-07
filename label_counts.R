library(data.table)
data.table(set.dir=Sys.glob("labels/*"))[, {
  samples.problems <- data.table(
    problem.dir=Sys.glob(file.path(set.dir, "samples/*/*/problems/*")))
  samples.problems[, labels.bed := file.exists(file.path(problem.dir, "labels.bed"))]
  samples.problems[labels.bed==TRUE, labels := {
    fread(paste("wc -l", paste(file.path(problem.dir, "labels.bed"), collapse=" ")))[-.N, V1]
  }]
  samples.problems[labels.bed==FALSE, labels := 0]
  samples.problems[, sample.dir := dirname(dirname(problem.dir))]
  samples.problems[, sample.id := basename(sample.dir)]
  samples.problems[, group.id := basename(dirname(sample.dir))]
  samples.problems[, problem.name := basename(problem.dir)]
  samples.problems[, problem.fac := factorChrom(problem.name)]

  ggplot()+
    theme(
      panel.margin=grid::unit(0, "lines"),
      axis.text.x=element_text(angle=90, hjust=1))+
    geom_tile(aes(
      problem.fac, sample.id, fill=labels),
      data=samples.problems)+
    coord_equal()+
    ##facet_grid(group.id ~ ., scales="free", space="free")+
    scale_fill_gradient(low="white", high="red")
  
}, by=list(set.dir)]
