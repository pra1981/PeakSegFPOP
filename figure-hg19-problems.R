file.name.dt <- data.table(file.name=c("hg19_problems.bed", "hg19_original.bed"))
problems <- file.name.dt[, {
  dt <- fread(file.name, colClasses=list(integer="V2"))
  str(dt)
  dt
}, by=file.name]
setnames(problems, c("file.name", "chrom", "problemStart", "problemEnd"))

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chrom ~ .)+
  geom_segment(aes(problemStart/1e6, file.name,
                   xend=problemEnd/1e6, yend=file.name),
               data=problems)+
  geom_point(aes(problemStart/1e6, file.name),
             shape=1,
               data=problems)
ggsave("figure-hg19-problems.png")
