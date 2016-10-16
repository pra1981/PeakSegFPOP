## $`labels/H3K4me3_PGP_immune/samples/tcell/McGill0107/problems/chr1:144810724-145833118`
##     n.data      penalty peaks     status fn fp errors
##  1: 110632      0.00000 49268 infeasible  0  3    Inf
##  2: 110632     69.17564  2592 infeasible  0  3    Inf
##  3: 110632   1247.15929   214 infeasible  0  2    Inf
##  4: 110632   2966.96827    87 infeasible  0  2    Inf
##  5: 110632   5298.72577    51   feasible  0  1      1
##  6: 110632   7441.95143    38 infeasible  0  1    Inf
##  7: 110632   8379.67743    34 infeasible  0  1    Inf
##  8: 110632   8565.77667    33 infeasible  1  1    Inf
##  9: 110632   8700.87092    32   feasible  1  0      1
## 10: 110632   8945.30898    30   feasible  1  0      1
## 11: 110632  11980.31315    26   feasible  1  0      1
## 12: 110632  77153.72998     5   feasible  1  0      1
## 13: 110632 284686.29266     1   feasible  1  0      1
## 14: 110632          Inf     0   feasible  3  0      3

problem.dir <- "weird-small"
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
coverage <- fread(coverage.bedGraph)
setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))

c("0", "69.1756437437981", "1247.15928799506", "2966.96827273658", 
"5298.72576637565", "7441.95142536979", "8379.67742862483", "8565.77666862542", 
"8700.87092360435", "8945.30898359216", "11980.3131528345", "77153.729978158", 
  "284686.292663906", "Inf")

penalty.vec <- c("2966.96827273658", "5298.72576637565", "7441.95142536979",
                 "8565.77666862542", 
"8700.87092360435")
seg.list <- list()
change.list <- list()
for(penalty in penalty.vec){
  segments.bed <- paste0(coverage.bedGraph, "_penalty=", penalty, "_segments.bed")
  penalty.segs <- fread(segments.bed)
  setnames(penalty.segs, c("chrom", "chromStart", "chromEnd", "status", "mean"))
  penalty.segs[.N, chromStart := coverage$chromStart[1]]
  penalty.change <- penalty.segs[, data.table(
    position=chromEnd[-1],
    diff.mean=diff(mean)
    )]
  penalty.num <- as.numeric(penalty)
  change.list[[penalty]] <- data.table(
    penalty.str=penalty,
    penalty.num,
    penalty.change)
  seg.list[[penalty]] <- data.table(
    penalty.str=penalty,
    penalty.num,
    penalty.segs)
}
seg <- do.call(rbind, seg.list)
change <- do.call(rbind, change.list)

ggplot()+
  ggtitle("Dashed vertical lines between equal segment means")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(penalty.num ~ .)+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=coverage)+
  geom_vline(aes(xintercept=position/1e3),
             linetype="dashed",
             size=1,
             data=change[diff.mean==0,])+
  scale_color_manual(
    values=c(
      background="black",
      peak="red"),
    limits=c(
      "peak",
      "background"))+
  xlab("position on chr1 (kb = kilo bases)")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean,
                   color=status),
               size=1,
               data=seg)

ggsave("figure-weird-small.png")

## $`labels/H3K36me3_TDH_other/samples/skeletalMuscleCtrl/McGill0036/problems/chr21:43005559-44632664`
##     n.data      penalty peaks     status fn fp errors
##  1: 116430      0.00000 51725 infeasible  0  6    Inf
##  2: 116430     98.72452  2933 infeasible  0  6    Inf
##  3: 116430   1653.65075   248 infeasible  0  6    Inf
##  4: 116430   3761.39828    89   feasible  0  3      3
##  5: 116430   7083.18830    49   feasible  0  1      1
##  6: 116430  10064.88208    33 infeasible  0  0    Inf
##  7: 116430  12325.49606    26 infeasible  0  0    Inf
##  8: 116430  13639.03836    25 infeasible  0  0    Inf
##  9: 116430  14094.37843    23   feasible  0  0      0
## 10: 116430  14109.75705    23   feasible  0  0      0
## 11: 116430  15298.22225    21   feasible  0  0      0
## 12: 116430 140005.79564     7   feasible  0  0      0
## 13: 116430 147871.34164     6   feasible  1  1      2
## 14: 116430 200427.38916     5   feasible  1  1      2
## 15: 116430 346614.13322     3   feasible  1  1      2
## 16: 116430          Inf     0   feasible  3  0      3

problem.dir <- "weird-small2"
coverage.bedGraph <- file.path(problem.dir, "coverage.bedGraph")
coverage <- fread(coverage.bedGraph)
setnames(coverage, c("chrom", "chromStart", "chromEnd", "count"))

c("0", "98.7245208358598", "1653.65075478066", "3761.39827561476", 
"7083.18829536253", "10064.882075448", "12325.4960601653", "13639.0383603886", 
"14094.3784297225", "14109.7570451573", "15298.2222457807", "140005.795637574", 
"147871.341638583", "200427.389160114", "346614.133222383", "Inf"
)
penalty.vec <- c("1653.65075478066", "3761.39827561476", 
                 "7083.18829536253", "10064.882075448",
                 "13639.0383603886", 
"14094.3784297225", "140005.795637574")
seg.list <- list()
change.list <- list()
for(penalty in penalty.vec){
  segments.bed <- paste0(coverage.bedGraph, "_penalty=", penalty, "_segments.bed")
  penalty.segs <- fread(segments.bed)
  setnames(penalty.segs, c("chrom", "chromStart", "chromEnd", "status", "mean"))
  penalty.segs[.N, chromStart := coverage$chromStart[1]]
  penalty.change <- penalty.segs[, data.table(
    position=chromEnd[-1],
    diff.mean=diff(mean)
    )]
  penalty.num <- as.numeric(penalty)
  change.list[[penalty]] <- data.table(
    penalty.str=penalty,
    penalty.num,
    penalty.change)
  seg.list[[penalty]] <- data.table(
    penalty.str=penalty,
    penalty.num,
    penalty.segs)
}
seg <- do.call(rbind, seg.list)
change <- do.call(rbind, change.list)

ggplot()+
  ggtitle("Dashed vertical lines between equal segment means")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(penalty.num ~ .)+
  geom_step(aes(chromStart/1e3, count),
            color="grey50",
            data=coverage)+
  geom_vline(aes(xintercept=position/1e3),
             linetype="dashed",
             size=1,
             data=change[diff.mean==0,])+
  scale_color_manual(
    values=c(
      background="black",
      peak="red"),
    limits=c(
      "peak",
      "background"))+
  xlab("position on chr21 (kb = kilo bases)")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean,
                   color=status),
               size=1,
               data=seg)
ggsave("figure-weird-small2.png")
