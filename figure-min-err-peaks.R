load("min.err.peaks.RData")
library(data.table)

wide <- dcast(
  min.err.peaks,
  experiment + sample + group + problem ~ interval.type,
  value.var=c(
    "min.peaks", "max.peaks", "n.infeasible",
    "min.incorrect.labels", "cpu.hours"))

ggplot()+
  geom_point(aes(
    min.incorrect.labels_problems,
    min.incorrect.labels_problems_infeasibleInf),
    data=wide)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

ggplot()+
  geom_point(aes(
    n.infeasible_problems,
    n.infeasible_problems_infeasibleInf),
    data=wide)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

only.new <- min.err.peaks[interval.type=="problems",]

hour.refs <- data.table(
  hjust=c(1, 1, 0),
  bases=c(1e8, 1e8, 1e7),
  hours=c(10, 60, 60*10)/60,
  label=c("10 minutes", "1 hour", "10 hours"))

ggplot()+
  geom_point(aes(cpu.hours, walltime.hours), data=only.new)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

ggplot()+
  geom_point(aes(log10(cpu.hours), log10(walltime.hours)),
             shape=1,
             data=only.new)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

ggplot()+
  ##ggtitle("target interval computation time")+
  scale_x_log10("segmentation problem size, mega bases (Mb)")+
  scale_y_log10("hours to compute target interval of penalty values with minimum label error")+
  geom_hline(aes(yintercept=hours), data=hour.refs, color="grey")+
  geom_text(aes(bases/1e6, hours, label=label, hjust=hjust),
            vjust=1.5,
            data=hour.refs, color="grey")+
  geom_point(aes(bases/1e6, walltime.hours, color=experiment),
             shape=1,
             data=only.new)
ggsave("~/projects/PeakSegFPOP-paper/figure-target-interval-time.pdf")

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample ~ experiment)+
  scale_color_gradient(low="black", high="red")+
  geom_segment(aes(bases, log10(min.peaks),
                   color=min.incorrect.labels,
                   xend=bases, yend=log10(max.peaks)),
               size=1,
               data=only.new)

only.new[, mid.peaks := (min.peaks+max.peaks)/2]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ experiment)+
  scale_fill_gradient(low="white", high="red")+
  geom_point(aes(log10(bases), log10(mid.peaks),
                 fill=min.incorrect.labels),
             shape=21,
             data=only.new)

incorrect.labels <- only.new[, list(
  sum=sum(min.incorrect.labels),
  mean=mean(min.incorrect.labels)
  ), by=.(experiment, problem, bases)]
setkey(incorrect.labels, bases)
prob.levs <- unique(paste(incorrect.labels$problem))

only.new[, problem.fac := factor(problem, prob.levs)]
only.new[, sample.id := sub("McGill0", "", sample)]
only.new[, status := ifelse({
  n.infeasible==0
}, "all.feasible", ifelse(
  n.feasible==0, "all.infeasible", "both"))]
ggplot()+
  theme_grey()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ ., scales="free", space="free")+
  scale_fill_gradient(low="white", high="red")+
  scale_linetype_manual(values=c(all.feasible=0, all.infeasible=1, both=2))+
  geom_tile(aes(sample.id, problem.fac,
                fill=min.incorrect.labels, linetype=status),
            color="black",
            size=1,
            data=only.new)

only.new[, problem.sample := paste(problem, sample)]
hg19.problems <- fread("hg19_problems.bed")
setnames(hg19.problems, c("chrom", "problemStart","problemEnd"))
hg19.problems[, problem := sprintf("%s:%d-%d", chrom, problemStart, problemEnd)]
hg19.problems[, chrom.fac := factorChrom(chrom)]
setkey(hg19.problems, problem)
setkey(only.new, problem)
labeled.problems <- hg19.problems[only.new]
stopifnot(nrow(labeled.problems)==nrow(only.new))
other.probs <- data.table(labeled.problems)
setkey(labeled.problems, sample)
setkey(other.probs, sample)
big.problems <- labeled.problems[other.probs, allow.cartesian=TRUE]
big.problems[, experiment := i.experiment]
only.new[, tooltip := paste0(
  format(bases, big.mark=",", scientific=FALSE), " bases, ",
  min.peaks,
  ifelse(min.peaks==max.peaks, "", paste0("-", max.peaks)),
  " peaks achieve ",
  min.incorrect.labels, " error",
  ifelse(min.incorrect.labels==1, "", "s"),
  " for ", problem.sample)]
library(animint)
break.vec <- rev(unique(as.integer(seq(
  0, max(only.new$min.incorrect.labels), l=5))))
viz <- list(
  title="PeakSegFPOP target intervals versus problem size", 
  heatmap=ggplot()+
    ggtitle("Label error heatmap, select sample and segmentation problem")+
    ylab("segmentation problem")+
    theme_grey()+
    theme_animint(height=350, width=900)+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(experiment ~ ., scales="free", space="free")+
    scale_fill_gradient(
      "incorrect labels",
      low="white", high="red",
      breaks=break.vec
    )+
    geom_tile(aes(sample.id, problem.fac, fill=min.incorrect.labels,
                  tooltip=tooltip,
                  showSelected=experiment,
                  clickSelects=problem.sample),
              data=only.new),
  genome=ggplot()+
    ggtitle("Select segmentation problem in hg19")+
    xlab("position on chromosome (mega bases)")+
    ylab("chromosome")+
    guides(color="none")+
    theme_animint(height=450)+
    geom_segment(aes(problemStart/1e6, chrom.fac,
                     xend=problemEnd/1e6, yend=chrom.fac),
                 color="grey",
                 data=hg19.problems)+
    geom_point(aes(problemStart/1e6, chrom.fac),
               color="grey",
               data=hg19.problems)+
    geom_segment(aes(i.problemStart/1e6, i.chrom.fac,
                     xend=i.problemEnd/1e6, yend=i.chrom.fac,
                     showSelected=problem.sample,
                     showSelected2=experiment,
                     clickSelects.variable="problem.sample",
                     clickSelects.value=i.problem.sample,
                     key=paste(i.problem.sample,i.experiment),
                     color=i.experiment),
                 ## validate_params=FALSE,
                 ## chunk_vars="problem.sample",
                 size=9,
                 alpha=0.7,
                 data=big.problems),
  ## scatter=ggplot()+
  ##   theme_bw()+
  ##   theme(panel.margin=grid::unit(0, "lines"))+
  ##   theme_animint(width=600)+
  ##   facet_grid(. ~ experiment)+
  ##   scale_fill_gradient(low="white", high="red")+
  ##   guides(fill="none")+
  ##   geom_segment(aes(
  ##     log10(min.peaks), log10(bases),
  ##     xend=log10(max.peaks), yend=log10(bases),
  ##     tooltip=tooltip,
  ##     showSelected=problem.sample),
  ##                size=4,
  ##                data=only.new)+
  ##   geom_point(aes(
  ##     log10(mid.peaks), log10(bases), 
  ##     tooltip=tooltip,
  ##     clickSelects=problem.sample,
  ##     fill=min.incorrect.labels),
  ##              shape=21,
  ##              alpha=0.6,
  ##              size=4,
  ##              data=only.new),
  compare=ggplot()+
    ggtitle("Problem size vs optimal peaks")+
    xlab("log10(peaks with minimum incorrect labels)")+
    ylab("log10(bases in segmentation problem)")+
    theme_animint(height=450)+
    geom_segment(aes(
      log10(min.peaks), log10(bases),
      tooltip=tooltip,
      xend=log10(max.peaks), yend=log10(bases),
      showSelected2=experiment,
      showSelected=problem.sample),
                 size=4,
                 data=only.new)+
    geom_point(aes(
      log10(mid.peaks), log10(bases), 
      tooltip=tooltip,
      clickSelects=problem.sample,
      color=experiment),
               alpha=0.6,
               size=4,
      data=only.new))
animint2dir(viz, "figure-min-err-peaks")

same.prob <- "chr11:96437584-134946516"
same.wide <- dcast(
  only.new[problem==same.prob,],
  sample ~ experiment,
  value.var=c("min.peaks", "mid.peaks", "max.peaks"))
ggplot()+
  ggtitle(paste(
    "Peaks with minimum incorrect labels in",
    same.prob))+
  theme_bw()+
  geom_abline(slope=1, intercept=0, color="grey")+
  geom_rect(aes(
    xmin=min.peaks_H3K4me3, xmax=max.peaks_H3K4me3,
    ymin=min.peaks_H3K36me3, ymax=max.peaks_H3K36me3),
    color="black",
    fill="grey",
    alpha=0.2,
    data=same.wide)+
  coord_equal()+
  scale_x_log10(
    "H3K4me3 peaks",
    limits=c(10, 1000),
    breaks=c(1, 10, 100, 1000))+
  scale_y_log10(
    "H3K36me3 peaks",
    limits=c(1, 100),
    breaks=c(1, 10, 100, 1000))
ggsave("~/projects/PeakSegFPOP-paper/figure-min-err-peaks-compare.pdf")
