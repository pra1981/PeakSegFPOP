load("min.err.peaks.RData")
library(data.table)

ggplot()+
  geom_point(aes(cpu.hours, walltime.hours), data=min.err.peaks)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

ggplot()+
  geom_point(aes(log10(cpu.hours), log10(walltime.hours)), data=min.err.peaks)+
  coord_equal()+
  geom_abline(slope=1, intercept=0, color="grey")

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(sample ~ experiment)+
  scale_color_gradient(low="black", high="red")+
  geom_segment(aes(bases, log10(min.peaks),
                   color=min.incorrect.labels,
                   xend=bases, yend=log10(max.peaks)),
               size=1,
               data=min.err.peaks)

min.err.peaks[, mid.peaks := (min.peaks+max.peaks)/2]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ experiment)+
  scale_fill_gradient(low="white", high="red")+
  geom_point(aes(log10(bases), log10(mid.peaks),
                 fill=min.incorrect.labels),
             shape=21,
             data=min.err.peaks)

incorrect.labels <- min.err.peaks[, list(
  sum=sum(min.incorrect.labels),
  mean=mean(min.incorrect.labels)
  ), by=.(experiment, problem, bases)]
setkey(incorrect.labels, bases)
prob.levs <- unique(paste(incorrect.labels$problem))

min.err.peaks[, problem.fac := factor(problem, prob.levs)]
min.err.peaks[, sample.id := sub("McGill0", "", sample)]
min.err.peaks[, status := ifelse({
  n.infeasible==0
}, "all.feasible", ifelse(
  n.feasible==0, "all.infeasible", "both"))]
ggplot()+
  theme_grey()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ ., scales="free", space="free")+
  scale_fill_gradient(low="white", high="red")+
  scale_linetype_manual(values=c(all.feasible=0, all.infeasible=1, both=2))+
  geom_tile(aes(sample.id, problem.fac, fill=min.incorrect.labels, linetype=status),
            color="black",
            size=1,
            data=min.err.peaks)

min.err.peaks[, problem.sample := paste(problem, sample)]
min.err.peaks[, tooltip := paste0(
  format(bases, big.mark=",", scientific=FALSE), " bases, ",
  min.peaks,
  ifelse(min.peaks==max.peaks, "", paste0("-", max.peaks)),
  " peaks achieve ",
  min.incorrect.labels, " error",
  ifelse(min.incorrect.labels==1, "", "s"),
  " for ", problem.sample)]
viz <- list(
  title="PeakSegFPOP target intervals versus problem size", 
  heatmap=ggplot()+
    theme_grey()+
    theme_animint(width=800)+
    theme(panel.margin=grid::unit(0, "lines"))+
    facet_grid(experiment ~ ., scales="free", space="free")+
    scale_fill_gradient(low="white", high="red")+
    geom_tile(aes(sample.id, problem.fac, fill=min.incorrect.labels,
                  tooltip=tooltip,
                  clickSelects=sample),
              data=min.err.peaks),
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
  ##                data=min.err.peaks)+
  ##   geom_point(aes(
  ##     log10(mid.peaks), log10(bases), 
  ##     tooltip=tooltip,
  ##     clickSelects=problem.sample,
  ##     fill=min.incorrect.labels),
  ##              shape=21,
  ##              alpha=0.6,
  ##              size=4,
  ##              data=min.err.peaks),
  compare=ggplot()+
    geom_segment(aes(
      log10(min.peaks), log10(bases),
      tooltip=tooltip,
      xend=log10(max.peaks), yend=log10(bases),
      showSelected=sample),
                 size=4,
                 data=min.err.peaks)+
    geom_point(aes(
      log10(mid.peaks), log10(bases), 
      tooltip=tooltip,
      clickSelects=sample,
      color=experiment),
               alpha=0.6,
               size=4,
               data=min.err.peaks))
library(animint)
animint2dir(viz, "figure-min-err-peaks")

