data.dt <- data.table(data=Sys.glob("new_models/*.RData"))
data.dt[, experiment := sub("_.*", "", basename(data))]
correct.peaks <- data.dt[, {
  load(data)
  correct.peaks
}, by=experiment]
correct.peaks[, sample := sub("/.*", "", sub(".*McGill0", "", problem.dir))]
correct.peaks[, bases := chromEnd-chromStart]
correct.peaks[, log10.bases := log10(bases)]
model <- correct.peaks[, list(
  mean=mean(log10.bases),
  sd=sd(log10.bases)
), by=experiment]
times <- 3
model[, upper.lim := mean + times*sd]
model[, lower.lim := mean - times*sd]
model[, upper.bases := 10^(upper.lim)]
model[, lower.bases := 10^(lower.lim)]
setkey(model, experiment)
setkey(correct.peaks, experiment)
model.peaks <- correct.peaks[model]
model.peaks[, {
  is.in <- lower.bases < bases & bases < upper.bases
  list(percent.in=sum(is.in)/.N)
}, by=experiment]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_histogram(aes(log10.bases, ..density..),
                 data=correct.peaks)+
  geom_vline(aes(xintercept=mean), data=model)+
  geom_vline(aes(xintercept=upper.lim), data=model)+
  geom_vline(aes(xintercept=lower.lim), data=model)   

setkey(model, experiment)
norm.density <- correct.peaks[, {
  log10.bases <- seq(min(log10.bases), max(log10.bases), l=100)
  m <- model[experiment]
  data.table(
    log10.bases,
    density=dnorm(log10.bases, m$mean, m$sd))
}, by=experiment]
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ .)+
  geom_histogram(aes(log10.bases, ..density..),
                 data=correct.peaks)+
  geom_vline(aes(xintercept=mean),
             color="red",
             size=2,
             data=model)+
  geom_line(aes(log10.bases, density), color="red", size=2, data=norm.density)+
  geom_tallrect(aes(xmin=lower.lim, xmax=upper.lim),
                fill="red",
                alpha=0.25,
                color=NA,
                data=model)+
  geom_text(aes(
    mean, 1,
    label=paste0(
      "mean=",
      format(as.integer(10^mean), scientific=FALSE, big.mark=','))),
            hjust=0,
            vjust=1,
            data=model)+
  geom_text(aes(
    upper.lim, 1,
    label=paste0(
      "limit=",
      format(as.integer(upper.bases), scientific=FALSE, big.mark=','))),
            hjust=0,
            vjust=1,
    data=model)
ggsave("~/projects/PeakSegFPOP-paper/figure-peak-size-model.pdf")

model.peaks[bases > upper.bases, table(sample, experiment)]

library(penaltyLearning)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(experiment ~ sample)+
  geom_tallrect(aes(xmin=lower.lim, xmax=upper.lim),
                alpha=0.5, data=model)+
  geom_histogram(aes(log10.bases, ..density..),
                 color="blue",
                 data=correct.peaks)+
  geom_vline(aes(xintercept=mean), data=model)

