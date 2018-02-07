library(data.table)
models.dt <- fread("target.intervals.models.csv")
iterations.dt <- fread("target.intervals.iterations.csv")

errors.dt <- models.dt[, list(
  min.errors=min(fp+fn)
  ), by=list(prob.dir)]
errors.dt[order(min.errors)]
setkey(errors.dt, prob.dir)
(err.tab <- table(errors.dt$min.errors))
err.tab/sum(err.tab)
## most errors 9 -- there are a lot of labels 14+23 possible errors.
models.dt[prob.dir=="H3K4me3_TDH_immune/samples/tcell/McGill0029/problems/chr3:93504854-194041961"]

## Alternate method for determining trivial intervals, look at last iteration.
last.iteration <- iterations.dt[, .SD[iteration==max(iteration)], by=prob.dir]
last.iteration[min.log.lambda==-Inf & max.log.lambda==Inf, table(iteration)]
last.not.trivial <- last.iteration[!(min.log.lambda==-Inf & max.log.lambda==Inf)]
setnames(last.not.trivial, c("prob.dir", "last.iteration", "last.min", "last.max"))
last.not.trivial[, table(last.iteration)]
last.not.trivial[last.iteration==5]
setkey(iterations.dt, prob.dir)
setkey(last.not.trivial, prob.dir)
not.trivial <- errors.dt[iterations.dt][last.not.trivial]

## time for each iteration.
seconds.dt <- models.dt[, list(
  max.seconds=max(seconds)
  ), by=list(prob.dir, iteration)]
setkey(seconds.dt, prob.dir, iteration)
setkey(not.trivial, prob.dir, iteration)
iteration.seconds <- seconds.dt[not.trivial]

f <- function(x){
  d <- diff(x)
  c(Inf, ifelse(is.finite(d), abs(d), Inf))
}
iteration.seconds[, `:=`(
  cum.seconds=cumsum(max.seconds),
  diff.min=f(min.log.lambda),
  diff.max=f(max.log.lambda)
  ), by=prob.dir]
iteration.seconds[, last.diff := abs(min.log.lambda-last.min)+abs(max.log.lambda-last.max)]
iteration.seconds[, cum.hours := cum.seconds / 60 / 60]
one.hour <- iteration.seconds[cum.hours < 1, .SD[which.max(iteration)], by=prob.dir]
table(one.hour$iteration)
iteration.seconds[iteration==1][order(cum.hours)]
models.dt[prob.dir=="/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/CTCF_TDH_ENCODE/samples/Input/MS026601/problems/chr4:75452279-190000000"]

iteration.seconds[prob.dir=="/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/H3K27ac_TDH_some/samples/kidney/MS040302/problems/chr1:30028082-91800000"]
models.dt[prob.dir=="/gs/project/mugqic/projects/thocking/blueprint/PeakSegFPOP-labels/H3K27ac_TDH_some/samples/kidney/MS040302/problems/chr1:30028082-91800000"]

## common point in models with many iterations is that they do not obtain 0 errors.
data.frame(models.dt[prob.dir %in% last.not.trivial[20 < last.iteration, prob.dir], data.table(
  iteration, fp, fn, peaks, penalty, bedGraph.lines)])

library(ggplot2)
ggplot()+
  geom_line(aes(
    iteration, log10(diff.min+diff.max), group=prob.dir),
            data=iteration.seconds)

ggplot()+
  theme_bw()+
  scale_color_gradient(low="grey", high="black")+
  geom_line(aes(
    iteration, log10(last.diff), group=prob.dir, color=min.errors),
            data=iteration.seconds)

ggplot()+
  theme_bw()+
  geom_line(aes(
    log10(cum.hours), log10(last.diff), group=prob.dir),
            data=iteration.seconds)

ggplot()+
  theme_bw()+
  geom_line(aes(
    log10(cum.hours), log10(last.diff), group=prob.dir),
            data=iteration.seconds[cum.hours<1])

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_wrap("min.errors")+
  geom_line(aes(
    iteration, log10(last.diff), group=prob.dir),
            data=iteration.seconds)
