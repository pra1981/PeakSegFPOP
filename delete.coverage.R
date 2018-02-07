#!/home/thocking/bin/Rscript
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V
#PBS -o /home/thocking/PeakSegFPOP/delete.coverage.out
#PBS -e /home/thocking/PeakSegFPOP/delete.coverage.err
#PBS -N DelCov
sample.dir.vec <- Sys.glob("~/PeakSegFPOP/labels/*/samples/*/*")
for(i in seq_along(sample.dir.vec)){
  sample.dir <- sample.dir.vec[[i]]
  cat(sprintf("%4d / %4d %s\n", i, length(sample.dir.vec), sample.dir))
  unlink(file.path(sample.dir, "problems", "*", "coverage.bedGraph"))
}
