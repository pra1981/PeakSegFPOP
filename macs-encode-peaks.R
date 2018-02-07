#!/home/thocking/bin/Rscript
#PBS -l nodes=1:ppn=4
#PBS -l walltime=24:00:00
#PBS -A bws-221-ae
#PBS -m ae
#PBS -M tdhock5@gmail.com
#PBS -V
#PBS -o /home/thocking/PeakSegFPOP/DFilter-all.out
#PBS -e /home/thocking/PeakSegFPOP/DFilter-all.err
#PBS -N DF-all
alignments.bam.vec <- grep("Input", Sys.glob("~/PeakSegFPOP/labels/H*_ENCODE/samples/*/*/alignments.bam"), invert=TRUE, value=TRUE)

for(bam.i in seq_along(alignments.bam.vec)){
  alignments.bam <- alignments.bam.vec[[bam.i]]
  cat(sprintf("%4d / %4d %s\n", bam.i, length(alignments.bam.vec), alignments.bam))
  sh <- ifelse(
    grepl("H3K36me3", alignments.bam),
    "macs2-broad-no-control-default.sh",
    "macs2-new-no-control-default.sh")
  cmd <- paste(sh, alignments.bam, dirname(alignments.bam))
  system(cmd)
}
