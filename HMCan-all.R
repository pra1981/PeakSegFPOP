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
alignments.bam.vec <- Sys.glob("~/PeakSegFPOP/labels/H*/samples/*/*/alignments.bam")
for(bam.i in seq_along(alignments.bam.vec)){
  alignments.bam <- alignments.bam.vec[[bam.i]]
  cat(sprintf("%4d / %4d %s\n", bam.i, length(alignments.bam.vec), alignments.bam))
  alignments.bed <- sub("bam$", "bed", alignments.bam)
  if(!file.exists(alignments.bed)){
    cmd <- paste("~/bin/bam2bed.sh", alignments.bam, alignments.bed)
    print(cmd)
    system(cmd)
  }
}


alignments.bed.vec <- Sys.glob("~/PeakSegFPOP/labels/H*/samples/*/*/alignments.bed")
for(bed.i in seq_along(alignments.bed.vec)){
  alignments.bed <- alignments.bed.vec[[bed.i]]
  cat(sprintf("%4d / %4d %s\n", bed.i, length(alignments.bed.vec), alignments.bed))
  sample.dir <- dirname(alignments.bed)
  group.dir <- dirname(sample.dir)
  samples.dir <- dirname(group.dir)
  set.dir <- dirname(samples.dir)
  set.name <- basename(set.dir)
  experiment <- sub("_.*", "", set.name)
  cmd <- paste0("~/bin/DFilter-", experiment, ".sh ", sample.dir)
  print(cmd)
  system(cmd)
}

