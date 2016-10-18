argv <- "test/PeakSegJoint-one"
argv <- "test/PeakSegJoint-overlapping"
argv <- "test/PeakSegJoint-two"

argv <- commandArgs(trailingOnly=TRUE)

library(data.table)
library(PeakSegJoint)

print(argv)
if(length(argv) == 0){
  stop("usage: convert_labels.R project_dir where there is
project_dir/samples/* directories and
project_dir/labels/* files")
}
project.dir <- normalizePath(argv, mustWork=TRUE)
labels.file.vec <- Sys.glob(file.path(project.dir, "labels", "*"))

g.pos.pattern <-
  paste0("(?<chrom>chr.+?)",
         ":",
         "(?<chromStart>[0-9 ,]+)",
         "-",
         "(?<chromEnd>[0-9 ,]+)",
         " ",
         "(?<annotation>[a-zA-Z]+)",
         "(?<sample_groups>.*)")
label.colors <- 
  c(noPeaks="246,244,191",
    peakStart="255,175,175",
    peakEnd="255,76,76",
    peaks="164,69,238")


str_match_named <- function
### Parse the first occurance of pattern from each of several subject
### strings using a named capture regular expression.
(subject.vec,
### character vector of subjects.
 pattern,
### named capture regular expression (character vector of length 1).
 type.list=NULL
### named list of functions to apply to captured groups.
 ){
  stopifnot(is.character(subject.vec))
  stopifnot(0 < length(subject.vec))
  stopifnot(is.character(pattern))
  stopifnot(length(pattern)==1)
  vec.with.attrs <- regexpr(pattern, subject.vec, perl=TRUE)
  no.match <- vec.with.attrs == -1 | is.na(subject.vec)
  capture.names <- names_or_error(vec.with.attrs)
  first <- attr(vec.with.attrs, "capture.start")
  first[no.match] <- NA
  last <- attr(vec.with.attrs, "capture.length")-1+first
  last[no.match] <- NA
  subs <- substring(subject.vec, first, last)
  m <- matrix(subs, length(subject.vec), length(capture.names),
              dimnames=list(names(subject.vec), capture.names))
  apply_type_funs(m, type.list)
### A data.frame with one row for each subject and one column for each
### capture group if type.list is a list of functions. Otherwise a
### character matrix. If subject.vec has names then they will be used
### for the rownames of the returned data.frame or character
### matrix. Otherwise if there is a group named "name" then it will
### not be returned as a column, and will instead be used for the
### rownames.
}

apply_type_funs <- function
### Convert columns of match.mat using corresponding functions from
### type.list.
(match.mat,
### character matrix (matches X groups).
 type.list
### named list of functions to apply to captured groups.
 ){
  stopifnot(is.character(match.mat))
  stopifnot(is.matrix(match.mat))
  if(is.null(rownames(match.mat)) && "name" %in% colnames(match.mat)){
    rownames(match.mat) <- match.mat[, "name"]
    match.mat <- match.mat[, colnames(match.mat) != "name", drop=FALSE]
  }
  if(is.list(type.list)){
    df <- data.frame(match.mat, stringsAsFactors=FALSE)
    for(col.name in names(type.list)){
      if(col.name %in% names(df)){
        type.fun <- type.list[[col.name]]
        df[[col.name]] <- type.fun(df[[col.name]])
      }
    }
    df
  }else{
    match.mat
  }
### If type.list is a list of functions, then return a data.frame
### whose columns are defined by calling the functions in type.list on
### the corresponding column of match.mat. Otherwise just return a
### character matrix. If match.mat does not already have rownames, and
### it has a column named "name", then that column will be used for
### the rownames, and that column will not be returned.
}

names_or_error <- function
### Extract capture group names. Stop with an error if there are no
### capture groups, or if there are any capture groups without names.
(vec.with.attrs
### Output from g?regexpr.
 ){
  capture.names <- attr(vec.with.attrs, "capture.names")
  if(!is.character(capture.names) || any(capture.names == "")){
    stop("pattern must contain named capture groups (?<name>subpattern)")
  }
  capture.names
### Character vector.
}

regions.by.file <- list()
regions.by.chunk.file <- list()
chunk.limits.list <- list()
bed.list <- list()
positive.regions.list <- list()
for(labels.file in labels.file.vec){
  ## there should be bigwig files in subdirectories under the same
  ## directory as labels.file.
  sample.dir.vec <- Sys.glob(file.path(project.dir, "samples", "*", "*"))
  sample.id <- basename(sample.dir.vec)
  sample.group.dir <- dirname(sample.dir.vec)
  sample.group <- basename(sample.group.dir)
  sample.df <- data.frame(sample.id, sample.group)
  samples.by.group <- split(sample.df, sample.df$sample.group)

  cat("Reading ", labels.file, "\n", sep="")
  labels.lines <- readLines(labels.file)
  is.blank <- labels.lines == ""
  chunk.id <- cumsum(is.blank)+1L
  label.df <- data.frame(chunk.id, line=labels.lines)[!is.blank, ]
  cat(length(unique(label.df$chunk.id)), " chunks, ",
      nrow(label.df), " label lines\n", sep="")

  ## Error checking.
  raw.vec <- paste(label.df$line)
  line.vec <- gsub(",", "", raw.vec)
  match.mat <- str_match_named(line.vec, g.pos.pattern)
  stopifnot(!is.na(match.mat[,1]))
  not.recognized <- ! match.mat[, "annotation"] %in% names(label.colors)
  if(any(not.recognized)){
    print(raw.vec[not.recognized])
    print(match.mat[not.recognized, , drop=FALSE])
    stop("unrecognized annotation")
  }
  match.df <-
    data.frame(chrom=match.mat[, "chrom"],
               chromStart=as.integer(match.mat[, "chromStart"]),
               chromEnd=as.integer(match.mat[, "chromEnd"]),
               annotation=match.mat[, "annotation"],
               sample.groups=match.mat[, "sample_groups"],
               chunk.id=label.df$chunk.id,
               stringsAsFactors=FALSE)
  match.by.chrom <- split(match.df, match.df$chrom)
  for(chrom in names(match.by.chrom)){
    chrom.df <- match.by.chrom[[chrom]]
    sorted <- chrom.df[with(chrom.df, order(chromStart, chromEnd)), ]
    same.as.next <- diff(sorted$chromStart) <= 0
    if(any(same.as.next)){
      bad.i <- which(same.as.next)
      print(sorted[c(bad.i, bad.i+1), ])
      stop("chromStart not increasing")
    }
    if(any(with(sorted, chromStart >= chromEnd))){
      print(sorted)
      stop("chromStart >= chromEnd")
    }
    overlaps.next <-
      with(sorted, chromStart[-1] < chromEnd[-length(chromEnd)])
    if(any(overlaps.next)){
      print(data.frame(sorted, overlaps.next=c(overlaps.next, FALSE)))
      stop("overlapping regions")
    }
  }

  ## determine total set of sample groups with positive=Peak
  ## annotations.
  stripped <- gsub(" *$", "", gsub("^ *", "", match.df$sample.groups))
  is.noPeaks <- stripped == ""
  match.df[is.noPeaks, "annotation"] <- "noPeaks"
  commas <- gsub(" +", ",", stripped)
  sample.group.list <- strsplit(commas, split=",")
  bed.list[[labels.file]] <- 
    data.frame(match.df[,c("chrom", "chromStart", "chromEnd")],
               name=paste0(match.df$annotation, ":", commas),
               score=0,
               strand=".",
               thickStart=match.df$chromStart,
               thickEnd=match.df$chromEnd,
               itemRgb=label.colors[paste(match.df$annotation)])
  names(sample.group.list) <- rownames(match.df)
  sample.group.vec <- unique(unlist(sample.group.list))
  cat("sample groups with peak annotations: ",
      paste(sample.group.vec, collapse=", "),
      "\n",
      sep="")

  ## Create some labeled regions for specific/nonspecific peaks.
  groups.up.vec <- sapply(sample.group.list, length)
  file.positive.regions <- match.df[0 < groups.up.vec,]
  input.has.peak <- grepl("Input", file.positive.regions$sample.groups)
  if(any(input.has.peak)){                         
    positive.regions.list[[labels.file]] <-
      with(file.positive.regions, data.frame(
        chrom, regionStart=chromStart, regionEnd=chromEnd,
        annotation, input.has.peak
        ))
  }

  match.by.chunk <- split(match.df, match.df$chunk.id)
  for(chunk.id in names(match.by.chunk)){
    ## Check that all regions are on the same chrom.
    chunk.df <- match.by.chunk[[chunk.id]]
    chunkChrom <- paste(chunk.df$chrom[1])
    if(any(chunk.df$chrom != chunkChrom)){
      print(chunk.df)
      stop("each chunk must span only 1 chrom")
    }
    regions.list <- list()
    for(ann.i in 1:nrow(chunk.df)){
      chunk.row <- chunk.df[ann.i, ]
      groups.up.vec <- sample.group.list[[rownames(chunk.row)]]
      is.observed <- sample.group.vec %in% groups.up.vec
      observed <- sample.group.vec[is.observed]
      not.observed <- sample.group.vec[!is.observed]
      to.assign <- list()
      ann <- chunk.row$annotation
      to.assign[observed] <- ann
      to.assign[not.observed] <- "noPeaks"
      for(sample.group in names(to.assign)){
        relevant.samples <- samples.by.group[[sample.group]]
        if(length(relevant.samples) == 0){
          glob.str <- file.path(sample.group.dir, "*")
          stop("no ", glob.str, " directories (but labels are present)")
        }
        annotation <- to.assign[[sample.group]]
        regions.list[[paste(ann.i, sample.group)]] <- 
          data.table(sample.id=paste(relevant.samples$sample.id),
                     sample.group,
                     chrom=chunk.row$chrom,
                     chromStart=chunk.row$chromStart,
                     chromEnd=chunk.row$chromEnd,
                     annotation)
      }
    }
    one.chunk <- do.call(rbind, regions.list)
    setkey(one.chunk, chromStart, chromEnd)
    file.and.chunk <- paste0(basename(labels.file), "-chunk", chunk.id)
    regions.by.chunk.file[[file.and.chunk]] <- one.chunk
    regions.by.file[[labels.file]][[chunk.id]] <- one.chunk
    chunk.limits.list[[file.and.chunk]] <- with(one.chunk, {
      data.frame(file.and.chunk,
                 chrom=chrom[1],
                 chromStart=min(chromStart),
                 chromEnd=max(chromEnd))
    })
  }
}
bed <- do.call(rbind, bed.list)
chunk.limits <- do.call(rbind, chunk.limits.list)
positive.regions <- do.call(rbind, positive.regions.list)
rownames(chunk.limits) <- NULL
rownames(bed) <- NULL
rownames(positive.regions) <- NULL

## Save positive regions for filtering final peaks.
## positive.regions.RData <- file.path(data.dir, "positive.regions.RData")
## save(positive.regions, file=positive.regions.RData)

## Save labels to bed file for viewing on UCSC.
bed.gz <- file.path(project.dir, "all_labels.bed.gz")
con <- gzfile(bed.gz, "w")
header <- 
  paste("track",
        "visibility=pack",
        "name=PeakSegJointLabels",
        'description="Visually defined labels',
        'in regions with and without peaks"',
        "itemRgb=on")
writeLines(header, con)
write.table(bed, con,
            row.names=FALSE,
            col.names=FALSE,
            quote=FALSE)
close(con)

limits.by.chrom <- split(chunk.limits, chunk.limits$chrom)
for(chrom in names(limits.by.chrom)){
  chrom.limits <- limits.by.chrom[[chrom]]
  ## Find overlapping chunks, and join them:
  clustered <- clusterPeaks(chrom.limits)
  limits.by.cluster <- split(clustered, clustered$cluster)
  chunks.per.cluster <- sapply(limits.by.cluster, nrow)
  not.ok <- 1 < chunks.per.cluster
  if(any(not.ok)){
    print(limits.by.cluster[not.ok])
    stop("chunks in different label files should not overlap")
  }
}

## Write labels to each sample.
all.regions <- do.call(rbind, regions.by.chunk.file)
regions.by.sample <- split(all.regions, all.regions[, paste0(sample.group, "/", sample.id)])
for(sample.path in names(regions.by.sample)){
  sample.labels <- regions.by.sample[[sample.path]]
  labels.bed <- file.path(project.dir, "samples", sample.path, "labels.bed")
  write.table(
    sample.labels[, .(chrom, chromStart, chromEnd, annotation)],
    labels.bed,
    quote=FALSE,
    row.names=FALSE,
    col.names=FALSE,
    sep="\t")
  cat("Wrote ", nrow(sample.labels),
      " labels to ", labels.bed,
      "\n", sep="")
}
