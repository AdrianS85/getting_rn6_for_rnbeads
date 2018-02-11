########################################################################################################################
## siteAnnotation.R
## created: 2012-10-26
## creator: Fabian Mueller
## ---------------------------------------------------------------------------------------------------------------------
## Functions dedicated to annotating loci relevant for bisulfite sequencing.
########################################################################################################################

## L I B R A R I E S ###################################################################################################
source("https://bioconductor.org/biocLite.R")
biocLite("BSgenome.Rnorvegicus.UCSC.rn5")

suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(RnBeads)) # get.cpg.stats, rnb.sort.regions

## G L O B A L S #######################################################################################################

## Chromosomes supported by the annotation packages
CHROMOSOMES <- list("hg19" = c(1:22, "X", "Y"), "mm9" = c(1:19, "X", "Y"), "mm10" = c(1:19, "X", "Y"),
                    "rn5" = c(1:20, "X"), "rn6" = c(1:20, "X", "Y")) ##TUTAJ DODA?EM rn6
for (assembly in names(CHROMOSOMES)) {
  CHROMOSOMES[[assembly]] <- paste0("chr", CHROMOSOMES[[assembly]])
  names(CHROMOSOMES[[assembly]]) <- CHROMOSOMES[[assembly]]
}
rm(assembly)

## Dinucleotide patterns to be annotated
NUCLEOTIDE.PATTERNS <- c("CG")
names(NUCLEOTIDE.PATTERNS) <- sapply(strsplit(NUCLEOTIDE.PATTERNS, ""), paste, collapse = "p")

## F U N C T I O N S ###################################################################################################

## get.genome.data
##
## Gets the specified genome.
##
## @param assembly Genome assembly of interest. Currently the only supported genomes are \code{"hg19"}, \code{"mm9"},
##                 \code{"mm10"} and \code{"rn5"}.
## @return Sequence data object for the specified assembly.
##
## @author Yassen Assenov
get.genome.data <- function(assembly) {
  if (assembly == "hg19") {
    suppressPackageStartupMessages(require(BSgenome.Hsapiens.UCSC.hg19))
    genome.data <- Hsapiens
  } else if (assembly == "mm9") {
    suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm9))
    genome.data <- Mmusculus
  } else if (assembly == "mm10") {
    suppressPackageStartupMessages(require(BSgenome.Mmusculus.UCSC.mm10))
    genome.data <- Mmusculus
  } else if (assembly == "rn5") {
    suppressPackageStartupMessages(require(BSgenome.Rnorvegicus.UCSC.rn5))
    genome.data <- Rnorvegicus
  } else if (assembly == "rn6") { ##TO DODA?EM
    suppressPackageStartupMessages(require(BSgenome.Rnorvegicus.UCSC.rn6)) 
    genome.data <- Rnorvegicus
  } else {
    stop("unsupported assembly")
  }
  return(genome.data)
}

########################################################################################################################

## rnb.update.sites
##
## Creates lists of Genomic ranges for genomic sites targeted by methylation.
##
## @param assembly            Genome assembly to use.
## @param cpgislands          Region annotation of CpG islands for the specified assembly.
## @param length.neighborhood Length, in base pairs, of the local neighborhood around the targeted dinucleotides. This
##                            neighborhood is used in calculating sequence-based statistics: CpG density and GC content.
## @return List of \code{\link{GRangesList}} objects. Every item is dedicated to a motif (e.g. CpGs) and groups
##         \code{\link{GRanges}} instances, one per chromosome.
## @author Fabian Mueller

XXX <- rnb.update.sites(assembly = "rn5", cpgislands = NULL, length.neighborhood = 100L) ##TUTEJ DODA?EM
rnb.update.sites <- function(assembly = "rn6", cpgislands = NULL, length.neighborhood = 100L) { ##TUTEJ ZMIENI?EM
  genome.data <- get.genome.data(assembly)
  chrom.lengths <- seqlengths(genome.data)[CHROMOSOMES[[assembly]]]
  sites <- list()
  pp.dnas <- DNAStringSet(NUCLEOTIDE.PATTERNS)
  for (i in names(NUCLEOTIDE.PATTERNS)){
    pp.p <- pp.dnas[[i]]
    pp.m <- reverseComplement(pp.p)
    curSites <- lapply(CHROMOSOMES[[assembly]], function(chrom) {
      matches.st <- lapply(list(pp.p, pp.m), function(x) {
        ranges(matchPattern(x, genome.data[[chrom]]))
      })
      cp.starts <- start(matches.st[[1]]) - length.neighborhood / 2L + 1L ##TO WYDAJE SI? NIE BY? INTEGEREM Z JAKIEGO? POWODU?
      cp.ends <- cp.starts + length.neighborhood - 1L ##TO PEWNIE TEZ NIE BEDZIE INTEGEREM
      cpg.stats <- suppressWarnings(get.cpg.stats(genome.data[[chrom]], cp.starts, cp.ends))
      matches.gr <- mapply(GRanges, seqnames = list(chrom), ranges = matches.st, strand = list("+", "-"),
                           "CpG" = list(cpg.stats[, "CpG"]), "GC" = list(cpg.stats[, "GC"]))
      matches.gr <- rnb.sort.regions(do.call(c, matches.gr))
      seqlevels(matches.gr) <- names(CHROMOSOMES[[assembly]])
      seqlengths(matches.gr) <- chrom.lengths
      matches.gr
    })
    sites[[i]] <- GRangesList(curSites)
    logger.status(c("Created site annotation for", i))
  }
  
  if (!is.null(cpgislands)) {
    sites <- rnb.update.site.annotation.with.cgistatus(sites, cpgislands)
    logger.status("Enriched sites with CpG island information")
  }
  return(sites)
}

########################################################################################################################

## rnb.update.site.annotation.with.cgistatus
##
## Enriches the annotation of genomic sites by adding columns showing their relation to CpG island status. This
## relation is one of \code{"Island"}, \code{"Shore"}, \code{"Shelf"}, \code{"Open Sea"}.
##
## @param sites            List of \code{GRangesList} objects, one per site type or probe annotation.
## @param cpgislands       \code{GRangesList} object specifying the locations of CpG islands.
## @param length.cgi.shore Length, in base pairs, of a CpG island shore.
## @param length.cgi.shelf Length, in base pairs, of a CpG island shelf.
## @return The modified \code{sites}.
##
## @author Fabian Mueller
rnb.update.site.annotation.with.cgistatus <- function(sites, cpgislands, length.cgi.shore = 2000L,
                                                      length.cgi.shelf = 2000L) {
  
  enrich.f <- function(ss, cgis) {
    cgi.relation <- rep.int("Open Sea", length(ss))
    if (length(cgis) != 0) {
      shores <- c(flank(cgis, length.cgi.shore - 1L, start = TRUE),
                  flank(cgis, length.cgi.shore - 1L, start = FALSE))
      shelves <- c(flank(shores, length.cgi.shelf - 1L, start = TRUE),
                   flank(shores, length.cgi.shelf - 1L, start = FALSE))
      cgi.relation[ss %in% shelves] <- "Shelf"
      cgi.relation[ss %in% shores] <- "Shore"
      cgi.relation[ss %in% cgis] <- "Island"
    }
    elementMetadata(ss)[, "CGI Relation"] <-
      factor(cgi.relation, levels = c("Open Sea", "Shelf", "Shore", "Island"))
    ss
  }
  
  chromosomes <- names(cpgislands)
  sites.enriched <- sites
  for (i in names(NUCLEOTIDE.PATTERNS)) {
    if (!(i %in% names(sites))) {
      next
    }
    if (!setequal(chromosomes, names(sites[[i]]))) {
      stop("incompatible sites and CGI ranges")
    }
    grl <- foreach(ss = as.list(sites[[i]][chromosomes]), cgis = as.list(cpgislands)) %dopar% enrich.f(ss, cgis)
    names(grl) <- chromosomes
    ## Fix issue with column being renamed to CGI.Relation
    grl <- unlist(GRangesList(grl))
    i.cgirelation <- which(colnames(mcols(grl)) == "CGI.Relation")
    if (length(i.cgirelation) != 0) {
      colnames(mcols(grl))[i.cgirelation] <- "CGI Relation"
    }
    grl <- GenomicRanges::split(grl, seqnames(grl))
    sites.enriched[[i]] <- grl
  }
  
  return(sites.enriched)
}