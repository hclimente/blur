#' Ideogram of a SNP module.
#'
#' @description Create a ideogram of a SNP module from \code{SConES} results using the Gviz package (Hahne and Ivanek, 2016).
#'
#' @param cones Results from \code{SConES}.
#' @param k Id of the module to plot.
#' @param genome Abbreviations of the genome to use: hg19 for human (default),  mm10 for mouse, etc. Argument to be passed to
#' \code{\link[Gviz]{IdeogramTrack}}.
#' @return An ideogram per chromosome showing the selected SNPs and the genes in the region.
#' @references Hahne F. and Ivanek R. (2016). "Statistical Genomics: Methods and Protocols." In Mathe E and Davis S (eds.), chapter
#' Visualizing Genomic Data Using Gviz and Bioconductor, pp. 335-351. Springer New York, New York, NY. ISBN 978-1-4939-3578-9, doi:
#' 10.1007/978-1-4939-3578-9_16, \url{http://dx.doi.org/10.1007/978-1-4939-3578-9_16}.
#' @export
plot_snp_module <- function(cones, k, genome = "hg19") {

  martini:::check_installed("GenomeInfoDb", "plot_snp_module")
  martini:::check_installed("GenomicRanges", "plot_snp_module")
  martini:::check_installed("Gviz", "plot_snp_module")
  martini:::check_installed("IRanges", "plot_snp_module")

  module <- subset(cones, module %in% k)

  by(module, module$chr, function(snps2plot) {

    snpRange <- GenomicRanges::GRanges(seqnames = paste0("chr", snps2plot$chr),
                                       ranges = IRanges::IRanges(start = snps2plot$pos,
                                                                 end = snps2plot$pos) )
    GenomeInfoDb::genome(snpRange) <- genome
    tsnp <- Gviz::AnnotationTrack(snpRange, name = "SNP", stacking ="dense")

    tideo <- Gviz::IdeogramTrack(genome = GenomeInfoDb::genome(snpRange),
                                 chromosome = names(GenomeInfoDb::genome(snpRange)))
    tseq <- Gviz::GenomeAxisTrack()

    tbio <- Gviz::BiomartGeneRegionTrack(genome = genome,
                                         chromosome = names(GenomeInfoDb::genome(snpRange)),
                                         start = head(snps2plot$pos, n = 1),
                                         end = tail(snps2plot$pos, n = 1),
                                         name = "Ensembl")

    Gviz::plotTracks(list(tideo, tseq, tbio, tsnp), showId = TRUE)
  })

  return(TRUE)

}

#' Converts a MAP data.frame to a BED data.frame
#'
#' @description Takes a map file and:
#'  \itemize{
#' \item{column 1: Used as the chromosome column in the BED file..}
#' \item{column 4: Used as start and end in the BED data.frame (as we work with SNPs).}
#' }
#'
#' @param map A MAP data.frame.
#' @return A BED data.frame.
map2bed <- function(map) {

  bed <- subset(map, select = c("chr", "pos"))
  colnames(bed) <- c("chr", "start")
  bed$chr <- paste0("chr", bed$chr)
  bed$chr <- ifelse(bed$chr == "chr23", "chrX", bed$chr)
  bed$end <- bed$start

  return(bed)

}
