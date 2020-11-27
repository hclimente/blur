#' Histogram of module sizes.
#'
#' @description Plot an histogram with the size of the modules bigger than min_size.
#'
#' @param cones martini results.
#' @param min_size Minimum size of the modules.
#' @return An histogram of the module size.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot geom_histogram aes labs
#' @importFrom magrittr %>%
#' @export
module_size <- function(cones, min_size = 1) {

  compute_mod_size(cones) %>%
    filter(size > min_size) %>%
    ggplot(aes(x = size)) +
    geom_histogram(bins = 10) +
    labs(x = "Module size", y = "Counts")

}

#' Density plots of module sizes and the median association score of their SNPs.
#'
#' @description Density plots of module sizes and the median association score of their SNPs.
#'
#' @param cones martini results.
#' @return An density plot of the module size and the median association score.
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot stat_binhex aes labs theme
#' @importFrom magrittr %>%
#' @export
module_size_association <- function(cones) {

  ggplot(compute_mod_size(cones), aes(x = size, y = C)) +
    stat_binhex() +
    labs(x = "Module size", y = "Median association score", fill = "# modules") +
    theme(legend.position="bottom")

}

#' Histogram of module sizes.
#'
#' @description Plot an histogram with the size of the modules bigger than min_size.
#'
#' @param cones martini results.
#' @return Size of each SNP and median association score.
#' @importFrom dplyr arrange filter group_by n summarise
#' @importFrom stats median
compute_mod_size <- function(cones) {

  cones %>%
    filter(selected) %>%
    group_by(module, chr) %>%
    summarise(size = n(), C = median(c)) %>%
    arrange(-size)

}

#' @importFrom dplyr mutate summarise select arrange left_join
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges findOverlaps IRanges
#' @importFrom S4Vectors queryHits subjectHits
#' @export
table_regions <- function(cones, snp2gene) {

  cones <- filter(cones, selected)

  if(! 'X' %in% cones$chr) {
    cytobands <- mutate(cytobands, chr = ifelse(chr == 'X', 23, chr))
  }

  snps <- with(cones, GRanges(seqnames = chr, IRanges(pos, width=1, names = snp)))
  bands <- with(cytobands, GRanges(seqnames = chr, IRanges(start, end, names = band)))
  olaps <- findOverlaps(snps, bands)
  cones <- cbind(cones[queryHits(olaps),], select(cytobands[subjectHits(olaps),],-chr))

  cones %>%
    left_join(snp2gene, by = 'snp') %>%
    group_by(chr, band, module) %>%
    summarise(region = unique(paste0(chr, band)),
              coords = paste(min(pos), max(pos), sep = '-'),
              genes = paste(unique(na.omit(gene)), collapse = ','),
              noGeneSnps = sum(is.na(gene)),
              bestSnp = snp[which.max(c)],
              numSnps = paste0(length(unique(snp)), ' (', noGeneSnps, ')')) %>%
    ungroup %>%
    arrange(chr, coords) %>%
    select(region, module, coords, bestSnp, genes, numSnps)

}
