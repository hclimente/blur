#' Ideogram of SConES results.
#'
#' @description Create a circular ideogram of the \code{SConES} results using the circlize package (Gu et al., 2014).
#'
#' @param cones Output from \code{SConES}.
#' @param genome Abbreviations of the genome to use: hg19 for human (default), mm10 for mouse, etc. Argument to be passed to
#' \code{\link[circlize]{circos.initializeWithIdeogram}} \code{species}.
#' @return A circular ideogram, including the manhattan plot, and the interactions between the selected SNPs.
#' @references Gu, Z., Gu, L., Eils, R., Schlesner, M., & Brors, B. (2014). circlize Implements and enhances circular visualization in R.
#' Bioinformatics (Oxford, England), 30(19), 2811-2. \url{https://doi.org/10.1093/bioinformatics/btu393}
#' @export
plot_ideogram <- function(cones, genome = "hg19") {

  martini:::check_installed("circlize", "plot_ideogram")

  circlize::circos.initializeWithIdeogram(species = genome)

  bed <- map2bed(cones)
  bed$c <- cones$c
  bed$selected <- cones$selected
  # order to put the selected snps in fron in the plot
  bed <- bed[with(bed, order(selected)),]

  circlize::circos.genomicTrackPlotRegion(bed,
                                          ylim = c(0, 1.1 * max(bed$c, na.rm = T)),
                                          panel.fun = function(region, value, ...) {
                                            # color according to selection/non-selection
                                            col = ifelse(value[[2]], "orange", "gray70")
                                            circlize::circos.genomicPoints(region, value, col = col, cex = 0.5, pch = 16)
                                          }, track.height = 0.3)

  # create links
  selected <- subset(cones, selected)

  regions <- by(selected, selected[,c("chr","module")], function(k) {
    data.frame(chr = paste0("chr", unique(k$chr)),
               module = unique(k$module),
               start = min(k$pos),
               end = max(k$pos))
  }) %>% do.call(rbind, .)

  links <- merge(regions, regions, by = "module")

  region1 <- subset(links, chr.x != chr.y, select = c("chr.x","start.x","end.x"))
  region2 <- subset(links, chr.x != chr.y, select = c("chr.y","start.y","end.y"))

  circlize::circos.genomicLink(region1, region2, col = sample(1:5, nrow(region1), replace = TRUE))

  circlize::circos.clear()

}
