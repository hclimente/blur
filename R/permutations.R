#' Get a LD-accounted gene score value
#'
#' @description Creates a score for genes by aggregating the value of the SNPs
#' mapped to the gene. The score is corrected by length of the gene and by
#' LD-structure.
#'
#' @param cones martini results.
#' @param snp2gene Dataframe with two columns: the id of the SNP, and the gene annotated to it.
#' @return A gene score based on the SNPs mapped to it.
#' @importFrom dplyr select arrange
#' @importFrom magrittr %>%
#' @export
cgp <- function(cones, snp2gene) {

  # get permutations
  cones <- arrange(cones, chr, pos)

  C <- cones$c
  permutations <- lapply(1:length(C) - 1, function(i) {
    C[unique(c((i+1):length(C), 1:i))]
  }) %>% do.call(cbind, .)

  # convert snp2gene to a list of boolean mask
  colnames(snp2gene) <- c("snp","gene")
  snp2gene <- by(snp2gene, snp2gene$gene, function(g){
    cones$snp %in% g$snp
  }) %>% as.list

  # get gene-scores based on the permutations
  L <- ncol(permutations)
  geneC <- lapply(snp2gene, function(i) {
    sum(apply(permutations[i,], 2, max, na.rm=T) > max(C[i], na.rm=T))
  }) %>% stack()
  colnames(geneC) <- c("l","gene")
  geneC <- select(geneC, gene, l)
  geneC$p <- 1 - geneC$l/(L + 1)

  return(geneC)

}
