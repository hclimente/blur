#' Plot a subnetwork
#'
#' @description Plot a subnetwork.
#'
#' @param gwas A GWAS experiment.
#' @param net A network of SNPs.
#' @param snps Logical indicating the rows of gwas$map that describe the subnetwork.
#' @return A representation of the causal subnetwork.
#' @importFrom ggnetwork ggnetwork geom_edges geom_nodes theme_blank
#' @importFrom ggplot2 ggplot aes
#' @importFrom magrittr %>%
#' @export
plot_subnetwork <- function(gwas, net, snps) {

  if (is.logical(snps) & length(snps) == nrow(gwas$map)) {
    snps <- gwas$map$snp.names[snps]
  } else if (class(snps) == "igraph.vs") {
    snps <- names(snps)
  } else if (is.character(snps)) {
    snps <- snps
  } else {
    stop("blur doesn't know how to deal with data of class", class(snps))
  }

  martini:::subnet(net, "name", snps) %>%
    ggnetwork %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(color = "grey50") +
    geom_nodes(aes(color = gene)) +
    theme_blank()

}

get_affected_genes <- function(cones, net) {

  cones <- cones %>% filter(selected)
  genes <- martini:::subvert(net, 'name', cones$snp)$gene
  genes <- table(genes)

  return(genes)

}

#' @importFrom dplyr data_frame filter mutate group_by summarize ungroup
#' @importFrom magrittr %>%
#' @importFrom stats na.omit
#' @export
get_affected_regions <- function(cones, net) {

  V <- martini:::subvert(net, 'name', cones$snp)
  snp2gene <- data_frame(snp = V$name, gene = V$gene)

  filter(cones, selected) %>%
    left_join(snp2gene, by = "snp") %>%
    group_by(module, chr) %>%
    summarize(start = min(pos),
              end = max(pos),
              genes = unique(gene) %>% na.omit %>% paste(collapse = ",") ) %>%
    ungroup

}
