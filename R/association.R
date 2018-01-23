plot_association <- function(cones, causal=FALSE) {

  cones <- mutate(cones, snptype = ifelse(selected, "Selected", "Non-selected"))

  if (any(causal)) {
    cones <- mutate(cones, snptype = ifelse(causal, paste(snptype, "causal"), paste(snptype, "non-causal")))
  }

  ggplot(cones, aes(x = snptype, y = c, fill = selected)) +
    geom_boxplot() +
    labs(x = "SNP type", y = "Association score")

}

compute_association <- function(gwas, snps = NULL) {

  X <- as(gwas$genotypes, "numeric")
  Y <- gwas$fam$affected

  if (is.logical(snps)) {
    X <- X[,snps]
  } else if (class(snps) == "igraph.vs") {
    snps <- gwas$map$snp.names %in% names(snps)
    X <- X[,snps]
  }

  association <- apply(X, 2, function(x){
    df <- data.frame(p = Y, g = x)
    suppressWarnings(chsq <- chisq.test(table(df)))
    data.frame(chi = chsq$statistic, p = chsq$p.value)
  }) %>%
    do.call(rbind, .) %>%
    mutate(snp = colnames(X))

  return(association)

}

