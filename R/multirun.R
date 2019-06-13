#' Consistency among different experiments
#'
#' @description Plot a heatmap displaying the consistency in the selection of different experiments.
#'
#' @param ... Named outputs of martini e.g. "AdditiveModel = cones_am, RecessiveModel = cones_re".
#' @return A heatmap with each of the SNPs selected in any of the experiments, and where they were selected.
#' @importFrom dplyr filter mutate starts_with
#' @importFrom tidyr gather
#' @importFrom ggplot2 element_blank ggplot geom_tile scale_fill_manual theme aes facet_grid labs
#' @importFrom magrittr %>%
#' @export
consistency <- function(...) {

  cones <- join_experiments(...)

  gather(cones, key = "experiment", value = "selected", starts_with("selected")) %>%
  mutate(selected = ifelse(is.na(selected), FALSE, selected),
         experiment = gsub("selected_", "", experiment),
         experiment = factor(experiment, levels = unique(experiment))) %>%
  ggplot(aes(x = experiment, y = as.character(pos), fill = ifelse(selected, "Yes", "No"))) +
    geom_tile() +
    labs(x = "Experiment", y = "Genomic position", fill = "Selected") +
    scale_fill_manual(values=c("Yes"="gray20", "No"="gray90")) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    facet_grid(chr ~ ., scales = "free")

}

#' Join cones
#'
#' @description Join the outputs of martini for different experiments
#'
#' @param ... Named outputs of martini e.g. '"Additive model" = cones_am, Recessive = cones_re'.
#' @return Dataframe containing all the experiments.
#' @importFrom dplyr full_join
#' @export
join_experiments <- function(...) {

  experiments <- list(...)
  labels <- names(experiments)
  for (i in 1:length(experiments)) {
    e <- subset(experiments[[i]], selected)
    colnames(e)[colnames(e) == 'c'] <- paste0("C_", labels[i])
    colnames(e)[colnames(e) == 'selected'] <- paste0("selected_", labels[i])
    colnames(e)[colnames(e) == 'module'] <- paste0("modules_", labels[i])

    if (i == 1) {
      cones <- e
    } else {
      cones <- full_join(cones, e, by = c("snp","chr","cm","pos","allele.1","allele.2"), all = T)
    }
  }

  cones <- cones %>%
    mutate(consistency = rowSums(.[grep("selected", names(.))], na.rm = TRUE)/length(experiments))

  return(cones)

}
