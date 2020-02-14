#' Generating rarefaction curves using Hill numbers
#'
#' These functions generate and plot rarefaction curves from a \code{\link{metabarlist}} object using the hill numbers framework (i.e. \eqn{^{q}D}), as well as Good's coverage index.
#'
#'
#' @param metabarlist       a \code{\link{metabarlist}} object
#' @param nboot              number of resampling events to estimate \eqn{^{q}D} at a given sequencing depth.
#' @param nsteps             number of steps between sample sizes for the rarefaction curves. Default is 10 steps.
#'
#' @return Function \code{hill_rarefaction} returns an object of class \code{"hill_rarefaction"}, which corresponds to a table of diversity indices for each pcr rarefied at each `nsteps`  sequencing depth, as well as the arguments `nboot` and `nsteps` to conduct the analysis.
#'
#' @details
#' \code{\link{hill_rarefaction}} builds a rarefaction analysis for each pcr of a \code{\link{metabarlist}} object using Hill numbers for q={0,1,2} (see Chao et al. 2014 for a review). These indices are equivalent to :
#' \itemize{
#' \item{Richness, for q=0}
#' \item{Exponential of the Shannon entropy, for q->1}
#' \item{Inverse of the Simpson index, for q=2}
#' }
#'
#' The function also returns the Good's coverage index (1-singletons/#reads). Note however that this index should be interpreted carefully in metabarcoding data:
#' #' \itemize{
#' \item{absolute singletons (across the whole metabarcoding dataset) are usually filtered out during the bioinformatic process (which is the case for the \code{\link{soil_euk}} data). The Good's coverage estimate returned here is only based on the number of singletons per pcr after this filtering process, so the true number of singletons is underestimated here.}
#' \item{This coverage index gives an assessment of the coverage of the amplicon diversity within a pcr: it includes remaining errors, etc.. The coverage of the genuine DNA fragment diversity in the biological sample is likely to be misestimated with this index.}
#' }
#'
#'
#' @references Chao, A., Chiu, C. H., & Jost, L. (2014). Unifying species diversity, phylogenetic diversity, functional diversity, and related similarity and differentiation measures through Hill numbers. Annual review of ecology, evolution, and systematics, 45, 297-324.
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # Create a subset of pcrs: only a subset of samples from the H20 plot
#'
#' soil_euk_h20 <- subset_metabarlist(soil_euk,
#'    table = "pcrs",
#'    indices = grepl("H20-[A-B]", rownames(soil_euk$pcrs)))
#'
#' # run rarefaction (use boot = 20 to limit computation time, should be more)
#' soil_euk_h20.raref <- hill_rarefaction(soil_euk_h20, nboot = 20, nsteps = 10)
#'
#' # plot the results
#' gghill_rarefaction(soil_euk_h20.raref)
#'
#' # plot the results while differenciating litter vs. soil samples
#' p <- gghill_rarefaction(soil_euk_h20.raref,
#'     group = soil_euk_h20$samples$Material[match(soil_euk_h20$pcrs$sample_id, rownames(soil_euk_h20$samples))])
#' p
#' p + scale_fill_manual(values = c("goldenrod4", "brown4", "grey")) +
#'   scale_color_manual(values = c("goldenrod4", "brown4", "grey")) +
#'   labs(color = "Material type")
#'
#' @author Lucie Zinger
#' @importFrom vegan diversity rrarefy
#' @import ggplot2
#' @import reshape2
#'
#' @describeIn hill_rarefaction Compute hill_rarefaction curves on a \code{\link{metabarlist}} object.
#' @export hill_rarefaction

hill_rarefaction <- function(metabarlist, nboot = 10, nsteps = 10) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads <- metabarlist$reads

    out <- do.call("rbind", lapply(rownames(reads), function(y) {
      # sample sizes
      if (sum(reads[y, ]) > 2000) {
        rsize <- c(
          5, 10, 100, 200, 500, 1000,
          round(seq(2000, sum(reads[y, ]), length.out = nsteps))
        )
      } else {
        rsize <- c(1, round(seq(2, sum(reads[y, ]), length.out = nsteps)))
      }

      # assess estimates at depths rsize
      out <- data.frame(
        pcr_id = y, reads = rsize,
        D0 = NA, D0.sd = NA,
        D1 = NA, D1.sd = NA,
        D2 = NA, D2.sd = NA,
        coverage = NA, coverage.sd = NA
      )
      for (s in 1:length(rsize)) {
        # rarefaction
        r <- lapply(1:nboot, function(z) {
          rtmp <- rrarefy(reads[y, ], rsize[s])
          rtmp[rtmp > 0]
        })
        out[s, "D0"] <- mean(sapply(r, function(z) length(z)))
        out[s, "D0.sd"] <- sd(sapply(r, function(z) length(z)))
        out[s, "D1"] <- exp(mean(sapply(r, function(z) diversity(z, index = "shannon"))))
        out[s, "D1.sd"] <- exp(sd(sapply(r, function(z) diversity(z, index = "shannon"))))
        out[s, "D2"] <- mean(sapply(r, function(z) diversity(z, index = "invsimpson")))
        out[s, "D2.sd"] <- sd(sapply(r, function(z) diversity(z, index = "invsimpson")))
        out[s, "coverage"] <- mean(sapply(r, function(z) 1 - length(which(z == 1)) / sum(z)))
        out[s, "coverage.sd"] <- sd(sapply(r, function(z) 1 - length(which(z == 1)) / sum(z)))
      }
      return(out)
    }))

    out <- list(samples = rownames(reads), hill_table = out, nboot = nboot, nsteps = nsteps)
    attr(out, "class") <- "hill_rarefaction"
    return(out)
  }
}


#' @describeIn hill_rarefaction Plot a object of class \code{"hill_rarefaction"}
#' @param hill_rar  an object of class \code{"hill_rarefaction"}.
#' @param group     a vector or factor giving the grouping of each pcr included in the \code{"hill_rarefaction"} object. Missing values will be treated as another group and a warning will be given. The elements should correspond to the pcrs included in the `hill_rar$samples` object. Default is `NULL` for no grouping.
#' @export gghill_rarefaction


gghill_rarefaction <- function(hill_rar, group = NULL) {
  if (class(hill_rar) != "hill_rarefaction") {
    stop("hill_rar must be a hill_rarefaction object")
  }

  hill_table <- hill_rar$hill_table
  b <- melt(hill_table[, -grep("\\.sd", colnames(hill_table))], id.var = c("reads", "pcr_id"))
  b$value.sd <- melt(hill_table[, grep("\\.sd|reads|pcr_id", colnames(hill_table))],
    id.var = c("reads", "pcr_id")
  )$value

  if (is.null(group)) {
    ggplot(b, aes(x = reads, y = value, group = pcr_id)) +
      geom_line() +
      geom_ribbon(aes(ymin = value - value.sd, ymax = value + value.sd), alpha = 0.3) +
      facet_wrap(~variable, scale = "free", ncol = 4) +
      labs(x = "#reads", y = "diversity / coverage estimate") +
      theme_bw() +
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  } else {
    if (length(group) != length(hill_rar$sample)) {
      stop("The names of each elements of group should be defined and correspond to
                the non duplicated names of pcr included in the hill_rar object.")
    }

    b$group <- group[match(b$pcr_id, hill_rar$samples)]

    ggplot(b, aes(x = reads, y = value, group = pcr_id)) +
      geom_line(aes(color = group)) +
      geom_ribbon(aes(ymin = value - value.sd, ymax = value + value.sd, fill = group), alpha = 0.3) +
      facet_wrap(~variable, scale = "free", ncol = 4) +
      labs(x = "#reads", y = "diversity / coverage estimate") +
      theme_bw() +
      guides(fill = FALSE) +
      theme(
        panel.grid = element_blank(),
        legend.position = "bottom"
      )
  }
}
