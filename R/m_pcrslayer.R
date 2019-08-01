#' Detecting dysfunctional PCRs
#'
#' Detecting dysfunctional PCRs based on PCR replicates.
#'
#'
#' @param metabarlist   a \code{\link{metabarlist}} object
#' @param replicates    a factor of sample names within which replicates belong (i.e. should be a prefix/suffix of pcr replicate names)
#' @param control       a character vector indicating the names of experimental controls against which PCR replicates are to be compared.
#' @param thresh.method a character indicating which method should be used to define the filtering threshold.
#' @param plot          a boolean indicating whether dissimilarity distribution should be plotted. Default is \code{TRUE}
#' @param wthn.btwn     an ouput from \code{pcr_within_between}
#' @param colvec        a grouping factor for coloring PCR replicates in \code{check_pcr_repl}
#' @param dyspcr        a character vector of the dysfunctional PCR identified by \code{pcrslayer}
#' @name pcrslayer
#'
#' @return a vector of dysfunctional PCRs
#'
#' @details
#'
#' The \code{pcrslayer} function identifies potential non-functional PCR reactions by comparing the dissimilarities in OTU composition within a biological sample (i.e. between PCR replicates, hereafter \emph{dw}) vs. between biological samples (hereafter \emph{db}). It relies on the assumption that PCR replicates from a same biological samples should be more similar than two different biological samples (\emph{dw} < \emph{db}). More specifically, the function consists in first constructing an average OTU community for each biological sample by averaging the OTUs abundances of PCR replicates from the same biological sample. Dissimilarities \emph{dw} are then defined as the pairwise Bray-Curtis dissimilarities between PCR replicates with their associated average OTU community. Dissimilarities \emph{db} correspond to the pairwise Bray-Curtis dissimilarities between average OTU communites from the different biological samples. A PCR replicate having a \emph{dw} above a given dissimilarity threshold \emph{tresh} is considered to be too distant from its associated average OTU community and are excluded from the analysis. The whole process is repeated iteratively until no more PCR are excluded from the analysis. If only one single PCR replicate is representative of a biological sample after this trimming, it is also considered as a dysfunctional PCR.
#'
#' The threshold \emph{tresh} is defined automatically with two alternative methods. Either it is the intersection of \emph{dw} and \emph{db} distributions (\code{tresh.method="interesect"}). Or it is the mode of the \emph{db} distribution (\code{tresh.method="mode"}).
#'
#' Function \code{check_pcr_thresh} enables visualization of \emph{dw} and \emph{db} distributions. Function \code{check_pcr_repl} enables visualization of PCR replicate dissimilarity patterns in a NMDS ordination and distance from their average OTU community.
#'
#' Function \code{pcr_control} is another way of detecting dysfunctional PCRs, and considers that any PCR replicate that is too similar to any control amplicon (blank or mock community) is dysfunctional. Note that this function will not be appropriate if one or more controls are contaminated with the DNA from biological samples (e.g. cross-contaminations). ### TO FINISH
#'
#' @examples
#' library(ggplot2)
#' data(soil_euk)
#' # define replicate factor
#' # Consider only biological samples
#' sample_subset <- subset_metabarlist(soil_euk, "pcrs", rownames(soil_euk$pcrs)[which(soil_euk$pcrs$type == "sample")])
#'
#' # first visualization
#' comp1 <- pcr_within_between(sample_subset$reads, replicates = sample_subset$pcrs$sample_id)
#' check_pcr_thresh(comp1, thresh.pcr = NULL)
#' # visualization of replicates through NMDS
#' nmds <- check_pcr_repl(sample_subset$reads,
#'   replicates = sample_subset$pcrs$sample_id,
#'   colvec = paste(sample_subset$samples$Habitat, sample_subset$samples$Material, sep = "|")
#' )
#' nmds + labs(fill = "sample type")
#'
#' # identify dysfunctional PCRs
#' bad_pcrs <- pcrslayer(sample_subset,
#'   thresh.method = "intersect",
#'   replicates = sample_subset$pcrs$sample_id
#' )
#'
#' nmds <- check_pcr_repl(sample_subset$reads,
#'   replicates = sample_subset$pcrs$sample_id,
#'   colvec = paste(sample_subset$samples$Habitat, sample_subset$samples$Material, sep = "|"),
#'   dyspcr = bad_pcrs
#' )
#' nmds + labs(fill = "sample type")
#'
#' # identify PCRs too close to negative controls ### TO FINISH
#' @author Lucie Zinger
#' @importFrom vegan decostand vegdist metaMDS
#' @export pcrslayer
#' @export pcr_within_between
#' @export check_pcr_thresh
#' @export check_pcr_repl
#'

pcrslayer <- function(metabarlist, replicates = metabarlist$pcrs$sample_id, thresh.method = "intersect", plot = T) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads_table <- metabarlist$reads

    if (nrow(reads_table) != length(replicates)) {
      stop("reads table and replicates must have the same length")
    }

    if (thresh.method != "intersect" & thresh.method != "mode") {
      stop('thresh.method should be one of "intersect" or "mode"')
    }

    subset_data <- data.frame(
      groups = replicates, replicating = TRUE,
      row.names = rownames(metabarlist$pcrs)
    )

    idx <- which(rowSums(reads_table) == 0)
    subset_data[idx, "replicating"] <- FALSE

    iteration <- 0
    repeat{
      iteration <- iteration + 1
      print(paste("Iteration", iteration))

      reads_table <- reads_table[
        rownames(subset_data[subset_data$replicating, ]),
      ]

      replicates <- subset_data[rownames(reads_table), "groups"]

      nb_bad_pcr <- length(rownames(subset_data[subset_data$replicating==F, ]))
      wthn_btwn <- pcr_within_between(reads_table, replicates)
      thresh_pcr <- pcr_threshold_estimate(wthn_btwn, thresh.method)
      if (plot == T) {
        check_pcr_thresh(wthn_btwn, thresh_pcr)
      }
      # ajoute à la liste les samples dont la valeur (dist par rapport au barycentre) est supérieure au thresh_pcr
      # et dont la valeur est égale à la valeur max des réplicats
      to_flag <- unname(unlist(lapply(wthn_btwn$pcr_intradist, function(y) {
        names(which(y > thresh_pcr & y == max(y)))
      })))
      subset_data[to_flag, "replicating"] <- FALSE

      # recherche les singletons dans la nouvelle matrice de reads
      reads_table <- reads_table[
        rownames(subset_data[subset_data$replicating, ]),
      ]
      replicates <- subset_data[rownames(reads_table), "groups"]
      # dans le cas où l'on n'a qu'un seul réplicat, on ajoute a bad_pcr les rownames du réplicat
      if (length(which(table(as.vector(replicates)) < 2)) != 0) {
        singletons <- sapply(
          names(which(table(as.vector(replicates)) < 2)),
          function(x) grep(x, rownames(reads_table))
        )
        singletons_ids <- rownames(reads_table)[unname(singletons)]
        subset_data[singletons_ids, "replicating"] <- FALSE
      }

      # stop the loop when none of replicat is added to the vector bad_pcr
      if (length(rownames(subset_data[subset_data$replicating==F,])) == nb_bad_pcr) {
        break
      }
    }
    return(subset_data)
  }
}


# pcr_controls = function(x, replicates, thresh.pcr, plot=T) {
#
#   x.n = decostand(x, MARGIN = 1, "total")
#   x.dist = vegdist(x.n, "bray")
#
#   x.dist.list = data.frame(t(combn(labels(x.dist),2)), dist=as.vector(x.dist))
#   x.dist.list$X1.type = replicates[match(x.dist.list$X1, rownames(x))]
#   x.dist.list$X2.type = replicates[match(x.dist.list$X2, rownames(x))]
#
#   #here intradist should be within control  while between should be between samples and controls
#   x.dist.list$comp = ifelse(x.dist.list$X1.type==x.dist.list$X1.type & x.dist.list$X1 %in% control, "control",
#                             ifelse(x.dist.list$X1.type==x.dist.list$X2.type &
#                                      !x.dist.list$X1.type %in% control, "sample", "between"))
#   out = list(pcr_intradist = x.dist.list$dist[which(x.dist.list$comp=="control")],
#              bar_dist = x.dist.list$dist[which(x.dist.list$comp=="between")])
#
#   thresh.pcr = pcr_threshold_estimate(out, thresh.pcr)
#   if(plot==T) {
#     check_pcr_thresh(out, thresh.pcr)
#   }
#   bad.all = unique(sort(unlist(x.dist.list[x.dist.list$comp=="between" &
#                                              x.dist.list$dist < thresh.pcr, c("X1", "X2")])))
#   bad_pcr = as.vector(bad.all[which(bad.all %in% rownames(x.n)[which(replicates==control)]==F)])
#   return(bad_pcr)
# }


pcr_within_between <- function(reads, replicates) {
  # barycentre calculation and intradist function

  # data standardization
  reads_stdt <- decostand(reads, MARGIN = 1, "total")

  bar <- t(sapply(by(reads_stdt, as.vector(replicates), colMeans), identity))

  # between barycentre distances
  bar_dist <- vegdist(bar, "bray")

  # within replicates distances
  pcr_intradist <- lapply(1:nrow(bar), function(x) {
    ind <- which(replicates == rownames(bar)[x])
    sapply(ind, function(y) {
      out <- vegdist(rbind(bar[x, ], reads_stdt[y, ]), "bray")
      names(out) <- rownames(reads_stdt)[y]
      out
    })
  })
  names(pcr_intradist) <- rownames(bar)
  return(list(bar_dist = bar_dist, pcr_intradist = pcr_intradist))
}

pcr_threshold_estimate <- function(wthn.btwn, thresh.method = "intersect") {
  dinter.max <- max(wthn.btwn$bar_dist)
  ddinter <- density(wthn.btwn$bar_dist, from = 0, to = 1)
  ddintra <- density(unlist(wthn.btwn$pcr_intradist), from = 0, to = 1)

  # assumption that each of them have a "unimodal" distribution
  dintra.mode <- ddintra$x[which.max(ddintra$y)]
  dinter.mode <- ddinter$x[which.max(ddinter$y)]

  # get the mixed distrib within both mode intervals and find intersection between the two
  # inter.y=ddinter$y[which(ddinter$x>dintra.mode & ddinter$x<dinter.mode & ddinter$y>ddintra$y)]
  # intra.y=ddintra$y[which(ddintra$x>dintra.mode & ddintra$x<dinter.mode & ddintra$y>ddinter$y)]
  # x = ddinter$x[which(ddinter$x>dintra.mode & ddinter$x<dinter.mode)]
  # s = c(intra.y, inter.y)
  # return(x[which.min(s)[1]])

  if (thresh.method == "intersect") {
    p <- which(ddintra$y - ddinter$y > 0 & ddinter$y > dintra.mode)
    out <- ddinter$x[p[length(p)]]
  } else {
    out <- dinter.mode
  }
  return(out)
}

check_pcr_thresh <- function(wthn.btwn, thresh.pcr) {
  d.bar <- density(wthn.btwn$bar_dist)
  d.intra <- density(unlist(wthn.btwn$pcr_intradist))
  plot(d.bar, lwd = 2, main = "Distances density")
  lines(d.intra, col = "green", lwd = 2)
  legend("topright", c("btwn samples", "wthn samples"),
    col = c("black", "green"), lwd = 2
  )
  abline(v = thresh.pcr)
}

check_pcr_repl <- function(x, replicates, colvec = NULL, dyspcr = NULL) {
  x.n <- decostand(x, "total")
  bar <- t(sapply(by(x.n, as.vector(replicates), colMeans), identity))
  all <- rbind(x.n, bar)

  mds <- metaMDS(all)
  d <- data.frame(scores(mds, display = "sites"))
  d$points <- ifelse(rownames(d) %in% rownames(bar), "bary", "samp")
  dsub <- d[d$points == "bary", ]
  d$NMDS1bary <- d$NMDS2bary <- NA
  d$NMDS1bary[d$points != "bary"] <- dsub$NMDS1[match(replicates, rownames(dsub))]
  d$NMDS2bary[d$points != "bary"] <- dsub$NMDS2[match(replicates, rownames(dsub))]

  d2 <- d[d$points != "bary", ]
  if (is.null(colvec) & is.null(dyspcr)) {
    ggplot(d2, aes(x = NMDS1, y = NMDS2)) +
      geom_point(shape = 21) + theme_bw() +
      geom_segment(data = d, aes(
        x = NMDS1bary, y = NMDS2bary,
        xend = NMDS1, yend = NMDS2
      ), colour = "grey")
  } else if (!is.null(colvec) & is.null(dyspcr)) {
    d$col <- NA
    d2$col <- colvec
    ggplot(d2, aes(x = NMDS1, y = NMDS2, fill = col)) +
      geom_point(shape = 21) + theme_bw() +
      geom_segment(data = d, aes(
        x = NMDS1bary, y = NMDS2bary,
        xend = NMDS1, yend = NMDS2
      ), colour = "grey")
  } else if (is.null(colvec) & !is.null(dyspcr)) {
    d$dyspcr <- NA
    d2$dyspcr <- ifelse(rownames(d2) %in% dyspcr, "dysfunctional", "functional")
    ggplot(d2, aes(x = NMDS1, y = NMDS2)) +
      geom_point(aes(shape = dyspcr)) + theme_bw() +
      geom_segment(data = d, aes(
        x = NMDS1bary, y = NMDS2bary,
        xend = NMDS1, yend = NMDS2
      ), colour = "grey") +
      scale_shape_manual(values = c(4, 21))
  } else {
    d$col <- NA
    d2$col <- colvec
    d$dyspcr <- NA
    d2$dyspcr <- ifelse(rownames(d2) %in% dyspcr, "dysfunctional", "functional")
    ggplot(d2, aes(x = NMDS1, y = NMDS2, fill = col)) +
      geom_point(aes(shape = dyspcr)) + theme_bw() +
      geom_segment(data = d, aes(
        x = NMDS1bary, y = NMDS2bary,
        xend = NMDS1, yend = NMDS2
      ), colour = "grey") +
      scale_shape_manual(values = c(4, 21))
  }
}
