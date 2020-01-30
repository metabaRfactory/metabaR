#' Detecting dysfunctional PCRs based on their reproducibility.
#'
#' Detecting dysfunctional PCRs based on their reproducibility.
#'
#'
#' @param metabarlist   a \code{\link{metabarlist}} object
#' @param replicates    a column name in the `pcrs`table corresponding to the sample names to which
#'                      pcr replicates belongs. Default is the `sample_id` column of the table `pcrs` from the
#'                      \code{\link{metabarlist}} object.
#' @param thresh.method a character indicating which method should be used to define the filtering threshold.
#' @param plot          a boolean indicating whether dissimilarity distribution should be plotted.
#'                      Default is \code{TRUE}
#' @param wthn.btwn     an ouput from \code{pcr_within_between}
#' @param groups        a column name in the `pcrs`table corresponding to a factor giving the groups for which the
#'                      graphical colors are drawn.
#' @param dyspcr        a character vector of the dysfunctional PCR identified by \code{pcrslayer}
#' @param control       a character vector indicating the names of controls against which PCR replicates are to be compared.
#'
#' @details
#'
#' The \code{pcrslayer} function identifies potential non-functional PCR reactions based on their reproducibility. It compares the dissimilarities in OTU composition within a biological sample (i.e. between PCR replicates, hereafter \emph{dw}) vs. between biological samples (hereafter \emph{db}). It relies on the assumption that PCR replicates from a same biological samples should be more similar than two different biological samples (\emph{dw} < \emph{db}). A PCR replicate having a \emph{dw} above a given dissimilarity threshold \emph{tresh} is considered to be too distant from its associated average OTU community and are excluded from the analysis. The whole process is repeated iteratively until no more PCR are excluded from the analysis. If only one single PCR replicate is representative of a biological sample after this trimming, it is also considered as a dysfunctional PCR.
#'
#' The \code{pcr_within_between} function computes dissimilarities in OTU composition within a biological sample \emph{dw} and between biological samples \emph{db}. It first consists inconstructing an average OTU community for each biological sample by averaging the OTUs abundances of PCR replicates from the same biological sample. Dissimilarities \emph{dw} are then defined as the pairwise Bray-Curtis dissimilarities between PCR replicates with their associated average OTU community. Dissimilarities \emph{db} correspond to the pairwise Bray-Curtis dissimilarities between average OTU communites from the different biological samples.
#'
#'
#' The threshold \emph{tresh} is defined automatically with two alternative methods. Either it is the intersection of \emph{dw} and \emph{db} distributions (\code{tresh.method="interesect"}). Or it is the mode of the \emph{db} distribution (\code{tresh.method="mode"}).
#'
#' Function \code{check_pcr_thresh} enables visualization of \emph{dw} and \emph{db} distributions. Function \code{check_pcr_repl} enables visualization of PCR replicate dissimilarity patterns in a NMDS ordination and distance from their average OTU community.
#'
#' Function \code{check_pcr_repl} enables visualization of dissimilarity patterns across all pcrs while showing pcr replicates centroidsthrough a Principal Coordinate Analysis (PCoA) based on Bray-Curtis dissimilarities.
#'
#' Function \code{pcr_control} is another way of detecting dysfunctional PCRs, and considers that any PCR replicate that is too similar to any control amplicon (blank or mock community) is dysfunctional. Note that this function will not be appropriate if one or more controls are contaminated with the DNA from biological samples (e.g. cross-contaminations). ### TO FINISH
#'
#' @return
#'
#' The \code{pcrslayer} function returns a vector of dysfunctional PCRs.
#'
#' The \code{pcr_within_between} function returns a list of dissimilarities \emph{dw} and \emph{db}.
#'
#' The \code{check_pcr_thresh} and  \code{check_pcr_repl} functions returns \code{ggplot} objects.
#'
#' @examples
#' library(ggplot2)
#' data(soil_euk)
#' # consider only biological samples
#' soil_euk_sub <- subset_metabarlist(soil_euk, "pcrs",
#'                                     rownames(soil_euk$pcrs)[which(soil_euk$pcrs$type == "sample")])
#'
#' # Visualization of within vs. between sample dissimilarities
#' soil_euk_sub_wb <- pcr_within_between(soil_euk_sub, "sample_id")
#' check_pcr_thresh(soil_euk_sub_wb, thresh.method = "intersect")
#'
#' # visualization of replicates through PCoA
#' ## create grouping factor according to habitat and material
#' soil_euk_sub$pcrs$habitat_material <- soil_euk_sub$pcrs$sample_id
#' idx = match(levels(soil_euk_sub$pcrs$habitat_material), rownames(soil_euk_sub$samples))
#' levels(soil_euk_sub$pcrs$habitat_material) <- paste(soil_euk_sub$samples$Habitat[idx],
#'                                                     soil_euk_sub$samples$Material[idx], sep = " | ")
#' ## vizualize dissimilarity patterns
#' mds <- check_pcr_repl(soil_euk_sub, replicates = "sample_id", groups = "habitat_material")
#' mds + labs(color = "sample type")
#'
#' # identify dysfunctional PCRs
#' bad_pcrs <- pcrslayer(soil_euk_sub, thresh.method = "intersect", replicates = "sample_id")
#'
#' # define a color vector that corresponds to the different habitat types. Should be named
#' mds <- check_pcr_repl(soil_euk_sub, replicates = "sample_id", groups = "habitat_material", dyspcr = bad_pcrs)
#' mds + labs(color = "sample type")
#'
#' # identify PCRs too close to negative controls ### TO FINISH
#' @author Lucie Zinger, Clement Lionnet
#' @importFrom vegan vegdist
#' @describeIn pcrslayer Detect dysfunctional PCRs in a \code{\link{metabarlist}} object.
#' @export pcrslayer

pcrslayer <- function(metabarlist, replicates, thresh.method = "intersect", plot = T) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads_table0 <- metabarlist$reads

    if (!replicates %in% colnames(metabarlist$pcrs)) {
      stop("replicates should be in colnames of the metabarlist$pcr table")
    }

    replicates0 <- metabarlist$pcrs[,replicates]

    if (thresh.method != "intersect" & thresh.method != "mode") {
      stop('thresh.method should be one of "intersect" or "mode"')
    }

    #vector of pcrs boolean: good = T bad = F
    good_pcrs <- rep(T, length(replicates0))
    names(good_pcrs) <- rownames(metabarlist$pcrs)

    #tag empty pcrs
    idx <- which(rowSums(reads_table0) == 0)
    good_pcrs[idx] <- FALSE

    iteration <- 0
    repeat{
      iteration <- iteration + 1
      print(paste("Iteration", iteration))

      #get good pcrs
      reads_table <- reads_table0[names(good_pcrs)[good_pcrs==T], ]

      replicates <- replicates0[which(good_pcrs==T)]

      nb_bad_pcr <- sum(good_pcrs==F)
      wthn_btwn <- pcr_within_between_internal(reads_table, replicates)
      thresh_pcr <- pcr_threshold_estimate(wthn_btwn, thresh.method)
      if (plot == T) {
        p = check_pcr_thresh(wthn_btwn, thresh.method)
        print(p)
      }

      #get names of pcrs not meeting quality criterion (max distance and distance > thresh_pcr)
      to_flag <- unname(unlist(sapply(wthn_btwn$pcr_intradist, function(y) {
        names(which(y > thresh_pcr & y == max(y)))
      })))
      good_pcrs[to_flag] <- FALSE

      # spot samples having only one representative (singelton)
      replicates_table <- as.data.frame.matrix(table(replicates0, good_pcrs))
      singleton_sample <- rownames(replicates_table)[replicates_table$`TRUE`<2]

      # update good_pcrs
      good_pcrs[replicates0 %in% singleton_sample] <- FALSE

      replicates <- replicates[which(good_pcrs==T)]

      # stop iterations when no more replicate spotted
      if(sum(good_pcrs==F) == nb_bad_pcr) {
        break
      }
    }
    return(names(which(good_pcrs==F)))
  }
}


#' @describeIn pcrslayer Computes a list of dissimilarities in OTU composition within a biological sample \emph{dw} and between biological samples \emph{db}.
#' @export pcr_within_between

pcr_within_between <- function(metabarlist, replicates) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads_table <- metabarlist$reads

    if (!replicates %in% colnames(metabarlist$pcrs)) {
      stop("replicates should be in colnames of the metabarlist$pcr table")
    }

    replicates <- metabarlist$pcrs[,replicates]

    out = pcr_within_between_internal(reads_table, replicates)
    return(out)
  }
}


#internal function pcr_within_between
pcr_within_between_internal <- function(reads, replicates) {
  # barycentre calculation and intradist function
  # data standardization
  reads_stdt <- reads/rowSums(reads)

  bar <- rowsum(reads_stdt, replicates)/as.vector(table(replicates))

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

#' @describeIn pcrslayer Vizualize \emph{dw} and \emph{db} dissimilarities and the threshold (defined automatically) above which pcr replicates are considered as too dissimilar.
#' @export check_pcr_thresh

check_pcr_thresh <- function(wthn.btwn, thresh.method = "intersect") {

  if(is.list(wthn.btwn) == F | all(names(wthn.btwn)==c("bar_dist", "pcr_intradist")) == F | length(wthn.btwn)!=2) {
    stop("wthn.btwn should be a list of length 2 and of names bar_dist and pcr_intradist,
         typically an ouput from pcr_within_between")
  }

  if (thresh.method != "intersect" & thresh.method != "mode") {
    stop('thresh.method should be one of "intersect" or "mode"')
  }

  d.bar <- density(wthn.btwn$bar_dist)
  d.intra <- density(unlist(wthn.btwn$pcr_intradist))

  d.out <- rbind(data.frame(d.bar[c("x","y")], distance="between samples"),
                data.frame(d.intra[c("x","y")], distance="between pcr replicates"))

  thresh.pcr <- pcr_threshold_estimate(wthn.btwn, thresh.method)

  p =
    ggplot(d.out, aes(x=x, y=y, color=distance)) +
      geom_line() + labs(x="distance", y="density", color="Distances") +
      theme_bw()
  if(is.null(thresh.pcr)) {p} else {p + geom_vline(xintercept = thresh.pcr, size=0.3, lty=2)}
}

#internal function pcr_threshold_estimate
pcr_threshold_estimate <- function(wthn.btwn, thresh.method = "intersect") {
  dinter.max <- max(wthn.btwn$bar_dist)
  ddinter <- density(wthn.btwn$bar_dist, from = 0, to = 1)
  ddintra <- density(unlist(wthn.btwn$pcr_intradist), from = 0, to = 1)

  # assumption that each of them have a "unimodal" distribution
  dintra.mode <- ddintra$x[which.max(ddintra$y)]
  dinter.mode <- ddinter$x[which.max(ddinter$y)]

  if (thresh.method == "intersect") {
    p <- which(ddintra$y - ddinter$y > 0 & ddinter$y > dintra.mode)
    out <- ddinter$x[p[length(p)]]
  } else {
    out <- dinter.mode
  }
  return(out)
}


#' @describeIn pcrslayer Vizualize pcrs dissimilarity patterns and pcr replicates centroids.
#' @export check_pcr_repl

check_pcr_repl <- function(metabarlist, replicates, groups = NULL, dyspcr = NULL) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads <- metabarlist$reads

    if (!replicates %in% colnames(metabarlist$pcrs)) {
      stop("replicates should be in colnames of the metabarlist$pcrs table")
    }
    replicates <- metabarlist$pcrs[,replicates]

    if (!is.null(groups)) {
      if(!groups %in% colnames(metabarlist$pcrs)) {
        stop("groups should be in colnames of the metabarlist$pcrs table")
      } else {
        groups <- metabarlist$pcrs[,groups]
      }}

    if (!is.null(dyspcr)) {
          dyspcr <- ifelse(rownames(reads) %in% dyspcr, "dyspcr", "0ok")
      }

    reads_stdt <- reads/rowSums(reads)
    bar <- rowsum(reads_stdt, replicates)/as.vector(table(replicates))
    all <- rbind(reads_stdt, bar)
    all.d <- vegdist(all, "bray")

    mds <- cmdscale(all.d, k = 2, eig = T, add = T)
    d <- data.frame(mds$points)
    d$points <- ifelse(rownames(d) %in% rownames(bar), "bary", "samp")

    d.new <- data.frame(d[d$points=="samp",])
    d.new$replicates <- replicates[match(rownames(d.new), rownames(reads))]
    d.new$bary_x <- d$X1[match(d.new$replicates, rownames(d))]
    d.new$bary_y <- d$X2[match(d.new$replicates, rownames(d))]
    d.new$groups <- 0
    d.new$dyspcr <- "0ok"

    if(!is.null(groups)) {d.new$groups = groups}
    if(!is.null(dyspcr)) {d.new$dyspcr = dyspcr}

    ggplot(d.new, aes(x=X1, y=X2, color=groups)) +
      geom_point(aes(shape = dyspcr)) + theme_bw() +
      scale_shape_manual(values = c(19,8), labels = c("good", "bad")) +
      geom_segment(aes(x=bary_x, y=bary_y, xend=X1, yend=X2), color="grey") +
      labs(x=paste("PCoA1 (", round(100*mds$eig[1]/sum(mds$eig),2), "%)", sep=""),
           y=paste("PCoA2 (", round(100*mds$eig[2]/sum(mds$eig),2), "%)", sep=""),
           shape="PCR type")
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
