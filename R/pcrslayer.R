#' Detecting PCRs replicate outliers.
#'
#' Detecting failed PCRs, i.e. PCR replicate outliers based on PCR similarity in composition of MOTUs.
#'
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param replicates    a vector corresponding to the sample names to which
#'                      pcr replicates belongs. Default is the `sample_id` column of the `pcrs` table
#'                      from the \code{metabarlist} object.
#' @param method        a character indicating which method should be used to identify PCR outliers. Can be
#'                      \code{"centroid"} or \code{"pairwise"}. Default is \code{"centroid"}.
#' @param FUN           a function for computing distances between replicates.
#'                      Default is Bray-Curtis distances on MOTU relative abundance table.
#' @param output_col    a character string for the column name in the `pcrs` table, in
#'                      which the result will be stored.
#'                      Default is "functional_pcr"
#' @param thresh.method a character indicating which method should be used to define the filtering
#'                      threshold. Can be \code{"intersect"} or \code{"mode"}.
#' @param plot          a boolean vector indicating whether dissimilarity distribution should be plotted.
#'                      Default is \code{TRUE}
#' @param wthn.btwn     an ouput from \code{pcr_within_between}
#' @param groups        a column name in the `pcrs` table corresponding to a factor giving the groups
#'                      from which the graphical colorscheme is drawn.
#' @param funcpcr       a boolean vector indicating whether PCRs are functional (\code{TRUE}) or not.
#' @param reads_table   a table `reads` from a \code{metabarlist} object
#'
#' @details
#'
#' The \code{pcrslayer} function identifies potential non-functional PCR reactions based on their reproducibility. It compares the dissimilarities in MOTU composition within a biological sample (i.e. between PCR replicates, hereafter \emph{dw}) vs. between biological samples (hereafter \emph{db}). It relies on the assumption that PCR replicates from the same biological sample should be more similar than two different biological samples (\emph{dw} < \emph{db}). Two methods for computing \emph{dw} and \emph{db} are available.
#'
#'\itemize{
#' \item{With method \code{"centroid"}, both \emph{dw} and \emph{db} distances are based on the centroid of samples. More specifically, a centroid community of each sample is built by computing the average MOTU abundances of the sample's pcr replicates. Then \emph{dw} is defined as the distance between pcr replicates and their corresponding centroid, while \emph{db} is defined as the distances between centroids of different samples. A PCR replicate having a \emph{dw} above a given dissimilarity threshold \emph{tresh} is considered as an outlier, i.e. too distant from its associated average MOTU community (\code{method="centroid"}). If only one single PCR replicate is representative of a biological sample after this trimming, it is also considered as a failed PCR. The process is repeated iteratively to recompute the sample centroid, as well as \emph{dw} and \emph{db} until no more PCRs are excluded from the analysis.}
#' \item{With method \code{"pairwise"}, pairwise distances between all PCR replicates are computed and then classified into \emph{dw} or \emph{db} depending on the pair considered. A PCR replicate having a an average \emph{dw} above a given dissimilarity threshold \emph{tresh} is considered as an outlier, i.e. too distant from its associated replicates (\code{method="pairwise"}). If only one single PCR replicate is representative of a biological sample after this trimming, it is also considered as a failed PCR. In this case, no iterations are done.}
#' }
#' For both methods, the user is free to chose their own distance metric and whether distances should be computed on relative abundances or true abundances, through the argument \code{FUN}. Two methods are currently pre-encoded:
#' \itemize{
#' \item{FUN_pcrdist_bray_freq computes Bray-Curtis distances on MOTU relative abundances}
#' \item{FUN_pcrdist_coa_freq computes Euclidean distances between PCRs from a Correspondance Analysis on MOTU relative abundances}
#' }.
#'
#' The \code{pcr_within_between} function is part of \code{pcrslayer}, and computes dissimilarities in MOTU composition within a biological sample \emph{dw} and between biological samples \emph{db} following either \code{method="centroid"} or \code{method="pairwise"} methods.
#'
#' The threshold \emph{tresh} is defined automatically with two alternative methods. Either it is the intersection of \emph{dw} and \emph{db} distributions (\code{tresh.method="intersect"}). Alternatively, it is the mode of the \emph{db} distribution (\code{tresh.method="mode"}).
#'
#' The \code{check_pcr_thresh} function enables visualisation of \emph{dw} and \emph{db} distributions. Function \code{check_pcr_repl} enables visualization of PCR replicate dissimilarity patterns in a NMDS ordination and distance from their average OTU community.
#'
#' The \code{check_pcr_repl} function enables visualisation of dissimilarity patterns across all PCRs while showing PCR replicate centroids through a Principal Coordinate Analysis (PCoA) based on Bray-Curtis dissimilarities.
#'
#' @return
#'
#' The \code{pcrslayer} function returns a metabarlist with a new boolean column vector of name
#' `output_col` in the `pcrs` table indicating whether the PCR is functional (\code{TRUE}) or
#'         or failed (\code{FALSE}).
#'
#' The \code{pcr_within_between} function returns a list of \emph{dw} and \emph{db} dissimilarities.
#'
#' The \code{check_pcr_thresh} and \code{check_pcr_repl} functions return \code{ggplot} objects.
#'
#' @seealso \code{\link{tagjumpslayer}}, \code{\link{contaslayer}} for other data curation procedures.
#'
#' @examples
#'
#'\donttest{
#'
#' library(ggplot2)
#' data(soil_euk)
#'
#' ## Consider only biological samples with # reads > 0
#' soil_euk_sub <- subset_metabarlist(soil_euk,
#'                                    "pcrs",
#'                                    soil_euk$pcrs$type == "sample" &
#'                                    rowSums(soil_euk$reads>0))
#'
#' ## Visualisation of within vs. between sample dissimilarities
#' soil_euk_sub_wb <- pcr_within_between(soil_euk_sub)
#' check_pcr_thresh(soil_euk_sub_wb, thresh.method = "intersect")
#'
#' ## Visualisation of replicates through PCoA
#' # create grouping factor according to habitat and material
#' soil_euk_sub$pcrs$habitat_material <- soil_euk_sub$pcrs$sample_id
#' idx <- match(levels(soil_euk_sub$pcrs$habitat_material), rownames(soil_euk_sub$samples))
#' levels(soil_euk_sub$pcrs$habitat_material) <- paste(soil_euk_sub$samples$Habitat[idx],
#'                                                     soil_euk_sub$samples$Material[idx],
#'                                                     sep = " | ")
#' # visualise dissimilarity patterns
#' mds <- check_pcr_repl(soil_euk_sub, groups = soil_euk_sub$pcrs$habitat_material)
#' mds + labs(color = "sample type")
#'
#' # identify failed PCRs
#' soil_euk_sub2 <- pcrslayer(soil_euk_sub)
#' tail(colnames(soil_euk_sub2$pcrs))
#' length(which(soil_euk_sub2$pcrs$functional_pcr == FALSE))
#'
#' # define a color vector that corresponds to the different habitat types. Should be named.
#' mds <- check_pcr_repl(soil_euk_sub2,
#'                       groups = soil_euk_sub2$pcrs$habitat_material,
#'                       funcpcr = soil_euk_sub2$pcrs$functional_pcr)
#' mds + guides(color = FALSE)
#'
#' }
#'
#' @author Lucie Zinger, Clement Lionnet, Fred Boyer
#' @importFrom vegan vegdist cca
#' @importFrom stats dist
#' @describeIn pcrslayer Detect failed PCRs, i.e. PCR outliers in a \code{metabarlist} object.
#' @export pcrslayer

pcrslayer <- function(metabarlist,
                      replicates = NULL,
                      method = "centroid",
                      FUN = FUN_pcrdist_bray_freq,
                      thresh.method = "intersect",
                      output_col = "functional_pcr",
                      plot = T) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {

    metabarlist0 <- metabarlist
    reads_table0 <- metabarlist$reads

    if (is.null(replicates)) {
      replicates <- metabarlist$pcrs$sample_id
    }

    if (length(replicates)!=nrow(metabarlist$reads)) {
      stop("`replicates` length should be equal to the number of rows of table `reads`")
    }

    replicates0 <- replicates

    if (!method  %in% c("centroid", "pairwise")) {
      stop("method should be one of `centroid` or `pairwise`")
    }

    if (!thresh.method %in% c("intersect", "mode")) {
      stop("thresh.method should be one of `intersect` or `mode`")
    }

    #vector of pcrs boolean: good = T bad = F
    good_pcrs <- rep(T, length(replicates0))
    names(good_pcrs) <- rownames(metabarlist$pcrs)

    #tag empty pcrs
    idx <- which(rowSums(reads_table0) == 0)
    good_pcrs[idx] <- FALSE

    if(method=="centroid") {

    iteration <- 0
    repeat{
      iteration <- iteration + 1
      print(paste("Iteration", iteration))

      #get good pcrs
      #reads_table <- reads_table0[names(good_pcrs)[good_pcrs==T], ]
      metabarlist <- subset_metabarlist(metabarlist0, "pcr", good_pcrs)

      replicates <- replicates0[which(good_pcrs==T)]

      nb_bad_pcr <- sum(good_pcrs==F)
      wthn_btwn <- pcr_within_between(metabarlist,
                                      replicates = replicates,
                                      method = method,
                                      FUN = FUN)
      thresh_pcr <- pcr_threshold_estimate(wthn_btwn, thresh.method)
      if (plot == T) {
        p = check_pcr_thresh(wthn_btwn, thresh.method)
        print(p + ggtitle(paste("Check PCR threshold - Iteration #", iteration)))
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
    } else if(method=="pairwise") {

      metabarlist <- metabarlist0
      replicates <- replicates0

      wthn_btwn <-pcr_within_between(metabarlist,
                                     replicates = replicates,
                                     method = method,
                                     FUN = FUN)
      thresh_pcr <- pcr_threshold_estimate(wthn_btwn, thresh.method)
      if (plot == T) {
        p = check_pcr_thresh(wthn_btwn, thresh.method)
        print(p + ggtitle(paste("Check PCR threshold - Iteration #", iteration)))
      }

      ### pb names of comparisons
      to_flag <- unname(unlist(sapply(wthn_btwn$pcr_intradist, function(y) {
        names(which(y > thresh_pcr & y == max(y)))
      })))
      good_pcrs[to_flag] <- FALSE

      #check for singeltons
      replicates_table <- as.data.frame.matrix(table(replicates0, good_pcrs))
      singleton_sample <- rownames(replicates_table)[replicates_table$`TRUE`<2]

      # update good_pcrs
      good_pcrs[replicates0 %in% singleton_sample] <- FALSE

    }

    #### warning if more than 20% of replicates are removed
    if (sum(good_pcrs==F) / length(replicates0) > 0.2) {
      warning("More than 20% of pcr replicates are removed !")
    }

    metabarlist0$pcrs[,output_col] <- good_pcrs
    return(metabarlist0)

    }
}


#' @describeIn pcrslayer Computes a list of dissimilarities in OTU composition within a biological sample \emph{dw} and between biological samples \emph{db}.
#' @importFrom utils combn
#' @importFrom stats xtabs median
#' @export pcr_within_between

pcr_within_between <- function(metabarlist,
                               replicates = NULL,
                               FUN = FUN_pcrdist_bray_freq,
                               method="centroid") {

  if (suppressWarnings(check_metabarlist(metabarlist))) {

    reads_table <- metabarlist$reads

  if(any(rowSums(reads_table)==0 )) {
      stop("you have empty rows in table `reads`")}

    if (is.null(replicates)) {
      replicates <- metabarlist$pcrs$sample_id
    }

    if (length(replicates)!=nrow(metabarlist$reads)) {
      stop("`replicates` length should be equal to the number of rows of table `reads`")
    }

    if (!method  %in% c("centroid", "pairwise")) {
      stop('method should be one of "centroid" or "pairwise"')
    }


    if(method=="centroid") {

      reads_stdt <- reads_table/rowSums(reads_table)
      bar <- rowsum(reads_stdt, replicates)/as.vector(table(replicates))

      # between barycentre distances

      if(substitute(FUN)=="FUN_pcrdist_coa_freq") {
        stop("FUN cannot be `FUN_pcrdist_coa_freq` if method is `centroid`")
      }

      bar_dist <- FUN(bar)

      # within replicates distances
      pcr_intradist <- lapply(1:nrow(bar), function(x) {
       ind <- which(replicates == rownames(bar)[x])
       sapply(ind, function(y) {
         out <- FUN(rbind(bar[x, ], reads_stdt[y, ]))
         names(out) <- rownames(reads_stdt)[y]
         out
         })
       })
      names(pcr_intradist) <- rownames(bar)

    } else if (method=="pairwise"){

      all.dist <- as.data.frame(t(combn(rownames(reads_table),2)))
      all.dist$dist <- FUN(reads_table)
      all.dist$S1 <- replicates[match(all.dist$V1, rownames(reads_table))]
      all.dist$S2 <- replicates[match(all.dist$V2, rownames(reads_table))]
      all.dist$type <- ifelse(as.vector(all.dist$S1)==as.vector(all.dist$S2), "intra", "inter")

      bar_dist = all.dist$dist[all.dist$type=="inter"]
      tmp = all.dist[all.dist$type=="intra",]
      pcr_intradist <- lapply(unique(tmp$S1), function(x) {
        df <- droplevels(tmp[tmp$S1==x,])
        df$V1 <- factor(df$V1, levels = (unique(c(levels(as.factor(df$V1)), levels(as.factor(df$V2))))))
        df$V2 <- factor(df$V2, levels = (unique(c(levels(as.factor(df$V1)), levels(as.factor(df$V2))))))
        out <- as.matrix(xtabs(dist ~ V2 + V1, data=df))
        out[upper.tri(out)] <- out[lower.tri(out)]
        colMeans(out)
      })
      names(pcr_intradist) <- unique(tmp$S1)
    }

    return(list(bar_dist = bar_dist, pcr_intradist = pcr_intradist))

  }
}

#' @describeIn pcrslayer Computes Bray-Curtis distances between pcrs standardised by the total number of reads.
#' @export FUN_pcrdist_bray_freq

FUN_pcrdist_bray_freq <- function(reads_table) {
  reads_std <- (reads_table)/rowSums(reads_table)
  distance_matrix <- vegdist(reads_std, method = "bray")
  return(distance_matrix)
}

#' @describeIn pcrslayer Computes eclidean distances from a PCoA ordination between pcrs standardised by the total number of reads.
#' @export FUN_pcrdist_coa_freq

FUN_pcrdist_coa_freq <- function(reads_table) {
  reads_std <- cca((reads_table)/rowSums(reads_table))
  distance_matrix <- dist(scores(reads_std)$sites)
  return(distance_matrix)
}


#' @describeIn pcrslayer Vizualize \emph{dw} and \emph{db} dissimilarities and the threshold (defined automatically) above which pcr replicates are considered as too dissimilar.
#' @importFrom stats density
#' @export check_pcr_thresh

check_pcr_thresh <- function(wthn.btwn, thresh.method = "intersect") {

  if(is.list(wthn.btwn) == F | all(names(wthn.btwn)==c("bar_dist", "pcr_intradist")) == F | length(wthn.btwn)!=2) {
    stop("wthn.btwn should be a list of length 2 and of names bar_dist and pcr_intradist,
         typically an ouput from pcr_within_between")
  }


  if (!thresh.method %in% c("intersect", "mode")) {
    stop('thresh.method should be one of "intersect" or "mode"')
  }

  d.bar <- density(wthn.btwn$bar_dist)
  d.intra <- density(unlist(wthn.btwn$pcr_intradist))

  d.out <- rbind(data.frame(d.bar[c("x","y")], distance="between samples"),
                data.frame(d.intra[c("x","y")], distance="within samples"))

  thresh.pcr <- pcr_threshold_estimate(wthn.btwn, thresh.method)

  p =
    ggplot(d.out, aes(x=.data$x, y=.data$y, color=.data$distance)) +
      geom_line() + labs(x="distance", y="density", color="Distances") +
      theme_bw()
  if(is.null(thresh.pcr)) {p} else {p + geom_vline(xintercept = thresh.pcr, size=0.3, lty=2)}
}

# function pcr_threshold_estimate
#' @describeIn pcrslayer Estimates the cutoff above which a pcr is considered as an outlier.
#' @export pcr_threshold_estimate
pcr_threshold_estimate <- function(wthn.btwn, thresh.method = "intersect") {
  dinter.max <- max(wthn.btwn$bar_dist)
  ddinter <- density(wthn.btwn$bar_dist, from = 0, to = 1)
  ddintra <- density(unlist(wthn.btwn$pcr_intradist), from = 0, to = 1)

  # assumption that each of them have a "unimodal" distribution
  dintra.mode <- ddintra$x[which.max(ddintra$y)]
  dinter.mode <- ddinter$x[which.max(ddinter$y)]

  if (thresh.method == "intersect") {
    p <- which(ddintra$y - ddinter$y > 0 & ddinter$x > dintra.mode)
    out <- ddinter$x[p[length(p)]]
  } else {
    out <- dinter.mode
  }
  return(out)
}


#' @describeIn pcrslayer Vizualize pcrs dissimilarity patterns and pcr replicates centroids.
#' @importFrom stats cmdscale
#' @export check_pcr_repl

check_pcr_repl <- function(metabarlist,
                           replicates=NULL,
                           groups = NULL,
                           funcpcr = NULL) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    reads <- metabarlist$reads

    if (is.null(replicates)) {
      replicates <- metabarlist$pcrs$sample_id
    }

    if (length(replicates)!=nrow(metabarlist$reads)) {
      stop("`replicates` length should be equal to the number of rows of table `reads`")
    }

    if (is.null(groups)) {
      stop("groups should be defined")
    }

    if (length(groups)!=nrow(metabarlist$reads)) {
      stop("`groups` length should be equal to the number of rows of table `reads`")
    }

    if (!is.null(funcpcr)) {

      if (length(funcpcr)!=nrow(metabarlist$reads)) {
        stop("`funcpcr` length should be equal to the number of rows of table `reads`")
      }

      if (!is.logical(funcpcr)) {
        stop("`funcpcr` should be logical")
      }

      funcpcr <- ifelse(funcpcr, "0ok", "dyspcr")

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
    d.new$funcpcr <- "0ok"

    if(!is.null(groups)) {d.new$groups = groups}
    if(!is.null(funcpcr)) {d.new$funcpcr = funcpcr}

    ggplot(d.new, aes(x=.data$X1, y=.data$X2, color=.data$groups)) +
      geom_point(aes(shape = funcpcr)) + theme_bw() +
      scale_shape_manual(values = c(19,8), labels = c("good", "bad")) +
      geom_segment(aes(x=.data$bary_x, y=.data$bary_y, xend=.data$X1, yend=.data$X2),
                   color="grey") +
      labs(x=paste("PCoA1 (", round(100*mds$eig[1]/sum(mds$eig),2), "%)", sep=""),
           y=paste("PCoA2 (", round(100*mds$eig[2]/sum(mds$eig),2), "%)", sep=""),
           shape="PCR type")
  }
}
