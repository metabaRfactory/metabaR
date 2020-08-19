#' Plotting a taxonomic tree
#'
#' Plots the taxonomic tree from a \code{metabarlist} object and maps an attribute onto it
#'
#'
#' @param metabarlist    a \code{metabarlist} object.
#'                       Should contain taxonomic information in the `motus` table.
#' @param taxo           a character string or vector of strings indicating the name of
#'                       the column (or group of columns) containing the full taxonomic
#'                       information in the `motus` table from the \code{metabarlist}
#'                       object.
#' @param sep.level      an optional character string to separate the terms.
#'                       Required only if `taxo` is a string. \code{NA} not allowed.
#' @param sep.info       an optional character string to separate the terms.
#'                       Required only if `taxo` is a string. \code{NA} not allowed.
#' @param thresh         a numeric indicating the relative abundance below which
#'                       taxon labels won't be plotted
#'
#' @name ggtaxplot
#'
#' @return a ggplot
#'
#' @details
#' This function allows users to visualise the full taxonomic tree for a set of samples and to map attributes onto the produced tree (e.g. number of reads per node/branches, nb of MOTUs, etc.). The taxonomic information should follow a standard structure across samples (e.g. standard taxonomy as in Genbank, SILVA or BOLD or with defined taxonomic levels if `taxo` is a vector) by decreasing level of taxonomic resolution: the function does not infer missing taxonomic ranks. The taxonomic information should contain a level that is common to all MOTU taxonomic assignments (common ancestor).
#'
#' @seealso \code{\link{taxoparser}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # Using a taxonomic path
#'
#' ## on all data. can take a while
#' ggtaxplot(soil_euk, "path", sep.level = ":", sep.info = "@")
#'
#' ## show only taxonomic labels if taxon has a relative abundance > 1e-3
#' ggtaxplot(soil_euk, "path", sep.level = ":", sep.info = "@", thresh = 1e-3)
#'
#' ## run on a particular clade e.g. here arthropoda, otherwise difficult to read
#' arthropoda <- subset_metabarlist(soil_euk,
#'   table = "motus",
#'   indices = grepl("Arthropoda", soil_euk$motus$path)
#' )
#' ### plot
#' ggtaxplot(arthropoda, "path", sep.level = ":", sep.info = "@")
#'
#'
#' # Using a taxonomic table
#' taxo.col <- c(
#'   "phylum_name", "class_name", "order_name",
#'   "family_name", "genus_name", "species_name"
#' )
#' ggtaxplot(arthropoda, taxo.col)
#' ## to do the same on the full soil_euk data,
#' # add a column "kingdom_name" with value "Eukaryota" in soil_euk$motus
#'
#' @author Lucie Zinger
#' @import ggplot2
#' @importFrom igraph graph V layout.reingold.tilford get.vertex.attribute get.data.frame
#' @export ggtaxplot

ggtaxplot <- function(metabarlist, taxo, sep.level, sep.info, thresh = NULL) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (!is.character(taxo)) {
      stop("`taxo` should be a character string or vector")
    }

    if (any(!taxo %in% colnames(metabarlist$motus))) {
      stop("`taxo` should be a character string or vector of strings present
           in the column names of the table `motus`")
    }

    if (length(taxo) == 1) {

      parse <- unname(taxoparser(metabarlist$motus[, taxo], sep.level, sep.info))
      path <- sapply(parse, toString)

    } else {
      path <- sapply(
        1:nrow(metabarlist$motus),
        function(x) toString(metabarlist$motus[x, taxo][!is.na(metabarlist$motus[x, taxo])])
      )
      # above trick is necessary to deal with intra NAs
      parse <- strsplit(path, ", ")
    }
    parse.mat <- do.call(rbind, lapply(parse, `length<-`, max(lengths(parse))))

    if (length(unique(parse.mat[,1]))>1) {
      stop("No common ancestor found: tree cannot be plotted")
    }

    edgelist <- NULL
    for (i in rev(2:ncol(parse.mat))) {
      idx <- which(!is.na(parse.mat[, i]))
      kid <- parse.mat[idx, i]
      parent <- parse.mat[idx, (i - 1)]
      kidfull <- apply(parse.mat[idx, 1:i, drop = F], 1, toString)
      parentfull <- apply(parse.mat[idx, 1:(i - 1), drop = F], 1, toString)

      edgelist <- rbind(
        edgelist,
        unique(cbind(
          parentfull,
          kidfull,
          parent, kid
        ))
      )
    }

    # edgelist[is.na(edgelist)] = "NA"
    # edgelist = edgelist[!duplicated(edgelist),] #should all be NAs.

    g <- igraph::graph.edgelist(edgelist[rev(1:nrow(edgelist)),
                                         c("parentfull", "kidfull")], directed = F)

    # Rename nodes with kid names
    igraph::V(g)$name2 <- ifelse(igraph::V(g)$name %in% edgelist[, "parent"],
      igraph::V(g)$name,
      edgelist[match(igraph::V(g)$name, edgelist[, "kidfull"]), "kid"]
    )


    # Assign proportion of MOTUs and reads to each node
    igraph::V(g)$motus <-
      sapply(V(g)$name, function(x) sum(grepl(x, path, fixed = T))) / nrow(metabarlist$motus)
    metabarlist$motus$count <- colSums(metabarlist$reads)
    igraph::V(g)$reads <-
      sapply(V(g)$name, function(x) sum(metabarlist$motus$count[grep(x, path, fixed = T)])) /
        sum(metabarlist$reads)

    # plots

    # nodes
    coords <- layout.reingold.tilford(g, root = 1, circular = F)
    colnames(coords) <- c("x", "y")
    vdf <- data.frame(as.data.frame(get.vertex.attribute(g)), coords)
    vdf <- vdf[which(vdf$motus != 0), ]

    # add segments
    edf <- get.data.frame(g)

    edf$from.x <- vdf$x[match(edf$from, as.vector(vdf$name))]
    edf$from.y <- vdf$y[match(edf$from, as.vector(vdf$name))]
    edf$to.x <- vdf$x[match(edf$to, as.vector(vdf$name))]
    edf$to.y <- vdf$y[match(edf$to, as.vector(vdf$name))]

    gp <-
      ggplot(data = vdf, aes(x = .data$x, y = .data$y,
                             size = .data$motus * 100, colour = .data$reads * 100)) +
      geom_segment(
        data = edf,
        aes(
          x = .data$from.x, xend = .data$to.x,
          y = .data$from.y, yend = .data$to.y
        ), size = 0.2, colour = "grey"
      ) +
      geom_point() +
      scale_color_viridis_c() +
      theme_void() +
      theme(legend.position = "bottom", legend.direction = "horizontal") +
      labs(color = "%reads", size = "%motus")

    if (is.null(thresh)) {
      gp <- gp + geom_text(aes(label = .data$name2), color = "darkgrey", show.legend = FALSE)
    } else {
      gp <- gp + geom_text(
        data = vdf[which(vdf$reads > thresh), ], aes(label = .data$name2),
        color = "darkgrey", show.legend = FALSE
      )
    }

    return(gp)
  }
}
