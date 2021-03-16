#' Including SILVAngs pipeline taxonomic annotations
#'
#' Importing and formatting taxonomic annotations obtained with the SILVAngs pipeline (\href{https://ngs.arb-silva.de/silvangs/}{https://ngs.arb-silva.de/silvangs/}) and including it in a \code{metabarlist} object.
#'
#'
#' @param metabarlist  a \code{metabarlist} object
#' @param silva.path   path to a table from the SILVAngs pipeline,
#'                     typically zipfile > ssu > exports > xxx---ssu---otus.csv
#' @param clust.path   path to a file from the SilvAngs pipeline indicating MOTU cluster membership,
#'                     typically zipfile > ssu > stats > sequence_cluster_map > data >
#'                     xxx---ssu---sequence_cluster_map---tmptaxo.clstr
#' @param taxonomy.path   path to a file containing the SILVA taxonomy. See below for the url
#' @name silva_annotator
#'
#' @return a \code{metabarlist} object with table `motus` including the taxonomic assignments from silva
#'
#' @details
#'
#' Users can be interested in using the SILVAngs pipeline for assigning a taxon to DNA sequences/MOTUs. Assuming that it is done on data already filtered (i.e. dereplicated, clustering etc.), resulting in one sequence per sequence or MOTU, one can then use the SILVAngs pipeline by setting all filtering parameters to "null" (i.e. returning to no filtration) and the taxonomic assignments parameters by default (or following the user's preferences). These taxonomic assignments are compiled in the data archives provided by SILVAngs, in which two files are important for the function `silva_annotator`:
#' \itemize{
#' \item{}{`zipfile>ssu>exports>xxx---ssu---otus.csv`: a csv file containing the taxonomic assignment for each OTU}
#' \item{}{`zipfile>ssu>stats>sequence_cluster_map>data>xxx---ssu---sequence_cluster_map---tmptaxo.clstr`: a file containing the mapping informations of sequences to their respective MOTUs. Even though the clustering is done before using the SILVAngs pipeline and that the clustering parameters in this pipeline are settled to 100\% identity of clustering, SILVAngs uses CDHit, which can group together OTUs sharing the same prefix/suffix. So it is common to retrieve less assignations than what was expected.}
#' }
#'
#' In addition, the function requires the SILVA taxonomy. It is available in the \code{metabaR} companion R package; https://github.com/metabaRfactory/metabaR_external_data.
#'
#' @seealso \code{\link{taxodecider}}
#'
#' @examples
#'
#' \donttest{
#'
#' dir <- tempdir()
#' url <- "https://raw.githubusercontent.com/metabaRfactory/metabaR_external_data/master/"
#'
#' silva_file <- "lit_euk---ssu---otus.csv"
#' silva_url <- paste(url, silva_file, sep="")
#' silva_path <- file.path(dir, silva_file)
#' download.file(silva_url, silva_path)
#'
#' clust_file <- "lit_euk---ssu---sequence_cluster_map---litiere_euk_cl97_agg_filt.clstr"
#' clust_url <- paste(url, clust_file, sep="")
#' clust_path <- file.path(dir, clust_file)
#' download.file(clust_url, clust_path)
#'
#' taxonomy_file <- "tax_slv_ssu_138.1.txt"
#' taxonomy_url <- paste(url, taxonomy_file, sep="")
#' taxonomy_path <- file.path(dir, taxonomy_file)
#' download.file(taxonomy_url, taxonomy_path)
#'
#' data(soil_euk)
#' soil_euk <- silva_annotator(
#'    metabarlist = soil_euk,
#'    silva.path = silva_path,
#'    clust.path = clust_path,
#'    taxonomy.path = taxonomy_path)
#'
#' library(ggplot2)
#'
#' ggplot(soil_euk$motus, aes(x=factor(1), fill=phylum_silva)) +
#'   geom_bar() + coord_polar("y") +
#'   theme_minimal() + labs(x=NULL, y=NULL) +
#'   theme(legend.position = "bottom")
#'
#' }
#'
#' @author Lucie Zinger, Anne-Sophie Benoiston
#' @importFrom seqinr read.fasta
#' @importFrom utils read.csv
#' @export silva_annotator

silva_annotator <- function(metabarlist, silva.path, clust.path, taxonomy.path) {

  if (suppressWarnings(check_metabarlist(metabarlist))) {
    motus <- metabarlist$motus

    silva <- read.csv(silva.path, header = T, sep = "\t", skip = 1, row.names = NULL)
    colnames(silva) <- c(colnames(silva)[-c(1, ncol(silva))], "classif_ncbi", "classif_silva")

    # taxo results formating
    tmp <- sapply(strsplit(as.vector(silva$classif_silva), "\\|"), function(x) x[4])


    if (dim(table(is.na(tmp))) > 1) {
      tmp[is.na(tmp)] <- "No blast hit" # format according to qiime
    }

    tmp.uniq <- unique(tmp)


    ### set taxonomy stuff
    # dictionnary of taxo ranks
    # ncbi-like
    taxolev.dict1 <- c(
      "superkingdom","superkingdom2", "superkingdom3", "kingdom", "subkingdom",
      "superphylum", "phylum", "subphylum", "infraphylum",
      "superclass", "class", "subclass", "infraclass",
      "superorder", "order", "cohort", "suborder",
      "subcohort", "infraorder", "parvorder",
      "superfamily", "family", "subfamily",
      "tribe", "subtribe",
      "genus", "subgenus", "section", "subsection", "series",
      "species group", "species subgroup", "species", "subspecies",
      "varietas", "forma",
      "no rank"
    )
    #load old version (incomplete)

    taxonomy <- read.csv(taxonomy.path, header = F, sep = "\t", row.names = NULL)
    taxonomy$V3 <- gsub("superkingdom", "superkingdom2", taxonomy$V3, fixed=T)
    taxonomy$V3 <- gsub("domain", "superkingdom", taxonomy$V3, fixed=T)
    taxonomy$V3 <- gsub("major_clade", "superkingdom3", taxonomy$V3, fixed=T)
    taxonomy$match <- sapply(strsplit(taxonomy$V1, ";"), function(x) x[length(x)])
    taxonomy <- taxonomy[taxonomy$match != "uncultured",]
    #Manual editing to get Metazoa, SARs, etc. info correctly
    taxonomy$V3[grep("Metazoa;$", taxonomy$V1)] <- "kingdom"
    taxonomy$V3[grep("Animalia;$", taxonomy$V1)] <- "superkingdom3"
    taxonomy$V3[grep("Archaeplastida;$", taxonomy$V1)] <- "superkingdom2" # plants
    taxonomy$V3[grep("Eukaryota;SAR;$", taxonomy$V1)] <- "superkingdom2" # part of protists

    tmp.uniq <- gsub("uncultured;$", "", tmp.uniq, perl=T)
    tmp.test <- sapply(strsplit(tmp.uniq, ";"), function(x) x[length(x)])

    if (any(gsub("No blast hit", "Bacteria", tmp.test, perl=T) %in% taxonomy$match == FALSE)) {
      stop(paste(
        "You probably use the wrong version of the SILVA taxonomy. The following taxa cannot be found:\n",
        toString(tmp.test[!gsub("No blast hit", "Bacteria", tmp.test, perl=T) %in% taxonomy$match])
      ),
      "\n\nPlease use the correct SILVA version, or rerun the SILVAngs pipeline on your data with the current SILVA version")
    }

    tmp1 <- do.call("rbind", lapply(strsplit(tmp.uniq, ";"), function(x) {
      if("SAR" %in% x){
        if("Bacteria" %in% x) {
          idx <- grep("Bacteria", taxonomy$V1, perl = T)
        } else {
          idx <- grep("Eukaryota", taxonomy$V1, perl = T)
        }
        taxonomy.sub <- taxonomy[idx,]
        taxonomy.sub$match[grep("Eukaryota;SAR;$", taxonomy.sub$V1)] <- "SARx"
        names(x) <- taxonomy.sub$V3[match(x ,taxonomy.sub$match)]
        if(x[which(x=="SAR")-1]=="Eukaryota"){
          names(x)[x=="SAR"] <- "superkingdom2"
        }
      } else {
        names(x) <- taxonomy$V3[match(x ,taxonomy$match)]
      }
      out <- rep(NA, length = length(taxolev.dict1))
      names(out) <- taxolev.dict1
      out[match(names(x), names(out))] <- x
      return(out)
    }))

    tmp1 <- tmp1[match(tmp, gsub("uncultured;$", "", tmp.uniq, perl=T)),]
    colnames(tmp1) <- paste(colnames(tmp1), "silva", sep = "_")

    tmp1 <- cbind.data.frame(silva[, c("cluster.acc", "cluster.id", "similarity", "X..sequences")],
                             tmp1,
                             lineage_silva = tmp)

    # propagate
    heads <- which(tmp1$X..sequences > 1)
    if (length(heads) > 0) {
      clust <- read.fasta(clust.path, as.string = T, forceDNAtolower = F)
      clust <- clust[match(tmp1$cluster.acc, names(clust))]

      tmp2 <- do.call("rbind", lapply(heads, function(x) {
        # existing seq
        h <- as.vector(unname(tmp1$cluster.acc[x]))
        tax <- tmp1[match(h, tmp1$cluster.acc), ]

        gp <- unname(gsub(" >", "", gsub("\\..*", "", unlist(strsplit(unlist(clust[x]), ","))[-1])))
        gp <- gp[grep(h, gp, invert = T)]
        tax2 <- vector("list", length = length(gp))
        for (i in 1:length(tax2)) {
          tax2[[i]] <- tax
        }
        data.frame(cluster.acc = gp, tax[rep(1, length(gp)), -1])
      }))

      colnames(tmp1) <- gsub(" ", "\\.", colnames(tmp1))
      tmp1 <- rbind.data.frame(tmp1, tmp2)
    }

    # merge with motus table
    idx <- match(rownames(motus), as.vector(tmp1$cluster.acc))
    out <- cbind.data.frame(motus, tmp1[idx, -c(1:2, 4)])
    metabarlist$motus <- out

    check_metabarlist(metabarlist)
    return(metabarlist)

  }}

