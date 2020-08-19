#' Metabarlist summary
#'
#' summary_metabarlist is a generic function used to produce summary statistics of a \code{metabarlist} object.
#'
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param method        type of summary to provide. Should match with `dataset`, `motus` or `pcrs`.
#'                      Default is `dataset`
#' @param groups        a grouping vector or factor of same number of rows than `motus` or `pcrs`
#'                      for which the summary should be done. Default is NULL
#'
#' @return The \code{summary_metabarlist} function returns basic summary statistics (nb of elements, reads and MOTUs in total or on average per pcrs) of a \code{metabarlist} object. The format of the value returned depends on the method used.
#'
#' @details
#'
#' \code{summary_metabarlist} returns basic summary statistics of a \code{metabarlist} object. The summary returned depends on the `method` used :
#' \itemize{
#' \item{`dataset`}{returns a list of two data.frames: dataset_dimension contains the dimensions of the full metabarlist object, dataset_statistics contains the number of reads, motus in pcrs and samples, as well as average and sd values of reads and motus per sample}
#' \item{`motus` or `pcrs`}{returns a data.frame similar to the dataset_statistics described above according to a grouping factor/vector for MOTUs or pcrs}
#' }
#'
#' @seealso \code{\link{summary}}
#' @examples
#'
#' data(soil_euk)
#'
#' ## Dataset summary
#' summary_metabarlist(soil_euk, method = "dataset")
#'
#' ## Data summary per control type (NA = samples)
#' summary_metabarlist(soil_euk, method = "pcrs",
#'     groups = soil_euk$pcrs$control_type)
#'
#' ## Data summary per phyla
#' summary_metabarlist(soil_euk, method = "motus",
#'     groups = soil_euk$motus$phylum_name)
#'
#' ## Data summary per Habitat (i.e. to get from soil_euk$samples).
#' # Here, NA values correspond to technical controls
#' summary_metabarlist(soil_euk, method = "pcrs",
#'     groups = soil_euk$samples$Habitat[match(soil_euk$pcrs$sample_id,
#'                                             rownames(soil_euk$samples))])
#'
#' @author Lucie Zinger
#' @importFrom stats sd xtabs aggregate
#' @export summary_metabarlist

summary_metabarlist <- function(metabarlist, method = "dataset", groups = NULL) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    extract_type_methods <- c("dataset", "motus", "pcrs")
    method <- match.arg(method, extract_type_methods)

    if (method == "dataset") {
      dataset_dimension <- data.frame(t(sapply(metabarlist, dim)))
      colnames(dataset_dimension) <- c("n_row", "n_col")

      idx <- which(metabarlist$pcrs$type == "sample")
      size <- rowSums(metabarlist$reads)
      rich <- rowSums(metabarlist$reads > 0)

      dataset_statistics <- data.frame(
        nb_reads = c(sum(metabarlist$reads), sum(metabarlist$reads[idx, ])),
        nb_motus = c(
          sum(colSums(metabarlist$reads) > 0),
          sum(colSums(metabarlist$reads[idx, ]) > 0)
        ),
        avg_reads = c(mean(size), mean(size[idx])),
        sd_reads = c(sd(size), sd(size[idx])),
        avg_motus = c(mean(rich), mean(rich[idx])),
        sd_motus = c(sd(rich), sd(rich[idx])),
        row.names = c("pcrs", "samples")
      )

      return(list(
        dataset_dimension = dataset_dimension,
        dataset_statistics = dataset_statistics
      ))

    } else {

      tab <- match.arg(method, extract_type_methods)

      if (is.null(groups)) {
        stop(paste("vector or factor `groups` should be provided for method '", tab, "'", sep=""))
      }

      if (length(groups) != nrow(metabarlist[[tab]])) {
        stop(paste("vector or factor `groups` should be of same length than the object
                   called with method '",
                   tab, "'", sep=""))
      }

      idx <- ifelse(is.na(as.vector(groups)), "NA", as.vector(groups))

      if (tab == "pcrs" | tab == "motus") {
        size <- if (tab == "pcrs") {
          rowSums(metabarlist$reads)
        } else {
          colSums(metabarlist$reads)
        }
        rich <- if (tab == "pcrs") {
          rowSums(metabarlist$reads > 0)
        } else {
          colSums(metabarlist$reads > 0)
        }

        dataset_statistics <- data.frame(
          nb_elements = if (tab == "pcrs") {
            as.vector(table(idx))
          } else {
            NA
          },
          nb_reads = as.vector(xtabs(size ~ idx)),
          nb_motus = if (tab == "pcrs") {
            as.vector(xtabs(rich ~ idx))
          } else {
            as.vector(table(idx))
          },
          avg_reads = aggregate(size, list(idx), mean)$x,
          sd_reads = aggregate(size, list(idx), sd)$x,
          avg_motus = aggregate(rich, list(idx), mean)$x,
          sd_motus = aggregate(rich, list(idx), sd)$x,
          row.names = levels(factor(idx))
        )
        if (tab == "pcrs") {
          return(dataset_statistics)
        } else {
          return(dataset_statistics[, -1])
        }
      }
    }
  }
}
