#' Generate a fasta file on a subset of sequences
#'
#' Generate a fasta file on a subset of OTUs from a \code{\link{TODEFINE}} object.
#'
#'
#' @param file    path/name of the output
#' @param id     a character vector containing identifiers for the OTUs of interest
#' @param seq   a character vector containing sequences for the OTUs of interest
#'
#' @name fasta_generator
#'
#' @return a fasta file
#'
#' @details
#' Creates a fasta file from a set of sequences of interest. Obviously, \code{id} and \code{seq} vectors should be in agreement. The header can easily be amended with further information with e.g. \code{paste}.
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # export in fasta the 10 most abundant MOTUs
#' idx <- order(soil_euk$motus$count, decreasing = T)[1:10]
#' fasta_generator(
#'   soil_euk, rownames(soil_euk$motus)[idx],
#'   "Dominants.fasta"
#' )
#'
#' # export in fasta the 10 most abundant MOTUs and their abundance
#' idx <- order(soil_euk$motus$count, decreasing = T)[1:10]
#' annotation <- data.frame(
#'   abundance = soil_euk$motus$count[idx],
#'   row.names = rownames(soil_euk$motus)[idx]
#' )
#' fasta_generator(
#'   soil_euk, rownames(soil_euk$motus)[idx],
#'   "Dominants.fasta",
#'   annotation
#' )
#' @author Lucie Zinger
#' @export fasta_generator

fasta_generator <- function(metabarlist, id = rownames(metabarlist$motus), output_file = NULL, annotation = NULL, annot_sep = " ") {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (is.null(output_file)) {
      stop("The output_file is required!")
    }

    if (!all(id %in% rownames(metabarlist$motus))) {
      stop("All identifiers must correspond exactly to rows of metabarlist$motus")
    }

    if (!is.null(annotation)) {
      if (!all(id %in% rownames(annotation))) {
        stop("All identifiers must correspond exactly to rows of annotation")
      }

      header_annot <- apply(annotation, 1, function(x) {
        paste(names(x), x, sep = "=", collapse = annot_sep)
      })
      header <- sapply(id, function(x) {
        bh <- paste(">", x, sep = "")
        paste(bh, header_annot[x], sep = annot_sep)
      })
    }
    else {
      header <- sapply(id, function(x) {
        paste(">", x, sep = "")
      })
    }

    if (file.exists(output_file)) {
      file.remove(output_file)
    }
    file.create(output_file)

    for (i in id) {
      cat(paste(header[i], metabarlist$motus[i, "sequence"], sep = "\n"),
        sep = "\n", file = output_file, append = T
      )
    }
  }
}


# fasta_generator <- function(id, seq, output_file = NULL) {
#   if (is.null(output_file)) {
#     stop("The output_file is required!")
#   }
#   if (length(id) != length(seq)) {
#     stop("id and seq should be of the same length")
#   }
#
#   if (file.exists(output_file)) {
#     file.remove(output_file)
#   }
#   file.create(output_file)
#
#   for (i in 1:length(id)) {
#     cat(paste(">", id[i], "\n", seq[i], "\n", sep = ""),
#       file = output_file, append = T
#     )
#   }
# }
