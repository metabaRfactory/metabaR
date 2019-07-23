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
#'   "Dominants.fasta", rownames(soil_euk$motus)[idx],
#'   soil_euk$motus$sequence[idx]
#' )
#'
#' # export in fasta the 10 most abundant MOTUs and their abundance
#' idx <- order(soil_euk$motus$count, decreasing = T)[1:10]
#' headers <- paste(rownames(soil_euk$motus)[idx],
#'   "abundance=", soil_euk$motus$count[idx],
#'   sep = ""
#' )
#' fasta_generator(
#'   "Dominants.fasta", headers,
#'   soil_euk$motus$sequence[idx]
#' )
#' @author Lucie Zinger
#' @export fasta_generator

fasta_generator <- function(file, id, seq) {
  if (length(id) != length(seq)) {
    stop("id and seq should be of the same length")
  }

  if (file.exists(file)) {
    file.remove(file)
  }
  file.create(file)

  for (i in 1:length(id)) {
    cat(paste(">", id[i], "\n", seq[i], "\n", sep = ""),
      file = file, append = T
    )
  }
}
