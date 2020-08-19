#' Generating a fasta file from a \code{metabarlist} object
#'
#' Generates a fasta file from a \code{metabarlist} object, typically for a subset of MOTUs of interest.
#'
#' @param metabarlist   a \code{metabarlist} object.
#' @param id            a character vector containing identifiers for the MOTUs of interest.
#'                      Default is the row names of table `motus`
#' @param output_file   the path/name of the output file.
#' @param annotation    a data frame containing the additional information to be added
#'                      in the sequence header. Each row name is a MOTU identifier and each
#'                      column header is the information key
#'                      (correspond to "key=" in the fasta header).
#'                      Default: NULL
#' @param annot_sep     the annotation separator. Default is ";"
#'
#' @name fasta_generator
#'
#' @return a fasta file
#'
#' @details
#'
#' Creates a fasta file from a set of sequences of interest. A single or multiple pieces of sequence information can be added to the sequence header via the `annotation` table. Using the default parameters, the resulting fasta file will have the following format:
#'
#'>Seq_id1 abundance=83079; GC=47
#'tcaatctcgtgtgactaaacgccacttgtccctctaagaagttacgccgacagaatgcgatcggcgaactatttagcaggctagagtctcgttcgttat
#'>Seq_id2 abundance=120520; GC=55
#'ctcaaacttccttggcctggaaggccatagtccctctaagaagctggccgcggagggtcacctccgcatagctagttagcaggctgaggtctcgttcgttaa
#'
#' @seealso \code{\link{metabarlist_generator}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' ## Export in fasta format the 10 most abundant MOTUs
#' idx <- order(soil_euk$motus$count, decreasing = TRUE)[1:10]
#' fasta_generator(soil_euk, rownames(soil_euk$motus)[idx],
#'   "Dominants.fasta"
#' )
#'
#' ## Export in fasta format the 10 most abundant MOTUs, their abundance and the GC content
#' # of the corresponding sequence
#' annotation <- data.frame(
#'   abundance = soil_euk$motus$count[idx],
#'   GC_content = soil_euk$motus$GC_content[idx],
#'   row.names = rownames(soil_euk$motus)[idx]
#' )
#'
#' dir <- tempdir()
#'
#' fasta_generator(
#'   soil_euk, rownames(soil_euk$motus)[idx],
#'    file.path(dir, "Dominants.fasta"),
#'   annotation
#' )
#' @author Lucie Zinger, Clément Lionnet
#' @export fasta_generator

fasta_generator <- function(metabarlist, id = rownames(metabarlist$motus), output_file = NULL, annotation = NULL, annot_sep = ";") {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    if (is.null(output_file)) {
      stop("output_file is required!")
    }

    if (!all(id %in% rownames(metabarlist$motus))) {
      stop(
        "values of id are not in the row names of table `motus`"
      )
    }

    if (!is.null(annotation)) {
      if (!all(id %in% rownames(annotation))) {
        stop(
          "elements of id (identifiers of the MOTUs to select) must be in the rownames
              of the 'annotation' data.frame"
        )
      }

      header_annot <- apply(annotation, 1, function(x) {
        paste(paste(" ", names(x), sep = ""),
              x,
              sep = "=",
              collapse = annot_sep)
      })

      header <- sapply(id, function(x) {
        bh <- paste(">", x, sep = "")
        paste(bh, header_annot[x], sep = "") #output as obitools.
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

