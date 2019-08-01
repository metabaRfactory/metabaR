#' Create a metabarlist object
#'
#' Formatting R tables to create a \code{\link{metabarlist}} object.
#'
#'
#' @param reads     MOTU abundance table. Rows and rownames of the table should correspond to PCRs and their names respectively. Columns and colnames should correspond to MOTUs and their names. Rownames in this table should correspond to PCR names respectively.
#' @param motus     MOTU characteristics table (e.g. taxonomy, sequence, etc.). Rows and rownames of the table should correspond to MOTUs and their names respectively, and the columns to their characteristics. Mandatory fields: `sequence`, i.e. the sequence representative of the MOTU.
#' @param pcrs      PCR characteristics table (e.g. tags, primers, plate wells, etc.). Rows and rownames of the table should correspond to PCRs and their names respectively, and the columns to their characteristics. Mandatory fields: (i) `sample_id`, i.e. the name of each biological sample. (ii) `type`, i.e. the type of PCR; can be `sample` or `control`. (iii) `control_type`, i.e. the type of control if applicable. Should be either: `NA` for samples, `extraction` for extraction negative controls, `pcr` for PCR negative controls, `sequencing` for sequencing negative controls (e.g. unused tag combinations), or `positive` for positive controls.
#' @param samples   Samples characteristics table. Rows and rownames of the table should correspond to biological samples and their names respectively, and the columns to their environnemental characteristics.
#'
#'
#' @name metabarlist_generator
#'
#' @return a \code{\link{metabarlist}} object
#'
#' @details
#'
#' This function formats R tables to create a \code{\link{metabarlist}} object. The four objects required are incorporated into a list of class \code{\link{metabarlist}}. Congruencies between all tables are tested internally with the \code{\link{check_metabarlist}} function.
#'
#' @examples
#'
#' # Create fake data
#' ## MOTUs abundance table
#'
#' reads <- matrix(sample(1:1000, 200, replace = TRUE), nrow = 10, ncol = 20)
#' rownames(reads) <- paste(c(rep(c("A", "B", "C"), each = 2), "ex", "pcr", "seq", "pos"),
#'   c(rep(c("r1", "r2"), 3), rep(0, 4)),
#'   sep = "_"
#' )
#' colnames(reads) <- sprintf("MOTU_%03d", 1:ncol(reads))
#'
#' ## MOTUs characteristics table
#' motus <- data.frame(
#'   fake_taxon = sample(letters, ncol(reads), replace = T),
#'   sequence = sapply(
#'     1:ncol(reads),
#'     function(x) paste(sample(c("a", "t", "c", "g"), 20, TRUE), collapse = "")
#'   ),
#'   row.names = colnames(reads), stringsAsFactors = FALSE
#' )
#'
#' ## PCR characteristics table
#' pcrs <- data.frame(
#'   sample_id = sapply(strsplit(rownames(reads), "_"), "[[", 1),
#'   type = c(rep("sample", 6), rep("control", 4)),
#'   control_type = c(rep(NA, 6), "extraction", "pcr", "sequencing", "positive"),
#'   fake_pcrplate = sample(c(1, 2), nrow(reads), replace = T),
#'   row.names = rownames(reads)
#' )
#'
#' ## Sample characteristics table
#' samples <- data.frame(
#'   fake_latitude = 1:3, fake_longitude = 1:3,
#'   row.names = unique(pcrs$sample_id[is.na(pcrs$control_type)])
#' )
#'
#' # Generate the metabarlist object
#'
#' test <- metabarlist_generator(reads = reads, motus = motus, pcrs = pcrs, samples = samples)
#' ## Warnings are returned because the PCR design (i.e. tags, primers, plate coordinates) is not defined in this example.
#' @seealso \code{\link{check_metabarlist}}
#'
#' @author Lucie Zinger & Clément Lionnet & Frédéric boyer
#' @export metabarlist_generator
#'

metabarlist_generator <- function(reads, motus, pcrs, samples) {
  reads <- reads[match(rownames(pcrs), rownames(reads)), ]
  reads[is.na(reads)] <- 0
  rownames(reads) <- rownames(pcrs)

  out <- list(
    reads = reads,
    motus = motus,
    pcrs = pcrs,
    samples = samples
  )

  attr(out, "class") <- "metabarlist"

  check_metabarlist(out)

  return(out)
}
