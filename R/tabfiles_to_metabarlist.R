#' Import 4 tables to create a metabarlist object
#'
#' Importing and formatting four files (corresponding to MOTU abundances, MOTU characteristics, PCR characteristics and sample characteristics respectively) to create a \code{\link{metabarlist}} object.
#'
#'
#' @param file_reads    path for the MOTU abundance table. Rows of the table should correspond to PCRs,
#'                      and columns should correspond to MOTUs. The first column of this table should
#'                      correspond to the names of the pcrs; the first line to the names of MOTUs.
#' @param file_motus    path for the MOTU characteristics table (e.g. taxonomy, sequence, etc.).
#'                      Rows of the table should correspond to MOTUs, and the columns to their
#'                      characteristics. The first column of this table should contain MOTU names.
#'                      Mandatory fields: `sequence`, i.e. the most abundant sequence of the MOTU.
#' @param file_pcrs     path for the pcrs characteristics table (e.g. tags, primers, plate wells, etc.).
#'                      Mandatory fields: (i) `sample_id`, i.e. the name of each sample. (ii) `type`,
#'                      i.e. the type of pcr; can be `sample` or `control`. (iii) `control_type`,
#'                      i.e. the type of control if applicable. Should be: `NA` for samples,
#'                      `extraction` for extraction negative controls, `pcr` for pcr negative controls,
#'                      `sequencing` for sequencing negative controls (e.g. unused tag combinations),
#'                      and `positive` for positive controls. The first column of this table should
#'                      correspond to the names of the pcrs.
#' @param file_samples  path for the sample characteristics table. The first column of this table
#'                      should contain the sample names.
#' @param files_sep     separator used to read the different tables. Should be the same in all tables.
#'                      Default is tabulation.
#'
#'
#' @name tabfiles_to_metabarlist
#'
#' @return a \code{\link{metabarlist}} object
#'
#' @details
#'
#' This function imports tabular data into \R to create a \code{\link{metabarlist}} object. The four files required are imported in \R, included into a list of class \code{\link{metabarlist}} with the \code{\link{metabarlist_generator}} function, and congruencies between all tables are tested with the \code{\link{check_metabarlist}} function.
#'
#' @examples
#'
#' soil_euk <- tabfiles_to_metabarlist(
#'   file_reads = system.file("extdata", "litiere_euk_reads.txt", package = "metabaR"),
#'   file_motus = system.file("extdata", "litiere_euk_motus.txt", package = "metabaR"),
#'   file_pcrs = system.file("extdata", "litiere_euk_pcrs.txt", package = "metabaR"),
#'   file_samples = system.file("extdata", "litiere_euk_samples.txt", package = "metabaR")
#' )
#'
#' @seealso \code{\link{check_metabarlist}}, \code{\link{metabarlist_generator}},
#'           \code{\link{obifiles_to_metabarlist}}, \code{\link{biomfiles_to_metabarlist}}
#'
#' @author Lucie Zinger & Clément Lionnet & Frédéric Boyer
#' @importFrom utils read.table read.csv2
#' @export tabfiles_to_metabarlist
#'

tabfiles_to_metabarlist <- function(file_reads, file_motus, file_pcrs, file_samples, files_sep = "\t") {
  if (!file.exists(file_reads)) {
    stop(paste("cannot open file_reads", file_reads,": No such file or directory"))
  }
  if (!file.exists(file_motus)) {
    stop(paste("cannot open file_reads", file_motus,": No such file or directory"))
  }
  if (!file.exists(file_pcrs)) {
    stop(paste("cannot open file_reads", file_pcrs,": No such file or directory"))
  }
  if (!file.exists(file_samples)) {
    stop(paste("cannot open file_reads", file_samples,": No such file or directory"))
  }


  reads <- as.matrix(read.csv2(file_reads,
    row.names = 1, h = T, sep = files_sep,
    check.names = F, stringsAsFactors = F
  ))
  motus <- read.table(file_motus,
    row.names = 1, h = T, sep = files_sep,
    check.names = F, stringsAsFactors = F
  )
  pcrs <- read.table(file_pcrs,
    row.names = 1, h = T, sep = files_sep,
    check.names = F, stringsAsFactors = F
  )
  samples <- read.table(file_samples,
    row.names = 1, h = T, sep = files_sep,
    check.names = F, stringsAsFactors = F
  )

  if (!all(rownames(reads) %in% rownames(pcrs))) {
    stop("cannot continue, rownames in reads are not part of rownames of pcrs")
  }

  reads <- reads[match(rownames(pcrs), rownames(reads)), ]
  reads[is.na(reads)] <- 0
  rownames(reads) <- rownames(pcrs)

  out <- metabarlist_generator(reads, motus, pcrs, samples)

  return(out)
}
