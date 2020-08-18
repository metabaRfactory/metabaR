#' Import BIOM and associated files to create a metabarlist object
#'
#' Imports and formats BIOM and associated files to create a \code{\link{metabarlist}} object.
#'
#'
#' @param file_biom      path for the \code{BIOM} file. This is either a JSON formatted file
#'                       (biom file format version 1) or a HDF5 formatted file (\code{BIOM} file format
#'                       version 2 and 2.1), as described in \link{http://biom-format.org/}.
#'                       This file should include at least MOTUs abundance data.
#'                       It may also store MOTUs and/or PCRs attributes data.
#'                       Mandatory fields for MOTUs and PCRs attributes data are described below.
#' @param file_samples   path for the sample characteristics table.
#'                       The first column of this table should contain the sample names.
#' @param file_pcrs      path for the PCRs characteristics table (e.g. tags, primers, plate wells, etc.),
#'                       if the \code{BIOM} file is missing these data. Mandatory fields:
#'                       (i) `sample_id`, i.e. the name of each sample. (ii) `type`,
#'                       i.e. the type of PCR; can be `sample` or `control`. (iii) `control_type`,
#'                       i.e. the type of control if applicable. Should be: 'NA' for samples,
#'                       `extraction` for extraction negative controls, `pcr` for pcr negative controls,
#'                       `sequencing` for sequencing negative controls (e.g. unused tag combinations),
#'                       and `positive` for positive controls. The first column of this table should
#'                       correspond to the names of the PCRs.
#' @param file_motus     path for the MOTUs characteristics table (e.g. taxonomy, sequence, etc.),
#'                       if the \code{BIOM} file is missing these data. Rows of the table should
#'                       correspond to MOTUs, and the columns to their characteristics.
#'                       Mandatory fields: 'sequence', i.e. the most abundant sequence of the MOTU.
#'                       The first column of this table should contain MOTU names.
#' @param ... other arguments to be pasted from \code{read.table}.
#'
#' @name biomfiles_to_metabarlist
#'
#' @return a \code{\link{metabarlist}} object
#'
#' @details
#'
#' @seealso \code{\link{check_metabarlist}}, \code{\link{metabarlist_generator}},
#'           \code{\link{obifiles_to_metabarlist}}, \code{\link{tabfiles_to_metabarlist}}
#'
#' @references
#' \link{http://biom-format.org/}
#'
#' @examples
#' add an example here, such as:
#'
#'soil_euk <- biomfiles_to_metabarlist(
#'  file_biom = system.file("extdata", "litiere_euk_reads_hdf5.biom", package = "metabaR"),
#'  file_motus = system.file("extdata", "litiere_euk_motus.txt", package = "metabaR"),
#'  file_pcrs = system.file("extdata", "litiere_euk_pcrs.txt", package = "metabaR"),
#'  file_samples = system.file("extdata", "litiere_euk_samples.txt", package = "metabaR"),
#'  sep = "\t")
#'
#' @author Anne-Sophie Benoiston, Lucie Zinger
#'
#' @import biomformat
#'
#' @export biomfiles_to_metabarlist

biomfiles_to_metabarlist <- function(file_biom, file_samples, file_pcrs = NULL, file_motus = NULL, ...) {
  if (!file.exists(file_biom)) {
    stop(paste("cannot open file_biom", file_biom,": No such file or directory"))
  }
  if (!file.exists(file_samples)) {
    stop(paste("cannot open file_samples", file_samples,": No such file or directory"))
    }

  biom <- suppressWarnings(read_biom(biom_file = file_biom))

  # reads
  reads <- t(as.matrix(biom_data(biom)))

  # motus
  if (is.null(observation_metadata(biom))) {
    if(missing(file_motus)) {
      stop("No metadata on MOTUs: a motus file is required")
    }
    else {
      if (!file.exists(file_motus)) {
        stop(paste("cannot open file_motus", file_motus,": No such file or directory"))
      }
      else {
        motus <- read.table(file_motus,
                            row.names = 1, h = T,
                            check.names = F, stringsAsFactors = F, ...)
      }
    }
  }
  else {
    motus <- observation_metadata(biom)
  }

  # pcrs
  if (is.null(sample_metadata(biom))) {
    if(missing(file_pcrs)) {
      stop("No metadata on PCRs: a pcrs file is required")
    }
    else {
      if (!file.exists(file_pcrs)) {
        stop(paste("cannot open file_pcrs", file_pcrs,": No such file or directory"))
      }
      else {
        pcrs <- read.table(file_pcrs,
                           row.names = 1, h = T,
                           check.names = F, stringsAsFactors = F, ...)
      }
    }
  }
  else {
    pcrs <- sample_metadata(biom)
    pcrs$control_type[pcrs$control_type == "NA"] <- NA
    pcrs$plate_no = as.numeric(pcrs$plate_no)
    pcrs$plate_col = as.numeric(pcrs$plate_col)
    pcrs$plate_row = as.factor(pcrs$plate_row)
    pcrs$tag_fwd = as.factor(pcrs$tag_fwd)
    pcrs$tag_rev = as.factor(pcrs$tag_rev)
  }

  # samples
  samples <- read.table(file_samples,
                        row.names = 1, h = T,
                        check.names = F, stringsAsFactors = F, ...
  )

  # check pcrs in reads present in pcrs table
  if (!all(rownames(reads) %in% rownames(pcrs))) {
    stop("cannot continue, rownames in reads are not part of rownames of pcrs")
  }

  out <- metabarlist_generator(reads, motus, pcrs, samples)

  return(out)
}
