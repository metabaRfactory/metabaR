#' Import OBITools and associated files to create a metabarlist object
#'
#' Importing and formatting OBITools and associated files (i.e. an \code{obitab} output file,
#' an \code{ngsfilter} file and table of sample characteristics respectively) to create
#' a \code{metabarlist} object (see Boyer et al. 2016).
#'
#'
#' @param file_obitab       path for the \code{obitab} output file. Rows of the table correspond to MOTUs,
#'                          and columns correspond to MOTU characteristics and counts across PCRs.
#'                          Mandatory fields: `sequence`; should be included when using \code{obitab}.
#' @param file_ngsfilter    path for the \code{ngsfilter} file. Rows of the table correspond to PCRs,
#'                          and the columns to their characteristics. Mandatory fields in the additional
#'                          information: (i) `sample_id`, i.e. the name of each sample. (ii) `type`, i.e.
#'                          the type of pcr; can be `sample` or `control`. (iii) `control_type`, i.e. the
#'                          type of control if applicable. Should be: `NA` for samples, `extraction` for
#'                          extraction negative controls, `pcr` for pcr negative controls, `sequencing`
#'                          for sequencing negative controls (e.g. unused tag combinations), and `positive`
#'                          for positive controls. The first column of this table should correspond to the
#'                          names of the pcrs.
#' @param file_samples      path for the sample characteristics table. The first column of this
#'                          table should contain the sample names.
#' @param ...               other arguments to be pasted from \code{read.table}
#'
#'
#' @name obifiles_to_metabarlist
#'
#' @return a \code{metabarlist} object
#'
#' @details
#'
#' This function imports OBITools outputs and related files into \R to create a \code{metabarlist} object. The three files required are imported in \R, included into a list of class \code{metabarlist} with the \code{\link{metabarlist_generator}} function, before congruencies between all tables are tested with the \code{\link{check_metabarlist}} function.
#'
#'
#' @seealso \code{\link{tabfiles_to_metabarlist}}, \code{\link{biomfiles_to_metabarlist}}
#'
#' @references Boyer, F., Mercier, C., Bonin, A., Le Bras, Y., Taberlet, P., & Coissac, E. (2016). obitools: a unix-inspired software package for DNA metabarcoding. Molecular Ecology Resources, 16(1), 176-182.
#'
#' @examples
#'
#' \donttest{
#'
#' dir <- tempdir()
#' url = "https://raw.githubusercontent.com/metabaRfactory/metabaR_external_data/master/"
#' 
#' obitab_file = "litiere_euk_cl97_agg_filt_tax.tab"
#' obitab_url = paste(url, obitab_file, sep="")
#' obitab_path <- file.path(dir, obitab_file)
#' download.file(obitab_url, obitab_path)
#'
#' ngsfilter_file = "ngsfilter_GWM-768.new_2.txt"
#' ngsfilter_url = paste(url, ngsfilter_file, sep="")
#' ngsfilter_path <- file.path(dir, "ngsfilter_GWM-768.new_2.txt")
#' download.file(ngsfilter_url, ngsfilter_path)
#'
#' samples_file = "Litiere_sample_list_2.txt"
#' samples_url = paste(url, samples_file, sep="")
#' samples_path <- file.path(dir, samples_file)
#' download.file(samples_url, samples_path)
#'
#' soil_euk <- obifiles_to_metabarlist(
#'   file_obitab = obitab_path,
#'   file_ngsfilter = ngsfilter_path,
#'   file_samples = samples_path,
#'   sep = "\t"
#' )
#'
#' }
#'
#' @seealso \code{\link{check_metabarlist}}, \code{\link{metabarlist_generator}}
#'
#' @author Lucie Zinger
#' @importFrom utils read.csv2
#' @export obifiles_to_metabarlist
#'

obifiles_to_metabarlist <- function(file_obitab, file_ngsfilter, file_samples, ...) {
  if (!file.exists(file_obitab)) {
    stop(paste("cannot open file_obitab", file_obitab, ": No such file or directory"))
  }
  if (!file.exists(file_ngsfilter)) {
    stop(paste("cannot open file_ngsfilter", file_ngsfilter, ": No such file or directory"))
  }
  if (!file.exists(file_samples)) {
    stop(paste("cannot open file_samples", file_samples, ": No such file or directory"))
  }

  obi <- read.csv2(file_obitab, header = T, check.names = F, stringsAsFactors = F, ...)

  # reads
  reads <- t(obi[, grep("sample\\:", colnames(obi))])
  rownames(reads) <- gsub("sample\\:", "", rownames(reads))
  colnames(reads) <- obi$id

  # motus
  motus <- obi[, grep("sample\\:", colnames(obi), invert = T)]
  rownames(motus) <- motus$id
  motus <- motus[, -match("id", colnames(motus))]

  # pcrs
  pcrs <- read_ngsfilter(file = file_ngsfilter, additional.sep = "=", ...)
  rownames(pcrs) <- pcrs$pcr_id
  pcrs <- pcrs[, -match("pcr_id", colnames(pcrs))]

  # samples
  samples <- read.csv2(file_samples, row.names = 1, header = T, check.names = F, stringsAsFactors = F, ...)

  # check pcrs in reads present in pcrs table
  if (!all(rownames(reads) %in% rownames(pcrs))) {
    stop("cannot continue, one or several rownames in the metabarlist table `reads`
         are not in the rownames of the metabarlist table `pcrs`")
  }

  # Add null lines for missing pcrs in reads
  reads <- reads[match(rownames(pcrs), rownames(reads)), ]
  reads[is.na(reads)] <- 0
  rownames(reads) <- rownames(pcrs)

  out <- metabarlist_generator(reads, motus, pcrs, samples)

  return(out)
}
