#' Import data to create a metabarlist object
#'
#' Importing and formatting files to create a \code{\link{metabarlist}} object.
#'
#'
#' @param file_reads    the path for the MOTU abundance table. Rows of the table should correspond to pcrs, and the columns to MOTUs. The first column of this table should correspond to the names of the pcrs.
#' @param file_motus    the path for the MOTU characteristics table (e.g. taxonomy, sequence, etc.). Rows of the table should correspond to MOTUs, and the columns to their characteristics. The first column of this table should contain MOTUs names. Mandatory fields: `sequence`, i.e. the most abundant sequence of the MOTU.
#' @param file_pcrs     the path for the pcrs characteristics table (e.g. tags, primers, plate wells, etc.). Mandatory fields: (i) `sample_id`, i.e. the name of each sample. (ii) `type`, i.e. the type of pcr; can be `sample` or `control`. (iii) `control_type`, i.e. the type of control if applicable. Should be: `NA` for samples, `extraction` for extraction negative controls, `pcr` for pcr negative controls, `sequencing` for sequencing negative controls (e.g. unused tag combinations), and `positive` for positive controls.
#' @param file_samples  the path for the sample characteristics table. The first column of this table shuld contain the sample names.
#' @param files_sep   separator used to read the different tables. Should be the same in all tables. Default is tabulation.
#'
#'
#' @name metabarlist_generator
#'
#' @return a \code{\link{metabarlist}} object
#'
#' @details
#'
#' This function aims at importing data into \R to create a \code{\link{metabarlist}} object. The four files required are imported in \R, included into a list of class \code{\link{metabarlist}}, and congruencies between all tables are tested with the \code{\link{check_metabarlist}} function.
#'
#' @examples
#'
#' @seealso \code{\link{check_metabarlist}}
#'
#' @author Lucie Zinger & Clément Lionnet
#' @export metabarlist_generator
#'

metabarlist_generator = function(file_reads, file_motus, file_pcrs, file_samples, files_sep = "\t") {
  if(!file.exists(file_reads))
    stop("file reads does not exist")
  if(!file.exists(file_motus))
    stop("file motus does not exist")
  if(!file.exists(file_pcrs))
    stop("file pcrs does not exist")
  if(!file.exists(file_samples))
    stop("file samples does not exist")


  reads = as.matrix(read.csv2(file_reads, row.names=1, h=T, sep=files_sep,
                              check.names = F, stringsAsFactors = F))
  motus = read.table(file_motus, row.names=1, h=T, sep=files_sep,
                     check.names = F, stringsAsFactors = F)
  pcrs = read.table(file_pcrs, row.names=1, h=T, sep=files_sep,
                    check.names = F, stringsAsFactors = F)
  samples = read.table(file_samples, row.names=1, h=T, sep=files_sep,
                       check.names = F, stringsAsFactors = F)

  if(!all(rownames(reads) %in% rownames(pcrs)))
    stop("cannot continue, rownames in reads are not part of rownames of pcrs")

  reads = reads[match(rownames(pcrs), rownames(reads)),]
  reads[is.na(reads)] = 0
  rownames(reads) = rownames(pcrs)

  out = list(reads = reads,
             motus = motus,
             pcrs = pcrs,
             samples = samples)

  attr(out, 'class') = "metabarlist"

  check_metabarlist(out)

  return(out)
}
