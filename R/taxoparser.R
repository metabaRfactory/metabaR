#' Parse taxonomic information
#'
#' Parse taxonomic information from full taxonomic paths
#'
#'
#' @param taxopath    a vector containing full taxonomic paths to parse
#' @param sep.level   a character string to separate the taxonomic levels in `taxopath`. NA character not allowed.
#' @param sep.info    a character string to separate taxonomic from taxorank information in `taxopath`. NA character not allowed.
#'
#' @name taxoparser
#'
#' @return a list of vectors containing parsed taxa as values and corresponding taxonomic ranks as value names.
#'
#' @details
#' The taxonomic path should include both taxa names AND their associated taxonomic rank (full names or abbreviations as in qiime or unite outputs). The function will use it together with separators by decreasing level of taxonomic resolution. The taxonomic information should follow a standard structure across samples (e.g. standard taxonomy as in Genbank, SILVA or BOLD by decreasing level of taxonomic resolution: the function does not infer missing taxonomic ranks).
#'
#' @seealso \code{\link{ggtaxplot}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' # Parse taxonomic path
#'
#' ## a ncbi-like type of full taxonomic path
#' taxoparsed <- taxoparser(taxopath = soil_euk$motus$path,
#'                          sep.info = "@",
#'                          sep.level= ":")
#'
#' ## a qiime/unite-like type of full taxonomic path.
#' arthropoda <- subset_metabarlist(soil_euk,
#'                                  table = "motus",
#'                                  indices = grepl("Arthropoda", soil_euk$motus$path))
#'
#' qiimepath <- apply(arthropoda$motus[,grep("[msry]_name", colnames(soil_euk$motus))], 1,
#'                   function(x) {
#'                      paste(sapply(1:length(x), function(y) {
#'                        paste(c("p", "c", "o", "f", "g", "s")[y], x[y], sep="_")
#'                        }), collapse =";")})
#'
#' taxoparsed <- taxoparser(taxopath = qiimepath,
#'                          sep.info = "_",
#'                          sep.level= ";")
#' @author Lucie Zinger
#' @export taxoparser

taxoparser <- function(taxopath, sep.level, sep.info) {
  if (!is.character(taxopath)) {
    stop("`taxo` should be a character string or vector")
  }
  if (is.null(sep.level) | is.null(sep.info)) {
    stop("`sep.level` and `sep.info` are required to parse the taxonomic information")
  }
  if (all(!grepl(sep.info, taxopath))) {
    stop("`sep.info` not found in taxonomic information: no parsing possible")
  }
  if (all(!grepl(sep.level, taxopath))) {
    stop("`sep.level` not found in taxonomic information: no parsing possible")
  }

  # dictionnary of taxo ranks
  # ncbi-like
  taxolev.dict1 <- c(
    "superkingdom", "kingdom", "subkingdom",
    "superphylum", "phylum", "subphylum",
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

  # qiime/unite like
  taxolev.dict2 <- c("k", "p", "c", "o", "f", "g", "s")

  # find were to parse (i.e. taxo level info before or after taxon name)
  parsid <- which(unlist(strsplit(unlist(strsplit(taxopath[1], sep.level)[1])[1], sep.info)) %in%
                    c(taxolev.dict1, taxolev.dict2))

  # set special characters
  spe <- c("$", "*", "+", ".", "?", "[", "]", "^", "{", "}", "|", "(", ")", "\\")
  if (sep.info %in% spe) sep.info <- paste("\\", sep.info, sep = "")
  if (sep.level %in% spe) sep.level <- paste("\\", sep.level, sep = "")

  #parse
  tmp <- sapply(taxopath, strsplit, sep.level)
  parse <- lapply(tmp, function(x) {
    tmp2 <- unlist(strsplit(x, sep.info))
    if (parsid == 2) {
      idx1 <- seq(parsid, sum(grepl(sep.info, x))*2, 2)-1
      idx2 <- seq(parsid, sum(grepl(sep.info, x))*2, 2)
    } else {
      idx1 <- seq(parsid, sum(grepl(sep.info, x))*2, 2)+1
      idx2 <- seq(parsid, sum(grepl(sep.info, x))*2, 2)
    }
    out <- tmp2[idx1]
    names(out) <- tmp2[idx2]
    out
  })
  parse
}
