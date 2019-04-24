#' Check that a list of tables is a well formed metabarlist
#'
#' Test that a list of tables containing information on MOTUs abundances, MOTUs characteristics, pcrs characteristics, and sample characteristics form a congruent \code{\link{metabarlist}}
#'
#' @param metabarlist a \code{metabarlist} object
#'
#' @name check_metabarlist
#'
#' @return TRUE or throws a stop
#'
#' @details
#'
#' Check for the properties awaited for a well formed \code{metabarlist} object:
#'
#' \itemize{
#' \item {is a list with four attributes named `reads`, `motus`, `pcrs` and `samples`}
#' \item {`reads` is a numeric matrix}
#' \item {`motus`, `pcrs` and `samples` are of type `data.frame`}
#' \item {`reads` have the exact same row names and columns as rownames of `motus` and `pcrs` (same order also)}
#' \item {`pcrs` have mandatory columns, i.e.  `sample_id`, `type` and `control_type`}
#' \item {values in `type` are properly defined, i.e. `sample` or `control`}
#' \item {values in and `control_type` are properly defined, i.e. `sequencing`, `pcr`, `extraction`, `positive`}
#'}
#'
#' The function issues a stop if these informations are not well recorded
#'
#' The function issue warnings if some tables lacks of non mandatory columns for the \code{metabaRffe} package to run, but mandatory for particular functions (e.g., \code{ggpcrplate})
#' \itemize{
#' \item {the column `sequence` for the `motus` table}
#' \item {the columns `tag_fwd`, `tag_rev`, `primer_fwd`, `primer_rev`, `plate_no`, `plate_col`, and `plate_row` for the `pcrs` table}
#' }
#'
#'
#' Moreover, the function issues a warning if any sample or MOTU is associated to a count of 0
#'
#' @examples
#'
#' data(soil_euk)
#'
#' check_metabarlist(soil_euk)
#'
#' @author Clément Lionnet & Frédéric Boyer & Lucie Zinger
#' @export check_metabarlist


check_metabarlist <- function(metabarlist) {

  if ( ! "metabarlist" %in% class(metabarlist)) {
    stop("Not a metabarlist class")
  }

  if ( ! (is.list(metabarlist))) {
    stop("metabarlist is not a list")
  }
  slots <- c("reads", "motus", "pcrs", "samples")
  if ( ! (all(slots %in% names(metabarlist)))) {
    stop(paste0("metabarlist has no tag(s): ", paste(slots[slots %in% names(metabarlist)], collapse=", ")))
  }


  if ( ! (is.matrix(metabarlist$reads) &&
          is.numeric(metabarlist$reads) &&
          ! any(is.na(metabarlist$reads) && all(metabarlist$reads>=0)))) {
    stop("metabarlist$reads must be a positive numeric matrix with no NA values")
  }

if(any(colnames(metabarlist$reads) %in% ""))
  stop("metabarlist$reads has empty column names")
if(any(rownames(metabarlist$reads) %in% ""))
  stop("metabarlist$reads has empty row names")
if(any(duplicated(colnames(metabarlist$reads))))
  stop("metabarlist$reads has duplicated column names")
if(any(duplicated(rownames(metabarlist$reads))))
  stop("metabarlist$reads has duplicated row names")


if ( ! (is.data.frame(metabarlist$motus) &&
        is.data.frame(metabarlist$pcrs) &&
        is.data.frame(metabarlist$samples))) {
  stop("metabarlist$motus, metabarlist$pcrs and metabarlist$samples must be data.frames")
}


if(any(colnames(metabarlist$motus) %in% ""))
  stop("metabarlist$motus has empty column names")
if(any(rownames(metabarlist$motus) %in% ""))
  stop("metabarlist$motus has empty row names")
if(any(duplicated(colnames(metabarlist$motus))))
  stop("metabarlist$motus has duplicated column names")
if(any(duplicated(rownames(metabarlist$motus))))
  stop("metabarlist$motus has duplicated row names")


if ( ! ((length(colnames(metabarlist$reads)) == length(rownames(metabarlist$motus))) &&
        (all(colnames(metabarlist$reads) == rownames(metabarlist$motus))))) {
  stop("columns of metabarlist$reads must correspond exactly to rows of metabarlist$motus")
}


if(any(colnames(metabarlist$pcrs) %in% ""))
  stop("metabarlist$pcrs has empty column names")
if(any(rownames(metabarlist$pcrs) %in% ""))
  stop("fmetabarlist$pcrs has empty row names")
if(any(duplicated(colnames(metabarlist$pcrs))))
  stop("metabarlist$pcrs has duplicated column names")
if(any(duplicated(rownames(metabarlist$pcrs))))
  stop("metabarlist$pcrs has duplicated row names")



if ( ! ((length(rownames(metabarlist$reads)) == length(rownames(metabarlist$pcrs))) &&
        (all(rownames(metabarlist$reads) == rownames(metabarlist$pcrs))))) {
  stop("rows of metabarlist$reads must correspond exactly to rows of metabarlist$pcrs")
}

if ( ! ('sequence' %in% colnames(metabarlist$motus) &&
        is.character(metabarlist$motus$sequence) &&
        all(! is.na(metabarlist$motus$sequence)))) {
  stop('metabarlist$motus$sequence must be defined and be exclusively characters')
}

pcrs_mandatory_cols = c('sample_id', 'type','control_type')
if ( ! (all(pcrs_mandatory_cols %in% colnames(metabarlist$pcrs)))) {
  stop("metabarlist$pcrs have mandatory columns: 'sample_id', 'type','control_type'")
}

if ( ! all(sort(unique(metabarlist$pcrs$type), na.last = T) %in% c('control','sample'))) {
  stop("metabarlist$pcrs$type must contain only 'control' and 'sample' values")
}

if ( ! all(is.na(metabarlist$pcrs$control_type) == (metabarlist$pcrs$type=='sample'))) {
  stop("metabarlist$reads$control_type must have 'NA' values for samples")
}

if ( ! all(is.na(metabarlist$pcrs$sample_id) != (metabarlist$pcrss$type=='sample'))) {
  stop("metabarlist$reads$control_type must no have 'NA' values for controls")
}

if ( ! all(sort(unique(metabarlist$pcrs$control_type[!is.na(metabarlist$pcrs$control_type)])) %in% c('extraction', 'pcr',  'positive', 'sequencing'))) {
  stop("metabarlist$pcrs$control_type must be either 'extraction', 'pcr', 'positive' or 'sequencing' for controls")
}


if(any(colnames(metabarlist$samples) %in% ""))
  stop("metabarlist$samples has empty column names")
if(any(rownames(metabarlist$samples) %in% ""))
  stop("metabarlist$samples has empty row names")
if(any(duplicated(colnames(metabarlist$samples)))) stop("metabarlist$samples has duplicated column names")
if(any(duplicated(rownames(metabarlist$samples)))) stop("metabarlist$samples has duplicated row names")


if ( ! (all(unique(metabarlist$pcrs$sample_id[metabarlist$pcrs$type=='sample']) %in% rownames(metabarlist$samples)))) {
  stop('All values in metabarlist$pcrs$sample_id should have a corresponding entry in metabarlist$samples')
}


cols_plate_design = c('tag_fwd', 'tag_rev', 'primer_fwd', 'primer_rev', 'plate_no', 'plate_col', 'plate_row')
if (! all(cols_plate_design %in% colnames(metabarlist$pcrs))) {
  warning(paste0("PCR plate design not properly recorded: ", paste(cols_plate_design[! cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '), " missing !\n"))
}
else {
  if (all(c('plate_no', 'plate_col', 'plate_row') %in% colnames(metabarlist$pcrs))) {
    if ( ! (is.numeric(metabarlist$pcrs$plate_no))) {
      stop("metabarlist$pcrs$plate_no must be numeric")
    }

    if ( ! (all(as.numeric(metabarlist$pcrs$plate_col) %in% 1:12))) {
      stop("metabarlist$pcrs$plate_col must correspond to numbers from 1 to 12")
    }

    if ( ! (all(metabarlist$pcrs$plate_row %in% c('A','B','C','D','E','F','G','H')))) {
      stop("metabarlist$pcrs$plate_row must be letters from A to H")
    }


    if ( ! (nrow(metabarlist$pcrs) == nrow(unique(metabarlist$pcrs[,c('plate_no', 'plate_col', 'plate_row')])))) {
      combi <- table(apply(metabarlist$pcrs[,c('plate_no', 'plate_col', 'plate_row')], MARGIN=1, FUN=function(x) paste(x, collapse=" ")))
      stop(paste0("Same combination of 'plate_no', 'plate_col', 'plate_row' found several times in metabarlist$pcrs (",
                  paste(names(combi)[combi>1],collapse=", "),
                  ")"))
    }
  }
  if (all(c('primer_fwd', 'primer_rev') %in% colnames(metabarlist$pcrs))) {
    if (nrow(unique(metabarlist$pcrs[,c('primer_fwd', 'primer_rev')]))>1) {
      warning('Several primer pairs described in metabarlist$pcrs')
    }
  }
  if (all(c('tag_fwd', 'tag_rev') %in% colnames(metabarlist$pcrs))) {
    if (! nrow(unique(metabarlist$pcrs[,c('tag_fwd', 'tag_rev')]))==nrow(metabarlist$pcrs)) {
      combi <- table(apply(metabarlist$pcrs[,c('tag_fwd', 'tag_rev')], MARGIN=1, FUN=function(x) paste(x, collapse=" ")))

      warning('Several tag pairs are non unique in metabarlist$pcrs (',
              paste(names(combi)[combi>1],collapse=", "),')')
    }
  }
}


if (any(rowSums(metabarlist$reads)==0)) {
  warning("Some pcrs have a number of reads of zero !")
}
if (any(colSums(metabarlist$reads)==0)) {
  warning('Some MOTUs have a number of reads of zero !')
}
TRUE

}
