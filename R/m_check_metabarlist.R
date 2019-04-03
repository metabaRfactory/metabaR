#' Test that the provided list is a well formed metabaRffe list
#'
#'
#' @param metabarlist a \code{metabarlist} object
#'
#' @name check_metabarlist
#'
#' @return TRUE or throws an stop
#'
#' @details
#'
#' Check for the properties awaited for a well formed metabaRffe list:
#' \itemize{
#' \item {is a list with three attributes named `reads`, `motus` and `pcrs`}
#' \item {`reads` is a numeric matrix}
#' \item {`motus`, `pcrs` and `samples` are of type `data.frame`}
#' \item {`reads` have the exact same row names and columns as rownames of `motus` and `pcrs` (same order also)}
#' \item {`pcrs` have mandatory columns, i.e.  `sample_id`, `type` and `control_type`}
#' \item {values in `type` are properly defined, i.e. `sample` or `control`}
#' \item {values in and `Control_type` are properly defined, i.e. `sequencing`, `pcr`, `extraction`, `positive`}
#'
#'
#' \item {issue warnings if}
#' \item {`motus` and `pcrs` have not mandatory columns for design description, i.e. tag_fwd, tag_rev, primer_fwd, primer_rev, plate_no, plate_col and plate_row for `pcrs`
#'   and sequence for `motus`}
#' }
#'
#' \item {issue an stop if these informations are not well recorded}
#'
#' Moreover, the function issues a warning if any sample or OTU is associated to a count of 0
#'
#' @examples
#'
#' data(soil_euk)
#'
#' check_metabarlist(soil_euk)
#'
#' @author Clément Lionnet & Frédéric Boyer
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


  if ( ! (is.matrix(metabarlist$reads) && is.numeric(metabarlist$reads) && ! any(is.na(metabarlist$reads) && all(metabarlist$reads>=0)))) {
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


if ( ! (is.data.frame(metabarlist$motus) && is.data.frame(metabarlist$pcrs) && is.data.frame(metabarlist$samples))) {
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


if ( ! ((length(rownames(metabarlist$reads)) == length(rownames(metabarlist$motus))) && (all(rownames(metabarlist$reads) == rownames(metabarlist$motus))))) {
  stop("lines of metabarlist$reads must correspond exactly to lines of metabarlist$motus")
}

if(any(colnames(metabarlist$pcrs) %in% ""))
  stop("metabarlist$pcrs has empty column names")
if(any(rownames(metabarlist$pcrs) %in% ""))
  stop("fmetabarlist$pcrs has empty row names")
if(any(duplicated(colnames(metabarlist$pcrs))))
  stop("metabarlist$pcrs has duplicated column names")
if(any(duplicated(rownames(metabarlist$pcrs))))
  stop("metabarlist$pcrs has duplicated row names")




if ( ! ((length(colnames(metabarlist$reads)) == length(rownames(metabarlist$pcrs)))  && (all(colnames(metabarlist$reads) == rownames(metabarlist$pcrs))))) {
  stop("columns of metabarlist$reads must correspond exactly to lines of metabarlist$pcrs")
}

if ( ! ('sequence' %in% colnames(metabarlist$motus) && is.character(metabarlist$motus$sequence) && all(! is.na(metabarlist$motus$sequence)))) {
  stop('metabarlist$motus$sequence must be defined and be exclusively character values')
}

pcrs_mandatory_cols = c('sample_id', 'type','control_type')
if ( ! (all(pcrs_mandatory_cols %in% colnames(metabarlist$pcrs)))) {
  stop("metabarlist$pcrs have mandatory columns: 'sample_id', 'type','control_type'")
}

if ( ! (sort(unique(metabarlist$pcrs$Type), na.last = T) == c('control','sample'))) {
  stop("metabarlist$pcrs$Type must contain only 'control' and 'sample' values")
}

if ( ! (is.na(metabarlist$reads$Control_type) == (metabarlist$reads$Type=='sample'))) {
  stop("metabarlist$reads$Control_type must have 'NA' values for samples")
}

if ( ! (is.na(metabarlist$reads$Sample_id) != (metabarlist$reads$Type=='sample'))) {
  stop("metabarlist$reads$Control_type must no have 'NA' values for controls")
}

if ( ! (sort(unique(metabarlist$reads$Control_type[!is.na(metabarlist$reads$Control_type)])) == c('extraction', 'pcr',  'positive', 'sequencing'))) {
  stop("metabarlist$reads$Control_type must no be either 'extraction', 'pcr',  'positive' or 'sequencing' for controls")
}


if(any(colnames(metabarlist$samples) %in% ""))
  stop("metabarlist$samples has empty column names")
if(any(rownames(metabarlist$samples) %in% ""))
  stop("metabarlist$samples has empty row names")
if(any(duplicated(colnames(metabarlist$samples)))) stop("metabarlist$samples has duplicated column names")
if(any(duplicated(rownames(metabarlist$samples)))) stop("metabarlist$samples has duplicated row names")



if ( ! (all(unique(metabarlist$reads$Sample_id[!is.na(metabarlist$reads$Sample_id)]) == rownames(metabarlist$samples)))) {
  stop('All values in metabarlist$reads$Sample_id should have a corresponding entry in metabarlist$samples')
}


cols_plate_design = c('tag_fwd', 'tag_rev', 'primer_fwd', 'primer_rev', 'plate_no', 'plate_col', 'plate_row')
if (! all(cols_plate_design %in% colnames(metabarlist$pcrs))) {
  warning(paste0("No properly recorded plate design: ", paste(cols_plate_design[! cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '), " missing !"))
}
else {
  if (all(c('plate_no', 'plate_col', 'plate_row') %in% metabarlist$pcrs)) {
    if ( ! (is.numeric(metabarlist$pcrs$plate_no))) {
      stop("metabarlist$pcrs$plate_no must be numeric")
    }

    if ( ! (all(as.numeric(metabarlist$pcrs$plate_col) %in% 1:12))) {
      stop("metabarlist$pcrs$plate_col must correspond to numbers to 1:12")
    }

    if ( ! (all(metabarlist$pcrs$plate_row) %in% c('A','B','C','D','E','F','G','H'))) {
      stop("metabarlist$pcrs$plate_row must be letters from A to H")
    }


    if ( ! (nrow(metabarlist$pcrs) == nrow(unique(metabarlist$pcrs[,c('plate_no', 'plate_col', 'plate_row')])))) {
      combi <- table(apply(metabarlist$pcrs[,c('plate_no', 'plate_col', 'plate_row')], MARGIN=1, FUN=function(x) paste(x, collapse=" ")))
      stop(paste0("the same combination of 'plate_no', 'plate_col', 'plate_row' found several times in metabarlist$pcrs (",
                  paste(names(combi)[combi>1],collapse=", "),
                  ")"))
    }
  }
  if (all(c('primer_fwd', 'primer_rev') %in% metabarlist$pcrs)) {
    if (nrow(unique(metabarlist$pcrs[,c('primer_fwd', 'primer_rev')]))>1) {
      warning('Several primer pairs described in metabarlist$pcrs')
    }
  }
  if (all(c('tag_fwd', 'tag_rev') %in% metabarlist$pcrs)) {
    if (! nrow(unique(metabarlist$pcrs[,c('tag_fwd', 'tag_rev')]))==nrow(metabarlist$pcrs)) {
      combi <- table(apply(metabarlist$pcrs[,c('tag_fwd', 'tag_rev')], MARGIN=1, FUN=function(x) paste(x, collapse=" ")))

      warning('Several tag pairs are non unique in metabarlist$pcrs (',
              paste(names(combi)[combi>1],collapse=", "),')')
    }
  }
}

if (any(rowSums(metabarlist$reads)==0)) {
  warning("Some samples have a count of zero !")
}
if (any(colSums(metabarlist$reads)==0)) {
  warning('Some OTUs have a count of zero !')
}
TRUE

}
