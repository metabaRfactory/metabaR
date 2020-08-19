#' Checks if a list of tables is a \code{\link{metabarlist}}
#'
#' Tests if a list of tables contains correct information on MOTU abundances, MOTU characteristics, PCR characteristics, and sample characteristics to be classed as a \code{\link{metabarlist}} object.
#'
#' @param metabarlist a \code{\link{metabarlist}} object
#'
#' @name check_metabarlist
#'
#' @return \code{check_metabarlist} returns TRUE if \code{metabarlist} fulfills of conditions to be of class \code{\link{metabarlist}}` (see details). If not, the function throws a stop.
#'
#' @details
#'
#' The function \code{check_metabarlist} checks if \code{metabarlist} has all the properties required to be a well formed \code{metabarlist} object:
#'
#' \itemize{
#' \item {\code{metabarlist} should be a list with four tables named
#'        `reads`, `motus`, `pcrs` and `samples`.}
#' \item {the `reads` table should be a numeric matrix.}
#' \item {the `motus`, `pcrs` and `samples` tables should be data frames.}
#' \item {the `reads` table row names should be exactly the same as row names
#'        in the `pcrs` table.}
#' \item {the `reads` table column names should be exactly the same as row names
#'        in the `motus` table.}
#' \item {the `pcrs` table should contain the mandatory columns required,
#'        i.e.  `sample_id`, `type` and `control_type`.}
#' \item {values in `type` should be properly defined:
#'        they can only be `sample` or `control`.}
#' \item {values in `control_type` should be properly defined:
#'        they can only be `sequencing`, `pcr`, `extraction`, `positive`, (or NA for samples).}
#' }
#'
#' The function will stop if these basic criteria are not met.
#'
#' \code{check_metabarlist} also issues warnings if tables lack columns which are non-mandatory for the \code{metabaR} package to run, but that are mandatory (or have a mandatory format) for specific functions (e.g., \code{\link{ggpcrplate}})
#' \itemize{
#' \item {the column `sequence` for the `motus` table.}
#' \item {the columns `tag_fwd`, `tag_rev`, `primer_fwd`, `primer_rev`,
#'        `plate_no`, `plate_col`, and `plate_row` for the `pcrs` table.}
#' }
#'
#'
#' In addition, the function issues a warning if any PCR or MOTU has a total count of 0 across the entire dataset.
#'
#' @seealso \code{\link{metabarlist_generator}}
#'
#' @examples
#'
#' data(soil_euk)
#' check_metabarlist(soil_euk)
#'
#' @author Clément Lionnet & Frédéric Boyer & Lucie Zinger
#' @export check_metabarlist


check_metabarlist <- function(metabarlist) {

  metabarlist.name = substitute(metabarlist)

  if (!"metabarlist" %in% class(metabarlist)) {
    stop(paste("class of", metabarlist.name, "is not 'metabarlist'"))
  }

  if (!is.list(metabarlist)) {
    stop(paste(metabarlist.name, "is not a list"))
  }

  slots <- c("reads", "motus", "pcrs", "samples")
  if (!all(slots %in% names(metabarlist))) {
    stop(paste0(
      metabarlist.name,
      " does not contain the following objects: ",
      "'",
      paste(slots[!slots %in% names(metabarlist)], collapse = "', '"),
      "'"
    ))
  }

  if (!(
    is.matrix(metabarlist$reads) &&
    is.numeric(metabarlist$reads) &&
    !any(is.na(metabarlist$reads)) &&
         all(metabarlist$reads >= 0)
  )) {
    stop(
      paste(
        "table `reads` in",
        metabarlist.name,
        "must be a positive numeric matrix with no NA values"
      )
    )
  }

  if (any(colnames(metabarlist$reads) %in% "")) {
    stop(paste("table `reads` in", metabarlist.name, "has empty column names"))
  }

  if (any(rownames(metabarlist$reads) %in% "")) {
    stop(paste("table `reads` in", metabarlist.name, "has empty row names"))
  }

  if (any(duplicated(colnames(metabarlist$reads)))) {
    stop(paste(
      "table `reads` in",
      metabarlist.name,
      "has duplicated column names"
    ))
  }

  if (any(duplicated(rownames(metabarlist$reads)))) {
    stop(paste(
      "table `reads` in",
      metabarlist.name,
      "has duplicated row names"
    ))
  }

  if (!(
    is.data.frame(metabarlist$motus) &&
    is.data.frame(metabarlist$pcrs) &&
    is.data.frame(metabarlist$samples)
  )) {
    stop(
      paste(
        "tables `motus`, `pcrs` and `samples` in",
        metabarlist.name,
        "must be data.frames"
      )
    )
  }

  if (any(colnames(metabarlist$motus) %in% "")) {
    stop(paste("table `motus` in", metabarlist.name, "has empty column names"))
  }

  if (any(rownames(metabarlist$motus) %in% "")) {
    stop(paste("table `motus` in", metabarlist.name, "has empty row names"))
  }

  if (any(duplicated(colnames(metabarlist$motus)))) {
    stop(paste(
      "table `motus` in",
      metabarlist.name,
      "has duplicated column names"
    ))
  }

  if (any(duplicated(rownames(metabarlist$motus)))) {
    stop(paste(
      "table `motus` in",
      metabarlist.name,
      "has duplicated row names"
    ))
  }

  if (!((length(colnames(metabarlist$reads)) == length(rownames(metabarlist$motus))) &&
        (all(
          colnames(metabarlist$reads) == rownames(metabarlist$motus)
        )))) {
    stop(
      paste(
        "columns of table `reads` in",
        metabarlist.name ,
        "must correspond exactly to rows of `motus`"
      )
    )
  }

  if (any(colnames(metabarlist$pcrs) %in% "")) {
    stop(paste("table `pcrs` in", metabarlist.name, "has empty column names"))
  }

  if (any(rownames(metabarlist$pcrs) %in% "")) {
    stop(paste("table `pcrs` in", metabarlist.name, "has empty row names"))
  }

  if (any(duplicated(colnames(metabarlist$pcrs)))) {
    stop(paste(
      "table `pcrs` in",
      metabarlist.name,
      "has duplicated column names"
    ))
  }

  if (any(duplicated(rownames(metabarlist$pcrs)))) {
    stop(paste("table `pcrs` in", metabarlist.name, "has duplicated row names"))
  }

  if (!((length(rownames(metabarlist$reads)) == length(rownames(metabarlist$pcrs))) &&
        (all(
          rownames(metabarlist$reads) == rownames(metabarlist$pcrs)
        )))) {
    stop(
      paste(
        "rows of table `reads` in",
        metabarlist.name,
        "must correspond exactly to rows of `pcrs`"
      )
    )
  }

  if (!(
    "sequence" %in% colnames(metabarlist$motus) &&
    is.character(metabarlist$motus$sequence) &&
    all(!is.na(metabarlist$motus$sequence))
  )) {
    stop(
      paste(
        "column `sequence` of table `motus` in",
        metabarlist.name,
        "must be defined and be exclusively characters"
      )
    )
  }

  pcrs_mandatory_cols <- c("sample_id", "type", "control_type")
  if (!(all(pcrs_mandatory_cols %in% colnames(metabarlist$pcrs)))) {
    stop(
      paste(
        "table `pcrs` in",
        metabarlist.name,
        "has mandatory columns: `sample_id`, `type`,`control_type`"
      )
    )
  }

  if (!all(sort(unique(metabarlist$pcrs$type), na.last = T) %in% c("control", "sample"))) {
    stop(
      paste(
        "column `type` of table `pcrs` in",
        metabarlist.name,
        "must contain only `control` or `sample` values"
      )
    )
  }

  if (!all(is.na(metabarlist$pcrs$control_type) == (metabarlist$pcrs$type == "sample"))) {
    stop(
      paste(
        "column `control_type` of table `pcrs` in",
        metabarlist.name,
        "must have 'NA' values for samples"
      )
    )
  }

  if (!all(is.na(metabarlist$pcrs$sample_id) != (metabarlist$pcrss$type == "sample"))) {
    stop(
      paste(
        "'NA' values are not allowed for controls in column `control_type` of table `pcrs` in",
        metabarlist.name
      )
    )
  }

  if (!all(sort(unique(metabarlist$pcrs$control_type[!is.na(metabarlist$pcrs$control_type)])) %in% c("extraction", "pcr", "positive", "sequencing"))) {
    stop(
      paste(
        "column `control_type` of table `pcrs` in",
        metabarlist.name,
        "must be either `extraction`, `pcr`, `positive` or `sequencing` for controls"
      )
    )
  }

  if (any(colnames(metabarlist$samples) %in% "")) {
    stop(paste(
      "table `samples` in",
      metabarlist.name,
      "has empty column names"
    ))
  }

  if (any(rownames(metabarlist$samples) %in% "")) {
    stop(paste("table `samples` in", metabarlist.name, "has empty row names"))
  }

  if (any(duplicated(colnames(metabarlist$samples)))) {
    stop(paste(
      "table `samples` in",
      metabarlist.name,
      "has duplicated column names"
    ))
  }

  if (any(duplicated(rownames(metabarlist$samples)))) {
    stop(paste(
      "table `samples` in",
      metabarlist.name,
      "has duplicated row names"
    ))
  }

  if (!(all(
    unique(metabarlist$pcrs$sample_id[metabarlist$pcrs$type == "sample"]) %in%
    rownames(metabarlist$samples)
  ))) {
    v <-
      unique(metabarlist$pcrs$sample_id[metabarlist$pcrs$type == "sample"])
    message(
      paste(
        v[!v %in% rownames(metabarlist$samples)],
        "from column `sample_id` of table `pcrs` in ",
        metabarlist.name,
        "are not found in the row names of table `samples`",
        collapse = ""
      )
    )
    stop(
      paste(
        "All values in column `sample_id` of table `pcrs` in",
        metabarlist.name,
        "should have a corresponding entry in table `samples`"
      )
    )
  }

  cols_plate_design <- c(
    "tag_fwd",
    "tag_rev",
    "primer_fwd",
    "primer_rev",
    "plate_no",
    "plate_col",
    "plate_row"
  )

  if (!all(cols_plate_design %in% colnames(metabarlist$pcrs))) {
    warning(
      paste0(
        "PCR plate design not properly provided: columns ",
        paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)],
              sep = ", "),
        " is missing in table `pcrs` of ",
        metabarlist.name,
        "!\n"
      )
    )
  } else {
    if (all(c("plate_no", "plate_col", "plate_row") %in% colnames(metabarlist$pcrs))) {
      if (!(is.numeric(metabarlist$pcrs$plate_no))) {
        warning(paste(
          "column `plate_no` of table `pcrs` in",
          metabarlist.name,
          "is not numeric"
        ))
      }

      if (!(all(as.numeric(metabarlist$pcrs$plate_col) %in% 1:12))) {
        stop(
          paste(
            "column `plate_col` of table `pcrs` in",
            metabarlist.name,
            "must contain numeric values ranging from 1 to 12"
          )
        )
      }

      if (!(all(
        metabarlist$pcrs$plate_row %in% c("A", "B", "C", "D", "E", "F", "G", "H")
      ))) {
        stop(
          paste(
            "column `plate_row` of table `pcrs` in",
            metabarlist.name,
            "must be letters ranging from A to H"
          )
        )
      }

      if (!(nrow(metabarlist$pcrs) ==
            nrow(unique(metabarlist$pcrs[, c("plate_no", "plate_col", "plate_row")])))) {
        combi <-
          table(apply(
            metabarlist$pcrs[, c("plate_no", "plate_col", "plate_row")],
            MARGIN = 1,
            FUN = function(x)
              paste(x, collapse = " ")
          ))
        stop(
          paste0(
            "Same combination(s) of 'plate_no', 'plate_col', 'plate_row' found
                  several times in table `pcrs` of ",
            metabarlist.name,
            " :",
            paste(names(combi)[combi > 1], collapse = ", ")
          )
        )
      }
    }

    if (all(c("primer_fwd", "primer_rev") %in% colnames(metabarlist$pcrs))) {
      if (nrow(unique(metabarlist$pcrs[, c("primer_fwd", "primer_rev")])) > 1) {
        warning(
          paste(
            "columns `primer_fwd` or `primer_rev` from table `pcrs` in",
            metabarlist.name,
            "contain several primer pairs"
          )
        )
      }
    }

    if (all(c("tag_fwd", "tag_rev") %in% colnames(metabarlist$pcrs))) {
      if (!nrow(unique(metabarlist$pcrs[, c("tag_fwd", "tag_rev")])) ==
            nrow(metabarlist$pcrs)) {
        combi <- table(apply(
          metabarlist$pcrs[, c("tag_fwd", "tag_rev")],
          MARGIN = 1,
          FUN = function(x)
            paste(x, collapse = " ")
        ))

        warning(
          paste(
            "columns `tag_fwd` and `tag_rev` from table `pcrs` in",
            metabarlist.name,
            "contain one or several tag pairs that are not unique:",
            paste(names(combi)[combi > 1], collapse = ", ")
          )
        )
      }
    }
  }

  if (any(rowSums(metabarlist$reads) == 0)) {
    warning(
      paste(
        "Some PCRs in",
        metabarlist.name,
        "have a number of reads of zero in table `reads`!"
      )
    )
  }

  if (any(colSums(metabarlist$reads) == 0)) {
    warning(
      paste(
        "Some MOTUsin",
        metabarlist.name,
        "have a number of reads of zero in table `reads`!"
      )
    )
  }

  TRUE
}
