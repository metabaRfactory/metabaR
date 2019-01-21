#' Plotting an amplicon attribute in the PCR plate design
#'
#' Plots an attribute of amplicons from a \code{\link{TODEFINE}} object in the correspondint 96-well PCR plate design
#'
#'
#' @param attr    a numeric vector of the amplicon's attribute to be plotted
#' @param plate_no  a numeric vector of the plate id(s) where the amplicons were amplified
#' @param plate_col  a numeric vector of the plate's column where the amplicons were amplified
#' @param plate_row  a character vector of letters corresponding to the plate's row where the amplicons were amplified
#' @param control_type  a factor indicating the amplicon type (controls type, sample types, etc.)
#'
#' @name ggplate
#'
#' @return a ggplot
#'
#' @details
#' Vizualizing some attributes of the amplicon in its experimental context, here its location in one or a set of 96-well PCR plates, can be useful to identify potential set of samples that exhibit particular features (e.g. low amount of reads).
#'
#' @examples
#'
#' data(soil_euk)
#' #plot number of reads per pcrs
#' seqdepth = rowSums(soil_euk$reads)
#' p = ggplate(attr = seqdepth, plate_no = soil_euk$pcrs$plate_no,
#'                   plate_col = soil_euk$pcrs$plate_col,
#'                   plate_row =  soil_euk$pcrs$plate_row,
#'                   control_type = soil_euk$pcrs$Control_type)
#' p + labs(size="# reads", fill="control type") +
#'   scale_fill_manual(values=c("brown", "pink", "cyan4", "red"),
#'                     na.value = "white")
#'
#' @author Lucie Zinger
#' @import ggplot2
#' @export ggplate


ggplate = function(attr, plate_no, plate_col, plate_row, control_type) {

  if (is.numeric(attr)==F)
    stop("attr should be numeric")
  if (length(grep("FALSE", as.numeric(plate_col) %in% 1:12)))
    stop("plate_col should range between 1 and 12")
  if (length(grep("FALSE", plate_row %in% LETTERS[1:8])))
    stop("plate_row should range between A and H")

  dat = data.frame(plate_no = plate_no, plate_col = as.numeric(plate_col),
                   plate_row = plate_row, attr=attr, control_type = control_type)

  ggplot(dat, aes(y=match(plate_row, LETTERS[1:8]), x=plate_col, size=attr)) +
    geom_raster(aes(fill=control_type)) +
    facet_wrap(~paste("Plate", plate_no), scale="free") + theme_bw() +
    scale_y_reverse(breaks = 1:8, labels=LETTERS[1:8]) +
    scale_x_continuous(breaks = 1:12) +
    geom_point() +
    labs(x=NULL, y=NULL)
}



