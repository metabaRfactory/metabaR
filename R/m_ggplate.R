#' Plotting an attribute in the PCR plate design
#'
#' Plots an attribute of amplicons from a \code{\link{TODEFINE}} object in the correspondint PCR plate design
#'
#'
#' @param x    a \code{\link{TODEFINE}} object
#'
#' @name ggplate
#'
#' @return a ggplot
#'
#' @details
#' TO WRITE
#'
#' @examples
#'
#' data(soil_euk)
#' #plot number of reads per pcrs
#' seqdepth = rowSums(soil_euk$reads)
#' ggplate(attr = seqdepth, plate_no = soil_euk$pcrs$plate_no,
#'                   plate_col = soil_euk$pcrs$plate_col,
#'                   plate_row =  soil_euk$pcrs$plate_row,
#'                   control_type = soil_euk$pcrs$Control_type)
#'
#' @author Fred Boyer, Lucie Zinger
#' @import ggplot2
#' @export ggplate


ggplate = function(attr, plate_no, plate_col, plate_row, type) {

  dat = data.frame(plate_no = plate_no, plate_col = plate_col,
                   plate_row = plate_row, attr=attr, type=type)

  ggplot(dat, aes(y=plate_row, x=plate_col, size=attr)) +
    geom_raster(aes(fill=type)) +
    facet_wrap(~plate_no, scale="free") + theme_bw() +
    scale_fill_hue(na.value="white") +
    geom_point() +
}



