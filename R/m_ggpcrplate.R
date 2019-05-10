#' Plotting a PCR attribute in the PCR plate design
#'
#' Plots an attribute of PCRs from a \code{\link{metabarlist}} object in the correspondint 96-well PCR plates design
#'
#'
#' @param metabarlist    a \code{\link{metabarlist}} object
#' @param table          the table where the information on which the plotting is based. Can be only `reads`,  or `pcrs`.
#' @param index          the name of the column name containing the information to be plotted. This information should be a numeric vector.
#'
#' @name ggpcrplate
#'
#' @return a ggplot
#'
#' @details
#' Vizualizing some attributes of the PCRs in their plate design context, i.e. their location in one or a set of 96-well PCR plates. This can be useful to identify potential problems (e.g. high amount of reads in one control due to cross contaminations with neighboring samples).
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #Plot the number of reads per pcrs
#' ggpcrplate(soil_euk, FUN=function(m){rowSums(m$reads)})
#'
#'
#' #Plot the number of reads of the most abundant MOTU
#' library(ggplot2)
#' ggpcrplate(soil_euk, FUN=function(m){m$reads[,which.max(colSums(m$reads))]})) +
#' labs(size="#reads of most\nabundant MOTU")
#'
#'
#' @author Lucie Zinger
#' @import ggplot2
#' @export ggpcrplate


ggpcrplate = function(metabarlist, legend_title="well_values", FUN = function(metabarlist) {rowSums(metabarlist$reads)}) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    function_values <- FUN(metabarlist)

    if (length(function_values) != nrow(metabarlist$pcrs)) {
      stop('provided information should have the length of pcrs')
    }

    if(!is.numeric(function_values))
      stop("provided information should be numeric")

    cols_plate_design <- c('plate_no', 'plate_col', 'plate_row')

    if (!all(cols_plate_design %in% colnames(metabarlist$pcrs)))
      stop("PCR plate design not properly provided: ",
           paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '),
           " missing !\n")

    plate_design <- metabarlist$pcrs[,c("plate_no", "plate_col", "plate_row", "control_type")]
    plate_design_levels <- c(levels(plate_design$control_type), "sample")
    plate_design$control_type <- factor(plate_design$control_type, levels=plate_design_levels)
    plate_design$control_type[is.na(plate_design$control_type)] <- "sample"
    plate_design$control_type <- factor(plate_design$control_type,
                                       levels=c("extraction", "pcr", "sequencing", "positive", "sample"))

    plate_design$well_values <- function_values
    plate_design$well_values[plate_design$well_values==0] <- NA

    ggplot(plate_design, aes(y=match(plate_row, LETTERS[1:8]), x=plate_col, size=well_values)) +
      geom_raster(aes(fill=control_type)) +
      facet_wrap(~plate_no, scale="free") + theme_bw() +
      scale_y_reverse(breaks = 1:8, labels=LETTERS[1:8]) +
      scale_x_continuous(breaks = 1:12) +
      scale_fill_manual(values=c("brown", "red", "pink", "cyan4", "white")) +
      geom_point(na.rm=TRUE) +
      labs(x=NULL, y=NULL, size=legend_title)
  }
}
