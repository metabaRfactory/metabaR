#' Plotting a PCR attribute in the PCR plate design
#'
#' Plots an attribute of PCRs from a \code{\link{metabarlist}} object into the corresponding 96-well PCR plate design
#'
#'
#' @param metabarlist    a \code{\link{metabarlist}} object.
#' @param FUN            a function which returns a vector containing the information to be plotted.
#'                       The vector should be numeric and of equal length to the number of
#'                       rows of the `reads` table.
#' @param legend_title   the title of legend containing the plotted information.
#'
#' @name ggpcrplate
#'
#' @return a ggplot
#'
#' @details
#' Visualising PCR attributes in their plate design context, i.e. their location in one or a set of 96-well PCR plates. This can be useful for identifying potential problems (e.g. high amount of reads in one control due to cross contaminations with neighboring samples).
#'
#' @seealso \code{\link{ggpcrtags}}
#'
#' @examples
#'
#' data(soil_euk)
#'
#' ## Plot the number of reads per pcr
#' ggpcrplate(soil_euk)
#'
#'
#' ## Plot the number of reads of the most abundant MOTU
#' ggpcrplate(soil_euk,
#'   legend_title = "#reads of most \nabundant MOTU",
#'   FUN = function(m) {
#'     m$reads[, which.max(colSums(m$reads))]
#'   }
#' )
#' @author Lucie Zinger, Clément Lionnet
#' @import ggplot2
#' @export ggpcrplate


ggpcrplate <- function(metabarlist, legend_title = "well_values",
                       FUN = function(metabarlist) {rowSums(metabarlist$reads)}) {
  if (suppressWarnings(check_metabarlist(metabarlist))) {
    function_values <- FUN(metabarlist)

    if (length(function_values) != nrow(metabarlist$pcrs)) {
      stop("provided information should have the length of pcrs")
    }

    if (!is.numeric(function_values)) {
      stop("provided information should be numeric")
    }

    cols_plate_design <- c("plate_no", "plate_col", "plate_row")

    if (!all(cols_plate_design %in% colnames(metabarlist$pcrs))) {
      stop(
        "PCR plate design not properly provided: ",
        paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)], sep = ", "),
        " missing !\n"
      )
    }

    # create a fictive dataframe to plot data with ggplot
    plate_design <- metabarlist$pcrs[, c("plate_no", "plate_col", "plate_row", "control_type")]
    # define factor level order for the legend
    plate_design_levels <- c("extraction", "pcr", "sequencing", "positive")
    plate_design$control_type <- factor(plate_design$control_type, levels = plate_design_levels)

    # add values of parameter's function to the fictive dataframe
    plate_design$well_values <- function_values
    plate_design$well_values[plate_design$well_values == 0] <- NA

    ggplot(plate_design, aes(y = match(plate_row, LETTERS[1:8]), x = plate_col, size = well_values)) +
      geom_raster(aes(fill = control_type), na.rm = TRUE) +
      facet_wrap(~plate_no, scale = "free") + theme_bw() +
      scale_y_reverse(breaks = 1:8, labels = LETTERS[1:8]) +
      scale_x_continuous(breaks = 1:12) +
      scale_fill_manual(values = c("brown", "red", "pink", "cyan4"), na.translate = FALSE) +
      geom_point(na.rm = TRUE) +
      labs(x = NULL, y = NULL, size = legend_title)
  }
}
