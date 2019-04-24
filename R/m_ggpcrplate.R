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
#' Vizualizing some attributes of the PCRs in their plate design context, i.e. their location in one or a set of 96-well PCR plates. This can be useful to identify potential problems (e.g. low amount of reads due to non-functional primer sets).
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #Plot the number of reads per pcrs
#'
#' ##Store the number of reads per pcrs in the pcrs table
#' soil_euk$pcrs$seq_depth = rowSums(soil_euk$reads)
#'
#' ##Plot the results
#' ggpcrplate(soil_euk, table = "pcrs", index = "seq_depth")
#'
#'
#' #Plot the number of reads of the most abundant MOTU
#'
#' ##Get the name of the most abundant MOTU
#' id = colnames(soil_euk$reads)[which.max(colSums(soil_euk$reads))]
#'
#' ##Plot the results
#' ggpcrplate(soil_euk, table = "reads", index = id) + labs(size="#reads of most\nabundant MOTU")
#'
#'
#' @author Lucie Zinger
#' @import ggplot2
#' @export ggpcrplate

ggpcrplate = function(metabarlist, table, index) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    cols_plate_design = c('plate_no', 'plate_col', 'plate_row')

    if (!all(cols_plate_design %in% colnames(metabarlist$pcrs)))
      stop("PCR plate design not properly provided: ",
           paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '),
           " missing !\n")

    extract_table_methods = c("reads", "pcrs")
    tab = match.arg(table, extract_table_methods)

    if(!index %in% colnames(metabarlist[[tab]]))
      stop("index is not in table")

    if(!is.numeric(metabarlist[[tab]][,index]))
      stop("selected information should be numeric")

    plate_design = metabarlist$pcrs[,c("plate_no", "plate_col", "plate_row", "control_type")]
    plate_design$control_type = factor(plate_design$control_type,
                                       levels=c("extraction", "pcr", "sequencing", "positive", NA))

    plate_design[index] = metabarlist[[tab]][,index]

    plate_design[index][plate_design[index]==0] = NA

    ggplot(plate_design, aes(y=match(plate_row, LETTERS[1:8]), x=plate_col, size=get(index))) +
      geom_raster(aes(fill=control_type)) +
      facet_wrap(~plate_no, scale="free") + theme_bw() +
      scale_y_reverse(breaks = 1:8, labels=LETTERS[1:8]) +
      scale_x_continuous(breaks = 1:12) +
      scale_fill_manual(values=c("brown", "red", "pink", "cyan4", "white")) +
      geom_point() +
      labs(x=NULL, y=NULL, size=index)

  }
}




