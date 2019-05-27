#' Plotting a PCR attribute in the tag pair design
#'
#' Plots an attribute of PCRs from a \code{\link{metabarlist}} object in the tag pair design
#'
#'
#' @param metabarlist    a \code{\link{metabarlist}} object
#'  @param FUN            a function which return a vector containing the information to be plotted. The vector should be a numeric vector which has the same length of table 'reads'.
#' @param legend_title   the title of legend containing the plotted information.
#' @param taglist        a character vector corresponding to the full list of tags ordered as they are used in the PCR plate scheme
#'
#' @name ggpcrtag
#'
#' @return a ggplot
#'
#' @details
#' Vizualizing some attributes of the PCRs in their tag pair design context, i.e. according to the tag combination used. This can be useful to identify potential problems (e.g. low amount of reads due to non-functional primer/tag sets).
#'
#' @examples
#'
#' data(soil_euk)
#'
#' #Plot the number of reads per pcrs
#' ggpcrtag(soil_euk)
#'
#' #Plot the number of reads of the most abundant MOTU
#' ggpcrtag(soil_euk, legend_title = "#reads of most\nabundant MOTU",
#'            FUN = function(m){m$reads[,which.max(colSums(m$reads))]})
#'
#'
#' @author Lucie Zinger
#' @importFrom cowplot insert_xaxis_grob insert_yaxis_grob ggdraw
#' @export ggpcrtag

ggpcrtag = function(metabarlist, legend_title = "well_values",
                    FUN = function(metabarlist) {rowSums(metabarlist$reads)},
                    taglist=as.character(unique(soil_euk$pcrs$tag_rev))) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    function_values <- FUN(metabarlist)

    if (length(function_values) != nrow(metabarlist$pcrs)) {
      stop('provided information should have the length of pcrs')
    }

    if(!is.numeric(function_values))
      stop("selected information should be numeric")

    cols_tag_design = c('tag_fwd', 'tag_rev')

    if (!all(cols_tag_design %in% colnames(metabarlist$pcrs)))
      stop("Tag pair design not properly provided: ",
           paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '),
           " missing !\n")

    tag_design = metabarlist$pcrs[,c("tag_fwd", "tag_rev", "control_type")]
    tag_design$control_type = factor(tag_design$control_type,
                                       levels=c("extraction", "pcr", "sequencing", "positive"))

    # add values of parameter's function to the fictive dataframe
    tag_design$well_values <- function_values
    tag_design$well_values[tag_design$well_values == 0] <- NA

    tag_design$tag_fwd = factor(tag_design$tag_fwd,
                                levels = rev(levels(tag_design$tag_fwd)[
                                  match(taglist, levels(tag_design$tag_fwd), nomatch=0)]))

    tag_design$tag_rev = factor(tag_design$tag_rev,
                                levels = levels(tag_design$tag_rev)[
                                  match(taglist, levels(tag_design$tag_rev), nomatch=0)])

   all =
     ggplot(tag_design, aes(y=tag_fwd, x=tag_rev, size=well_values)) +
      geom_raster(aes(fill=control_type), na.rm = TRUE) +
      theme_bw() +
      scale_fill_manual(values=c("brown", "red", "pink", "cyan4"), na.translate = FALSE) +
      geom_point(na.rm = TRUE) +
      labs(size=legend_title) +
      scale_x_discrete(expand=c(0.02,0)) +
      scale_y_discrete(expand=c(0.02,0)) +
      theme(axis.text.x = element_text(angle=45, h=1),
            legend.position = "bottom",
            legend.box = "vertical")

   fwd =
     ggplot(tag_design, aes(x=tag_fwd, y=well_values)) +
       geom_boxplot(lwd=0.5, na.rm = TRUE) +
       theme_bw() +
       stat_summary(fun.y=median, geom="line", aes(group=1), color="cyan2", lwd=1, na.rm = TRUE) +
       labs(y=legend_title, x=NULL) +
       scale_x_discrete(expand=c(0.02,0)) +
       scale_y_discrete(expand=c(0.02,0)) +
       coord_flip() +
       theme(axis.text.y = element_blank())

   rev =
     ggplot(tag_design, aes(x=tag_rev, y=well_values)) +
       geom_boxplot(lwd=0.5, na.rm = TRUE) +
       theme_bw() +
       stat_summary(fun.y=median, geom="line", aes(group=1), color="cyan2", lwd=1, na.rm = TRUE) +
       labs(y=legend_title, x=NULL) +
       scale_x_discrete(expand=c(0.02,0)) +
       scale_y_discrete(expand=c(0.02,0)) +
       theme(axis.text.x = element_blank())

   out = insert_xaxis_grob(all, rev, grid::unit(.2, "null"), position = "top")
   out = insert_yaxis_grob(out, fwd, grid::unit(.2, "null"), position = "right")
   ggdraw(out)
  }
}

