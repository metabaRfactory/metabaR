#' Plotting a PCR attribute in the tag pair design
#'
#' Plots an attribute of PCRs from a \code{\link{metabarlist}} object in the tag pair design
#'
#'
#' @param metabarlist    a \code{\link{metabarlist}} object
#' @param table          the table containing the information used for plotting. Can be either `reads`,  or `pcrs`.
#' @param index          the column name of the column containing the information to be plotted. This information should be a numeric vector.
#' @param taglist        a character vector corresponding to the full list of tags ordered as they are used in the PCR plate scheme
#'
#' @name ggpcrtag
#'
#' @return a ggplot
#'
#' @details
#' Vizualizing attributes of the PCRs in their tag pair design context, i.e. according to the tag combination used. This can be useful to identify potential problems (e.g. low amount of reads due to non-functional primer/tag sets).
#'
#' @examples
#'
#' data(soil_euk)
#'
#' taglist = c("acacacac", "acagcaca", "gtgtacat", "tatgtcag", "tagtcgca",
#'             "tactatac", "actagatc", "gatcgcga", "cgctctcg", "gtcgtaga",
#'             "gtcacgtc", "gactgatg", "agactatg", "gcgtcagc", "tgacatca",
#'             "acatgtgt", "gtacgact", "atgatcgc", "acgacgag", "catcagtc",
#'             "atcagtca", "tctactga", "gatgatct", "ctgcgtac", "agcgacta",
#'             "tcagtgtc", "actctgct", "atatagcg", "ctatgcta", "tcgcgctg",
#'             "agcacagt", "tagctagt", "agtgctac", "cgtataca", "cgagtcgt",
#'             "cacatgat")
#'
#' #Plot the number of reads per pcr
#'
#' ##Store the number of reads per pcr in the pcrs table
#' soil_euk$pcrs$seq_depth = rowSums(soil_euk$reads)
#'
#' ##Plot the results
#' ggpcrtag(soil_euk, table = "pcrs", index = "seq_depth", taglist = taglist)
#'
#'
#' #Plot the number of reads of the most abundant MOTU
#'
#' ##Get the name of the most abundant MOTU
#' id = colnames(soil_euk$reads)[which.max(colSums(soil_euk$reads))]
#'
#' ##Plot the results
#' ggpcrtag(soil_euk, table = "reads", index = id, taglist = taglist)
#'
#'
#' @author Lucie Zinger
#' @importFrom cowplot insert_xaxis_grob insert_yaxis_grob ggdraw
#' @export ggpcrtag

ggpcrtag = function(metabarlist, table, index, taglist) {

  if(suppressWarnings(check_metabarlist(metabarlist))) {

    cols_tag_design = c('tag_fwd', 'tag_rev')

    if (!all(cols_tag_design %in% colnames(metabarlist$pcrs)))
      stop("Tag pair design not properly provided: ",
           paste(cols_plate_design[!cols_plate_design %in% colnames(metabarlist$pcrs)], sep=', '),
           " missing !\n")

    extract_table_methods = c("reads", "pcrs")
    tab = match.arg(table, extract_table_methods)

    if(!index %in% colnames(metabarlist[[tab]]))
      stop("index is not in table")

    if(!is.numeric(metabarlist[[tab]][,index]))
      stop("selected information should be numeric")

    tag_design = metabarlist$pcrs[,c("tag_fwd", "tag_rev", "control_type")]
    tag_design$control_type = factor(tag_design$control_type,
                                       levels=c("extraction", "pcr", "sequencing", "positive", NA))

    tag_design[index] = metabarlist[[tab]][,index]

    tag_design[index][tag_design[index]==0] = NA

    tag_design$tag_fwd = factor(tag_design$tag_fwd,
                                levels = rev(levels(tag_design$tag_fwd)[
                                  match(taglist, levels(tag_design$tag_fwd), nomatch=0)]))

    tag_design$tag_rev = factor(tag_design$tag_rev,
                                levels = levels(tag_design$tag_rev)[
                                  match(taglist, levels(tag_design$tag_rev), nomatch=0)])

   all =
     ggplot(tag_design, aes(y=tag_fwd, x=tag_rev, size=get(index))) +
      geom_raster(aes(fill=control_type)) +
      theme_bw() +
      scale_fill_manual(values=c("brown", "red", "pink", "cyan4", "white")) +
      geom_point() +
      labs(size=index) +
      scale_x_discrete(expand=c(0.02,0)) +
      scale_y_discrete(expand=c(0.02,0)) +
      theme(axis.text.x = element_text(angle=45, h=1),
            legend.position = "bottom",
            legend.box = "vertical")
   fwd =
     ggplot(tag_design, aes(x=tag_fwd, y=get(index))) +
       geom_boxplot(lwd=0.5) +
       theme_bw() +
       stat_summary(fun.y=median, geom="line", aes(group=1), color="cyan2", lwd=1) +
       labs(y=index, x=NULL) +
       scale_x_discrete(expand=c(0.02,0)) +
       scale_y_discrete(expand=c(0.02,0)) +
       coord_flip() +
       theme(axis.text.y = element_blank())

   rev =
     ggplot(tag_design, aes(x=tag_rev, y=get(index))) +
       geom_boxplot(lwd=0.5) +
       theme_bw() +
       stat_summary(fun.y=median, geom="line", aes(group=1), color="cyan2", lwd=1) +
       labs(y=index, x=NULL) +
       scale_x_discrete(expand=c(0.02,0)) +
       scale_y_discrete(expand=c(0.02,0)) +
       theme(axis.text.x = element_blank())


   out = insert_xaxis_grob(all, rev, grid::unit(.2, "null"), position = "top")
   out = insert_yaxis_grob(out, fwd, grid::unit(.2, "null"), position = "right")
   ggdraw(out)
  }
}

