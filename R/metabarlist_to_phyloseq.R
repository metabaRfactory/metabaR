#' Converts a metabarlist object in a phyloseq object
#'
#' Converts a \code{metabarlist} object in a \code{phyloseq} object for downstream diversity analyses.
#'
#' @param metabarlist       a \code{metabarlist} object aggregated at sample level
#' @param tax_colnames      a character vector with the column names to keep from the \code{motus} table of the \code{metabarlist} object
#' @param include_controls  logical: if \code{TRUE}, controls are added in the \code{phyloseq} \code{sample_data} object, as well as the columns of the \code{pcrs} table. MOTUs abundance in controls should be present in the \code{motus} table The default (\code{FALSE}) is not to include controls. In that case, the columns of the \code{pcrs} table are added to the \code{sample_data} object and the \code{motus} table should not include abundances in controls.
#'
#' @name metabarlist_to_phyloseq
#'
#' @return A \code{phyloseq} object with three components: \code{otu_table}, \code{sample_data} and \code{taxonomy_table}
#'
#' @details
#'
#' The function \code{metabarlist_to_phyloseq} converts a \code{metabarlist} object into a \code{phyloseq} object. A \code{phyloseq} object can contain several component data, such as \code{otu_table}, \code{sample_data} and \code{taxonomy_table}. The \code{otu_table} is equivalent to the \code{reads} matrix, the \code{sample_data} is equivalent to the \code{samples} data frame and the \code{taxonomy_table} corresponds to the taxonomic informations that may be contained in the \code{motus} data frame. Information from the \code{pcrs} table is included in the \code{sample_data} component.
#'
#' @references \url{https://joey711.github.io/phyloseq/}
#'
#' @seealso \code{\link{subset_metabarlist}}, \code{\link{aggregate_pcrs}}, the \code{phyloseq} package
#'
#' @examples
#'
#' data(soil_euk) # load soil_euk data
#' soil_euk_sub <- subset_metabarlist(soil_euk, "pcrs", soil_euk$pcrs$type == "sample")
#' soil_euk_agg <- aggregate_pcrs(soil_euk_sub)
#' tax_colnames <- c("phylum_name", "class_name", "order_name", "family_name", "genus_name",
#'  "species_name")
#' soil_euk_physeq <- metabarlist_to_phyloseq(soil_euk_agg, tax_colnames)
#' soil_euk_physeq
#'
#' soil_euk_agg <- aggregate_pcrs(soil_euk)
#' soil_euk_physeq <- metabarlist_to_phyloseq(soil_euk_agg, tax_colnames,
#'  include_controls = TRUE)
#' soil_euk_physeq
#'
#' @author Anne-Sophie Benoiston
#'
#' @importFrom phyloseq otu_table sample_data tax_table phyloseq
#'
#' @export metabarlist_to_phyloseq

metabarlist_to_phyloseq <- function(metabarlist, tax_colnames, include_controls = FALSE){
  if(suppressWarnings(check_metabarlist(metabarlist))){
    if(!include_controls){
      if(!setequal(rownames(metabarlist$pcrs), rownames(metabarlist$samples))){
        stop("The pcrs and samples tables should have the same row names and not include controls.")
      }
    }
    else if(setequal(rownames(metabarlist$pcrs), rownames(metabarlist$samples))){
      warning("The pcrs and samples tables have the same row names, so it seems that you are not including controls in the phyloseq object.")
    }
  otu <- otu_table(metabarlist$reads, taxa_are_rows = FALSE)
  sample <- merge(metabarlist$samples, metabarlist$pcrs, by="row.names", all=TRUE)
  rownames(sample) <- sample$Row.names
  sample$Row.names <- NULL
  sample <- sample_data(sample)
  tax <- tax_table(as.matrix(metabarlist$motus[,tax_colnames]))
  return(phyloseq(otu, sample, tax))
  }
}

