#' Converts a metabarlist object in a phyloseq object
#'
#' Converts a \code{metabarlist} object in a \code{phyloseq} object for downstream diversity analyses.
#'
#' @param metabarlist   a \code{metabarlist} object
#' @param tax_colnames  a character vector with the column names to keep from the \code{motus} table of the \code{metabarlist} object
#'
#' @name metabarlist_to_phyloseq
#'
#' @return A \code{phyloseq} object with three components: \code{otu_table}, \code{sample_data} and \code{taxonomy_table}
#'
#' @details
#'
#' The function \code{metabarlist_to_phyloseq} converts a \code{metabarlist} object into a \code{phyloseq} object. A \code{phyloseq} object can contain several component data, such as \code{otu_table}, \code{sample_data} and \code{taxonomy_table}. The \code{otu_table} is equivalent to the \code{reads} matrix, the \code{sample_data} is equivalent to the \code{samples} data frame and the \code{taxonomy_table} corresponds to the taxonomic informations that may be contained in the \code{motus} data frame.
#' The sample names in the \code{reads} matrix (i.e. rownames) should be the the same as in the \code{samples} data frame.
#'
#' @seealso the \code{phyloseq} package
#'
#' @examples
#'
#' data(soil_euk) # load soil_euk data
#' soil_euk <- subset_metabarlist(soil_euk, "pcrs", metabarlist$pcrs$type == "sample") # select only samples, i.e. delete controls
#' soil_euk_agg <- aggregate_pcrs(soil_euk) # aggregate pcrs at sample level
#' tax_colnames <- c("phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name") # create vector with column names to keep from the motus table
#' soil_euk_physeq <- metabarlist_to_phyloseq(soil_euk, tax_colnames)
#' soil_euk_physeq
#'
#' @author Anne-Sophie Benoiston
#'
#' @importFrom phyloseq otu_table sample_data tax_table phyloseq
#'
#' @export metabarlist_to_phyloseq

metabarlist_to_phyloseq <- function(metabarlist, tax_colnames){
  if(suppressWarnings(check_metabarlist(metabarlist))){
    otu <- otu_table(metabarlist$reads, taxa_are_rows = FALSE)
    sample <- sample_data(metabarlist$samples)
    tax <- tax_table(as.matrix(metabarlist$motus[,tax_colnames]))
    return(phyloseq(otu, sample, tax))
  }
}

