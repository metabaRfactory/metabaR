#' Soil eukaryote DNA metabarcoding data from French Guiana
#'
#' Data from a DNA metabarcoding experiment aiming at assessing the diversity of soil eukaryotes in two tropical forests of French Guiana.
#'
#' @details
#' Samples were collected in 2 sampling sites in contrasting habitats:
#' \itemize{
#'   \item{}{The Mana site located in a white sand forest, characterized highly oligotrophic soils and particular tree species adaptated to the local harsh conditions.}
#'   \item{}{The Petit Plateau site is located in the pristine rainforest (Nouragues natural reserve) characterized by soils richer in clay and organic matter.}
#'  }
#' At each sites, 16 samples (mesh grid size of 20 m) were collected in 1 plot of 1 ha.
#' For each sampling points, two types of compartments:
#' \itemize{
#'   \item{}{Soil: composite sample of 5 soil cores}
#'   \item{}{Litter: sampling of XXXX m3 of litter}
#'   \item{}{=>A total of 64 DNA extracts + 4 extraction blank controls (1 per site and compartment)}
#' }
#' For each DNA extracts, a short region of the 18S rRNA (Taberlet et al. 2018) was
#' amplfied by PCR in quadruplicate, following the protocol described in Zinger et al. 2019.
#' The resulting amplicons were pooled and sequenced on a Illumina HiSeq plateform, using the
#' the paired-end technology.
#'
#' The total experiment hence includes the following amplicon types:
#' \itemize{
#'   \item{}{4 replicates / samples (n = 256)}
#'   \item{}{4 replicates for the 4 extraction blank controls (n = 16)}
#'   \item{}{4 replicates for 8 PCR blank controls (n = 32)}
#'   \item{}{4 replicates for 12 sequencing blank controls (n = 48)}
#'   \item{}{4 replicates for 8 positive controls (plant DNA from 16 species; n = 32)}
#'   \item{}{=> 384 amplicons in total}
#' }
#'
#' The retrieved data were then processed using the OBITools (Boyer et al. 2016) and SUMACLUST (Mercier et al. 2013) packages. Briefly, paired-end reads were assembled, assigned to their respective samples/marker and dereplicated. Low-quality sequences (containing Ns, shorter than 50 bp or singletons) were excluded; the remaining ones were clustered into operational taxonomic units (OTUs) using SUMACLUST at a sequence similarity threshold of 0.97. The representative sequence of each OTU (most abundant one) was assigned to taxonomic clade using a databased built from the EMBL (release 136) with the ecoPCR program (Ficetola et al., 2010).
#'
#'
#' The data `soil_euk` is a \code{metabarlist} containing four tables
#'
#' \itemize{
#' \item{}{`reads`: a numeric matrix containing the MOTUs abundances (expressed as a number of reads) for each pcr (i.e. technical replicates of both biological samples and positive and negative controls)}
#' \item{}{`motus`: a dataframe containing the MOTUs characteristics (e.g. taxonomy, sequence) for each MOTUs)}
#' \item{}{`pcrs`: a dataframe containing information on each pcr (e.g. control type, pcr wells, etc.)}
#' \item{}{`samples`: a dataframe containing information on each environmental sample (e.g. habitat type, etc.)}
#' }
#'
#' @docType data
#'
#' @usage data(soil_euk)
#'
#' @format An object of class \code{metabarlist}; see \code{\link[metabaRffe]{metabarlist}}.
#'
#' @keywords dataset
#'
#' @references Boyer, F., Mercier, C., Bonin, A., Le Bras, Y., Taberlet, P., & Coissac, E. (2016). obitools: a unix‐inspired software package for DNA metabarcoding. Molecular ecology resources, 16(1), 176-182.
#' @references Ficetola, G. F., Coissac, E., Zundel, S., Riaz, T., Shehzad, W., Bessière, J., ... & Pompanon, F. (2010). An in silico approach for the evaluation of DNA barcodes. BMC genomics, 11(1), 434.
#' @references Mercier, C., Boyer, F., Bonin, A., & Coissac, E. (2013, November). SUMATRA and SUMACLUST: fast and exact comparison and clustering of sequences. In Programs and Abstracts of the SeqBio 2013 workshop. Abstract (pp. 27-29).
#' @references Taberlet, P., Bonin, A., Zinger, L., & Coissac, E. (2018). Environmental DNA: For Biodiversity Research and Monitoring. Oxford University Press.
#' @references Zinger, L., Taberlet, P., Schimann, H., Bonin, A., Boyer, F., De Barba, M., ... & Chave, J. (2019). Body size determines soil community assembly in a tropical forest. Molecular ecology, 28(3), 528-543.
#'
#'
#' @source \href{LINK WEB TO DEFINE}{DNA metabarcoding data archive}
#'
#' @examples
#' data(soil_euk)
#' lapply(soil_euk, dim)
"soil_euk"
