## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE-------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("metabaRfactory/metabaRffe")

## ----loadpackage---------------------------------------------------------
library(metabaRffe) # modify the name once we'll all agree on that

## ----help----------------------------------------------------------------
?soil_euk

## ----soil_euk_data-------------------------------------------------------
data(soil_euk) 
summary_metabarlist(soil_euk)

## ----namesex1------------------------------------------------------------
colnames(soil_euk$pcrs)

## ----namesex2------------------------------------------------------------
colnames(soil_euk$samples)

## ----import, eval=F------------------------------------------------------
#  soil_euk <- tabfiles_to_metabarlist(file_reads = "litiere_euk_reads.txt",
#                                      file_motus = "litiere_euk_motus.txt",
#                                      file_pcrs = "litiere_euk_pcrs.txt",
#                                      file_samples = "litiere_euk_samples.txt")

## ----diag1_readsmotus----------------------------------------------------
#compute the number of reads per pcr
soil_euk$pcrs$nb_reads <- rowSums(soil_euk$reads)
#compute the number of motus per pcr
soil_euk$pcrs$nb_motus <- rowSums(soil_euk$reads>0)

## ----diag1_boxplotreadsmotus, warning=F, message=F, fig.width=7, fig.height=4----
#load requested package for plotting
library(ggplot2)
library(reshape2)

#create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- melt(soil_euk$pcrs[,c("control_type", "nb_reads", "nb_motus")])

ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

## ----diag2_readsMOTUs, message=F, warning=F, fig.width=5, fig.height=4----
#Using the nb_reads and nb_motus defined previously in the soil_euk$pcrs table
ggplot(soil_euk$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")

## ----ggpcrplate1, warning=F, message=F, fig.width=7, fig.height=5--------
ggpcrplate(soil_euk, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

## ----ggpcrplate2, warning=F, message=F, fig.width=7, fig.height=5, eval=F----
#  ggpcrplate(soil_euk)

## ----ggpcrtag, warning=F, message=F, fig.width=7, fig.height=7-----------
#Here the list of all tag/indices used in the experiment is available in the column "tag_rev" of the soil_euk$pcrs table
tag.list <- as.character(unique(soil_euk$pcrs$tag_rev))
ggpcrtag(soil_euk, legend_title = "# of reads per PCR", FUN = function(m) {rowSums(m$reads)},
                     taglist = tag.list) 

## ----subset, message=F, warning=F----------------------------------------
#get the samples names from the H20 plot
h20_id <- rownames(soil_euk$pcrs)[grep("H20-[A-B]", rownames(soil_euk$pcrs))]

#subset the data
soil_euk_h20 <- subset_metabarlist(soil_euk, table = "pcrs", indices = h20_id)

#check results
summary_metabarlist(soil_euk_h20)

## ----hillraref, message=F, warning=F-------------------------------------
soil_euk_h20.raref = hill_rarefaction(soil_euk_h20, nboot = 20, nsteps = 10)
head(soil_euk_h20.raref$hill_table)

## ----gghill, message=F, warning=F, fig.width=7, fig.height=2.5-----------
gghill_rarefaction(soil_euk_h20.raref) 

## ----gghill2, message=F, warning=F, fig.width=7, fig.height=3------------
#define a vector containing the Material info for each pcrs 
material <- soil_euk_h20$samples$Material[match(soil_euk_h20$pcrs$sample_id,
                                               rownames(soil_euk_h20$samples))]

#use of gghill_rarefaction requires a vector with named pcrs
material <- setNames(material,rownames(soil_euk_h20$pcrs))

#plot
p <- gghill_rarefaction(soil_euk_h20.raref, group=material)
p + scale_fill_manual(values = c("goldenrod4", "brown4", "grey")) +
    scale_color_manual(values = c("goldenrod4", "brown4", "grey")) +
    labs(color="Material type")

## ----contaslayer, message=F, warning=F, fig.width=6, fig.height=3--------
#Define a vector containing the extraction negative control names
ext.controls <- rownames(soil_euk$pcrs)[which(soil_euk$pcrs$control_type=="extraction")]

#now detect contaminants within your metabarlist
conta.ext <- contaslayer(soil_euk, controls = ext.controls)

#tag these sequences in the metabarlist object
soil_euk$motus$flagged_motus <- NA
soil_euk$motus[conta.ext, "flagged_motus"]  = "extraction contaminant"

## ----contaslayer2, message=F, warning=F, echo=F, fig.width=6, fig.height=3----
library(kableExtra)
dt <- soil_euk$motus[!is.na(soil_euk$motus$flagged_motus), c("count", "best_identity.order_filtered_embl_r136_noenv_EUK", "path", "sequence")]
dt$best_identity.order_filtered_embl_r136_noenv_EUK <- round(dt$best_identity.order_filtered_embl_r136_noenv_EUK )
colnames(dt) <- c("total # reads", "similarity to ref", "full taxonomic path", "sequence")

rownames(dt) <- NULL

kable(dt[order(dt[,1], decreasing = TRUE)[1:10],]) %>%
  kable_styling(bootstrap_options= c("striped", "hover", "condensed"), 
                font_size = 8, full_width = F)

## ----contaslayer3, message=F, warning=F, fig.width=7, fig.height=6-------
#Identify the most common contaminant
max.conta = conta.ext[which.max(soil_euk$motus[conta.ext, "count"])]

#Distribution of the most abundant contaminant in the PCR plate design in terms of relative abundance
ggpcrplate(soil_euk, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

## ----contaslayer4, message=F, warning=F, fig.width=4, fig.height=3-------
#compute relative abundance of all pcr contaminants together 
a <- data.frame(conta.relab = rowSums(soil_euk$reads[,!is.na(soil_euk$motus$flagged_motus)]) / 
                                    rowSums(soil_euk$reads))
#add information on control types
a$control_type <- soil_euk$pcrs$control_type[match(rownames(a), rownames(soil_euk$pcrs))]

ggplot(a, aes(x=control_type, y=conta.relab, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") + 
  theme_bw() + 
  scale_y_log10()

## ----contaslayer5, message=F, warning=F,  fig.width=4, fig.height=3------
#set a minimum value of total contaminant relative abundance
thresh <- 1e-1
#and flagg
soil_euk$pcrs$flagged_pcrs = ifelse(a$conta.relab[match(rownames(soil_euk$pcrs), rownames(a))]>thresh, 
                             "too much conta pcr", NA)

