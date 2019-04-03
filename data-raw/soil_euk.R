# import obitab output
tmp = read.csv("data-raw/litiere_euk_cl97_agg_filt_tax.tab", h=T, sep="\t", check.names = F)

# #Get reads table
# reads = t(tmp[,grep("sample\\:", colnames(tmp))])
# rownames(reads) = gsub("sample\\:", "", rownames(reads))
# colnames(reads) = tmp$id
# reads = cbind(rownames(reads), reads)
# colnames(reads)[1] = "pcr_id"
# write.table(reads, file = "data-raw/litiere_euk_reads.txt",
#              sep="\t", quote=F, row.names = F)

reads = as.matrix(read.csv2("data-raw/litiere_euk_reads.txt", row.names=1, h=T, sep="\t",
                             check.names = F))


# Get motus table
# motus = tmp[,grep("sample\\.", colnames(tmp), invert = T)]
# rownames(motus) = motus$id
# motus = motus[,c("id", "count", "GC_content", "seq_length",
#                  "best_identity.order_filtered_embl_r136_noenv_EUK",
#                  "taxid_by_db.order_filtered_embl_r136_noenv_EUK",
#                  "phylum_name", "class_name", "order_name", "family_name",
#                  "genus_name", "species_name",
#                  "rank", "scientific_name", "path",
#                  "sequence")]
#write.table(motus, file = "data-raw/litiere_euk_motus.txt",
#             sep="\t", quote=F, row.names = F)

motus = read.table("data-raw/litiere_euk_motus.txt", row.names=1, h=T, sep="\t")


all(colnames(reads) == rownames(motus))


# Get pcr table
pcrs = read.table("data-raw/PCR_information.tsv", sep="\t", h=T)
rownames(pcrs) = as.vector(pcrs$pcr_id)

all(rownames(reads) %in% rownames(pcrs))

pcrs = pcrs[,-match("pcr_id", colnames(pcrs))]
pcrs = data.frame(pcr_id = rownames(pcrs), pcrs)
write.table(pcrs, file = "data-raw/litiere_euk_pcrs.txt",
                         sep="\t", quote=F, row.names = F)

reads = reads[match(rownames(pcrs), rownames(reads)),]
reads[is.na(reads)] = 0
rownames(reads) = rownames(pcrs)

# Get samples table and make it as an input table for the import data function
# tmp = read.table("data-raw/Litiere_sample_list.txt", h=T, row.names=1, sep="\t")
# #format to obtain only samples description
# samples = unique(tmp[,1:8])
# #remove controls
# samples = samples[which(is.na(samples$Control_type)),-c(9:10)]
# #rename samples table
# rownames(samples) = substr(rownames(samples), 1, nchar(rownames(samples))-3)
# samples = data.frame(sample_id = rownames(samples), samples)
# #export table
# write.table(samples, file = "data-raw/litiere_euk_samples.txt",
#             sep="\t", quote=F, row.names = F)
#

samples = read.table("data-raw/litiere_euk_samples.txt", h=T, row.names = 1, sep="\t")

all(rownames(samples) %in% pcrs$sample_id)

dim(motus)
dim(reads)
dim(pcrs)
dim(samples)
nlevels(droplevels(pcrs$sample_id[pcrs$type=="sample"]))
all(rownames(reads) == rownames(pcrs))
all(colnames(reads) == rownames(motus))


soil_euk = list(reads = reads,
                motus = motus,
                pcrs = pcrs,
                samples = samples)

devtools::use_data(soil_euk)
