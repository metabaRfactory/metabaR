samples = read.table("Litiere_sample_list.txt", h=T, row.names=1, sep="\t")
tmp = read.csv("litiere_euk_cl97_agg_filt_tax.tab", h=T, sep="\t")

reads = t(tmp[,grep("sample\\.", colnames(tmp))])
rownames(reads) = gsub("sample\\.", "", rownames(reads))
colnames(reads) = tmp$id

motus = tmp[,grep("sample\\.", colnames(tmp), invert = T)]
rownames(motus) = motus$id
motus = motus[,c("count", "GC_content", "seq_length",
                 "best_identity.order_filtered_embl_r136_noenv_EUK",
                 "taxid_by_db.order_filtered_embl_r136_noenv_EUK",
                 "phylum_name", "class_name", "order_name", "family_name",
                 "genus_name", "species_name",
                 "rank", "scientific_name",
                 "sequence")]

samp = samples[match(rownames(reads), rownames(samples)),]


soil_euk = list(reads = reads,
                motus = motus,
                samples = samp)
devtools::use_data(soil_euk)
