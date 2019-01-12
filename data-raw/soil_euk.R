#obitab output
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
                 "rank", "scientific_name", "path",
                 "sequence")]

#samples characteristics
samples = read.table("Litiere_sample_list.txt", h=T, row.names=1, sep="\t")

#plate locations in ngsfilter
ngsfilt = read.table("ngsfilter_GWM-768.new.txt", h=F, sep="\t")
colnames(ngsfilt) = c("project", "amplicon", "tagcombo", "primerF", "primerR", "info")

idx = match(rownames(samples), gsub("-", "_", as.vector(ngsfilt$amplicon)))
samples$project = ngsfilt[idx,"project"]
samples$tag_fwd = sapply(strsplit(as.vector(ngsfilt[idx,"tagcombo"]), "\\:"), "[[", 1)
samples$tag_rev = sapply(strsplit(as.vector(ngsfilt[idx,"tagcombo"]), "\\:"), "[[", 2)
samples$primer_fwd = ngsfilt[idx,"primerF"]
samples$primer_rev = ngsfilt[idx,"primerR"]
samples$plate_no = sapply(strsplit(gsub("F @ position=|;", "", ngsfilt[idx,"info"]), "_"), "[[", 1)
samples$plate_col = gsub("[A-Z]", "", sapply(strsplit(gsub("F @ position=|;", "", ngsfilt[idx,"info"]), "_"), "[[", 2))
samples$plate_row = gsub("[0-9]", "", sapply(strsplit(gsub("F @ position=|;", "", ngsfilt[idx,"info"]), "_"), "[[", 2))

samp = samples[match(rownames(reads), rownames(samples)),]


soil_euk = list(reads = reads,
                motus = motus,
                pcrs = samp)
devtools::use_data(soil_euk)
