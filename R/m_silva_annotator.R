#' Integrating SILVANgs pipeline taxonomic annotations
#'
#' Importing and formatting SILVANgs pipeline taxonomic annotation for a \code{\link{TODEFINE}} object.
#'
#'
#' @param x            a motus dataframe from a \code{\link{TODEFINE}} object
#' @param silva.path   path to a table from the SILVANgs pipeline, typically zipfile>ssu>exports>xxx---ssu---otus.csv
#' @param clust.path   path to a file from the SilvaNgs pipeline indicating otu cluster membership zipfile>ssu>stats>sequence_cluster_map>data>xxx---ssu---sequence_cluster_map---tmptaxo.clstr
#' @name silva_annotator
#'
#' @return a motus dataframe for a \code{\link{TODEFINE}} object including the taxonomic assignemnts from silva
#'
#' @details
#'
#' To amend
#'
#' @examples
#'
#' @author Lucie Zinger
#' @importFrom seqinr read.fasta
#' @export silva_annotator
#'

silva_annotator = function(motus, silva.path, clust.path) {
  silva = read.csv(silva.path, h=T, sep="\t", skip = 1, row.names = NULL)
  colnames(silva) = c(colnames(silva)[-c(1,ncol(silva))], "classif_ncbi", "classif_silva")

  #taxo formating
  tmp = sapply(strsplit(as.vector(silva$classif_silva), "\\|"), function(x) x[4])

  if(dim(table(is.na(tmp)))>1) {
    tmp[is.na(tmp)] = "No blast hit" #format according to qiime
  }

  taxorank = c("superkingdom_silva", "kingdom_silva", "phylum_silva", "class_silva",
               "order_silva", "family_silva", "genus_silva")

  tmp1 = do.call("rbind", lapply(strsplit(tmp, ";"), function(x) {
    out = rep(NA, length(taxorank))
    names(out) = taxorank
    if(x[1]=="Eukaryota") {
      out["superkingdom_silva"] = "Eukaryota"
      if(x[2]=="Opisthokonta") {
        if(length(grep("Aphelidea", x))!=0) {
          out["class_silva"] = ifelse(length(x)!=grep("Aphelidea", x), x[grep("Aphelidea", x)], NA)
          out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                        x[length(x)]=="uncultured", NA, x[length(x)])
        }
        if(length(grep("Holozoa", x))!=0 & length(grep("Metazoa", x))==0) {
          out["order_silva"] = ifelse(length(x)!=grep("Holozoa", x), x[grep("Holozoa", x)+1], NA)
          out["family_silva"] = ifelse(length(grep("dae$", x))!=0, x[grep("dae$", x)], NA)
          out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                        x[length(x)]=="uncultured", NA, x[length(x)])
        }
        if(length(grep("Metazoa", x))!=0) {
          out["kingdom_silva"] = "Metazoa"
          out["phylum_silva"] = ifelse(length(grep("Bilateria", x))!=0, x[grep("Bilateria", x)+1],
                                       ifelse(length(grep("Porifera", x))!=0, "Porifera",
                                              ifelse(length(grep("Cnidaria", x))!=0, "Cnidaria", NA)))
          a = paste("^", out["phylum_silva"], "$", sep="")
          if(length(grep("Porifera|Cnidaria", x))!=0) {
            out["class_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+1], NA)
            out["order_silva"] = ifelse(length(x)>(grep(a, x)+1), x[length(x)], NA)
            #might be class/subclass or families
          }
          if(length(grep("Nematoda|Mollusca|Nemertea|Tardigrada|Rotifera|
                         Entoprocta|Bryozoa|Platyhelminthes|Annelida", x))!=0) {
            out["class_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+1], NA)
            out["order_silva"] = ifelse(length(x)>(grep(a, x)+1), x[length(x)], NA)
            #might be class/subclass or families
        }
          if(length(grep("Arthropoda|Brachiopoda", x))!=0) {
            out["class_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+2], NA)
            out["order_silva"] = ifelse(length(x)>(grep(a, x)+2), x[length(x)], NA)
            #might be class/subclass for certain groups (e.g. collembola)
          }
          if(length(grep("Gastrotricha|Echinodermata", x))!=0) {
            out["order_silva"] = ifelse(length(x)>(grep(a, x)+1), x[length(x)], NA)
            #might be families
          }
          if(length(grep("Chordata", x))!=0) {
            out["class_silva"] = ifelse(length(x)>(grep(a, x)+1), x[length(x)], NA)
          }}
        if(x[3]=="Nucletmycea") {
          if(length(grep("Fungi", x))!=0) {
            out["kingdom_silva"] = "Fungi"
            out["phylum_silva"] = ifelse(length(grep("mycota$", x))!=0, x[grep("mycota$", x)], NA)
            out["class_silva"] = ifelse(length(grep("mycetes$", x))!=0, x[grep("mycetes$", x)], NA)
            out["order_silva"] = ifelse(length(grep("ales$", x))!=0, x[grep("ales$", x)], NA)
            out["family_silva"] = ifelse(length(grep("ceae$", x))!=0, x[grep("ceae$", x)], NA)
            out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                          x[length(x)]=="uncultured", NA, x[length(x)])
          } else {
            out["family_silva"] = ifelse(length(grep("ae$", x))!=0, x[grep("ae$", x)], NA)
            out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                          x[length(x)]=="uncultured", NA, x[length(x)])
          }}
      } else if(length(grep("Archaeplastida", x)!=0)){
        if(length(grep("Chloroplastida", x)!=0)){
          out["kingdom_silva"] = "Viridiplantae"
          if(length(grep("Streptophyta", x))!=0){
            out["phylum_silva"] = ifelse(length(grep("Streptophyta", x))!=0,
                                         x[grep("Streptophyta", x)], NA)
            out["order_silva"] = ifelse(length(grep("ales$", x))!=0, x[grep("ales$", x)], NA)
            out["family_silva"] = ifelse(length(grep("ceae$", x))!=0, x[grep("ceae$", x)], NA)
            out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                          x[length(x)]=="uncultured", NA, x[length(x)])
          } else {
            out["phylum_silva"] = ifelse(length(grep("Chlorophyta", x))!=0,
                                         x[grep("Chlorophyta", x)], NA)
            out["class_silva"] = ifelse(length(grep("ceae$", x))!=0, x[grep("ceae$", x)], NA)
            out["order_silva"] = ifelse(length(grep("ales$", x))!=0, x[grep("ales$", x)], NA)
            out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                          x[length(x)]=="uncultured", NA, x[length(x)])
          }} else {
            out["class_silva"] = ifelse(length(grep("ceae$|phyta$", x))!=0,
                                        x[grep("ceae$|phyta$", x)], NA)
            out["order_silva"] = ifelse(length(grep("ales$", x))!=0, x[grep("ales$", x)], NA)
            out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                          x[length(x)]=="uncultured", NA, x[length(x)])
          }
      } else if(length(grep("^SAR$", x))!=0) {
        out["kingdom_silva"] = "SAR" #not kingdom but it helps
        out["phylum_silva"] = ifelse(length(x)!=grep("^SAR$", x),
                                     x[grep("^SAR$", x)+1], NA)
        if(out["phylum_silva"]=="Alveolata") {
          out["phylum_silva"] = ifelse(x[grep("Alveolata", x)+1]=="uncultured", NA,
                                       x[grep("Alveolata", x)+1])
        }
        if(is.na(out["phylum_silva"])==F){
          a = paste("^", out["phylum_silva"], "$", sep="")
          #actually no rank, but it helps for SARs
          out["class_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+1], NA)
          a = paste("^", out["class_silva"], "$", sep="")
          if(length(grep(a, x))==0) {
            out["order_silva"] = NA
          } else {
            out["order_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+1], NA)
          }
          out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                        x[length(x)]=="uncultured", NA, x[length(x)])
        }} else if(length(grep("Amoebozoa", x))!=0) {
          out["phylum_silva"] = "Amoebozoa"
          a = paste("^", out["phylum_silva"], "$", sep="")
          out["class_silva"] = ifelse(length(x)!=grep(a, x), x[grep(a, x)+1], NA)
          out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                        x[length(x)]=="uncultured", NA, x[length(x)])
          #mess
        } else {
          out["phylum_silva"] = x[2]
          out["genus_silva"] = ifelse(x[length(x)] %in% out | x[length(x)]=="Incertae Sedis" |
                                        x[length(x)]=="uncultured", NA, x[length(x)])
        }} else {
          out["superkingdom_silva"] = out["kingdom_silva"] = x[1]
          out[3:(2+length(x)-1)] = x[2:length(x)]
        }
    out
  }))

  tmp1 = as.data.frame(tmp1)
  if(length(grep("uncultured", tmp1$genus_silva))>1) {
    tmp1$genus_silva = as.vector(tmp1$genus_silva)
    idx = grep("uncultured", tmp1$genus_silva)
    tmp1$genus_silva[idx] = NA
  }

  tmp1 = cbind.data.frame(silva[,c("cluster.acc", "cluster.id", "similarity", "X..sequences")],
                          tmp1, lineage_silva=tmp)
  #propagate
  heads = which(tmp1$X..sequences>1)
  if(length(heads)>0) {

    clust = read.fasta(clust.path, as.string = T, forceDNAtolower = F)
    clust = clust[match(tmp1$cluster.acc, names(clust))]

    tmp2 = do.call("rbind", lapply(heads, function(x) {
      #existing seq
      h = as.vector(unname(tmp1$cluster.acc[x]))
      tax = tmp1[match(h, tmp1$cluster.acc),]

      gp = unname(gsub(" >", "", gsub("\\..*", "", unlist(strsplit(unlist(clust[x]), ","))[-1])))
      gp = gp[grep(h, gp, invert = T)]
      tax2 = vector("list", length = length(gp))
      for(i in 1:length(tax2)) {tax2[[i]]=tax}
      data.frame(cluster.acc = gp, tax[rep(1, length(gp)),-1])
    }))

    tmp1 = rbind.data.frame(tmp1, tmp2)
  }

  #merge with motus table
  idx = match(rownames(motus), as.vector(tmp1$cluster.acc))
  out = cbind.data.frame(motus, tmp1[idx,-c(1:2,4)])

  return(out)
}

