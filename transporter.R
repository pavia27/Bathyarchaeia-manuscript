  rm(list = ls(all.names = TRUE))
  library(egg)
  library(tidyverse)
  library(reshape2)
  library(pheatmap)
  a<-read.csv("Transporter.csv")
  b<- data.frame(Clade = colnames(a)[7:99], Genome_Name = unlist(a[1, 7:99])) %>%
    mutate(Clade = stringr::str_replace(Clade, "\\..*", ""))
  b<- b %>% `rownames<-`(.[,2]) %>% select(-Genome_Name)
  c<-read.csv("Transporter.csv", skip = 1)
  d<-c[,6:99] %>% `rownames<-`(.[,1]) %>% select(-gene)
  d[is.na(d)]<-0
  e<-c[,c(6,4)] %>% `rownames<-`(.[,1]) %>% select(-gene)
  i<-list(
    Clade = c(BC42="#66418A",BC41="#3B498E",BC40="#277DA1",BC39="#43AA8B",BC38="#7FC96B",BC36="#FFD043",BC15="#F8961E",BC3="#F3722C",BC1="#E63232"),
    pathway=c("Sulfur_phosphate"="#dadaeb","Nitrogen"="#bcbddc","Iron"="#9e9ac8","Metal"="#756bb1","amino_acid"="#54278f"))
  pheatmap(d,cluster_cols=T,cluster_rows=T,show_rownames=T,color=c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252","#252525","#000000"),clustering_method = "ward.D",fontsize_col=4,fontsize_row=4,border_color=NA,annotation_legend=T,cellwidth=4,cellheight=4,annotation_col=b,annotation_row = e,annotation_colors =i)
