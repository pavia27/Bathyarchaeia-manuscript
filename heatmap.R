  #full
  rm(list = ls(all.names = TRUE))
  require(gplots)
  library(RColorBrewer)
  aai <-read.csv("aai_scores.upper_triangle.csv", row.names = 1, header = TRUE)
  a<-aai
  aai<-aai[,-c(1,2)]
  aai[is.na(aai)] <- 100
  aai_plot <- as.matrix(aai)
  breaks <- c(0, 45, 65, 95, 100)
  col = c("#2c7bb6","#abd9e9","#fdae61","#d7191c")
  hclustmethod <- function(x) hclust(x, method="average")
  heatmap.2(aai_plot, keysize=1.5, key.xlab="AAI", key.title = "", col=col, trace="none", margins = c(6, 6), cexRow = 0.1, cexCol = 0.1,dendrogram='both',density.info="none",Rowv = T, Colv = T,sepcolor="white",sepwidth=c(0,0),breaks = breaks,offsetRow = 0.000001,offsetCol = 0.000001,ColSideColors = a$BC,hclustfun=hclustmethod)
  #legend("left",legend = unique(a$TO),col = unique(a$BC),lty= 1,lwd = 5,cex=.7)
  library(pheatmap)
  #sub
  rm(list = ls(all.names = TRUE))
  require(gplots)
  library(RColorBrewer)
  aai <-read.csv("aai_scores.upper_triangle_sub.csv",header = TRUE)
  a<-aai
  a$BC<-a$Bathyarchaeota_Clade
  a$Bathyarchaeota_Clade[grepl("BC15",a$Bathyarchaeota_Clade)]<-"#F8961E"
  a$Bathyarchaeota_Clade[grepl("BC36",a$Bathyarchaeota_Clade)]<-"#FFD043"
  a$Bathyarchaeota_Clade[grepl("BC38",a$Bathyarchaeota_Clade)]<-"#7FC96B"
  a$Bathyarchaeota_Clade[grepl("BC39",a$Bathyarchaeota_Clade)]<-"#43AA8B"
  a$Bathyarchaeota_Clade[grepl("BC40",a$Bathyarchaeota_Clade)]<-"#277DA1"
  a$Bathyarchaeota_Clade[grepl("BC41",a$Bathyarchaeota_Clade)]<-"#3B498E"
  a$Bathyarchaeota_Clade[grepl("BC42",a$Bathyarchaeota_Clade)]<-"#66418A"
  a$Bathyarchaeota_Clade[grepl("BC1",a$Bathyarchaeota_Clade)]<-"#E63232"
  a$Bathyarchaeota_Clade[grepl("BC3",a$Bathyarchaeota_Clade)]<-"#F3722C"
  rownames(aai)<-aai[,1]
  aai<-aai[,-c(1:3)]
  aai[is.na(aai)] <- 100
  aai_plot <- as.matrix(aai)
  breaks <- c(0, 45, 63, 95, 100)
  col = c("#2c7bb6","#abd9e9","#fdae61","#d7191c")
  heatmap.2(aai_plot, keysize=1.5, key.xlab="AAI", key.title = "", col=col, trace="none", margins = c(6, 6), cexRow = 0.1, cexCol = 0.1,dendrogram='both',density.info="none",Rowv = T, Colv = T,sepcolor="white",sepwidth=c(0,0),breaks = breaks,offsetRow = 0.000001,offsetCol = 0.000001,ColSideColors = a$Bathyarchaeota_Clade)
  #legend("left",legend = unique(a$BC),col = unique(a$Bathyarchaeota_Clade),lty= 1,lwd = 5,cex=.7)  

