install.packages("Seurat")
library(dplyr)
library(patchwork)
library(Seurat)
library(tidyverse)
library(gridExtra)
#load seurat data 
cell_data <- readRDS("/projectnb/bf528/users/frazzled/project_4/programmer/pbmc_seurat.rds")
cell_data2<- readRDS("/projectnb/bf528/users/frazzled/project_4/programmer/pbmc_seurat2.rds")
#find markers 
#previous format: cell_markers <- FindAllMarkers(cell_data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) #%>% 
#group_by(cluster) %>% 
#top_n(n = 5, wt = avg_log2FC) This filters the cluster 

cell_markers <- FindAllMarkers(cell_data, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.1) 
cell_markers2 <- FindAllMarkers(cell_data2, only.pos = TRUE,min.pct = 0.1,logfc.threshold = 0.1) 
#min = 0.5 logfc.threshold =0.1

#take the gene col and export as csv to upload into ensembl to obtain corresponding gene names and give back to programmer to load into
#seurat object 
#e_genes <- cell_markers["gene"]
#write_csv(e_genes,"Ensembl_GeneID")

#create csv files of all diff genes and diff genes expressed with adj p value <0.05 on find all markers 
marker_genes2<- write_csv(cell_markers2,"Marker_Genes")
padj_marker_genes2 <- cell_markers[cell_markers2$p_val_adj<0.05,]
padj_genes2<- write_csv(padj_marker_genes2,"Padj_Marker_Genes")

#run UMAP 
UMAP_markers <-RunUMAP(cell_data2, dims = 1:10) #label by cell type 

UMAP <- DimPlot(UMAP_markers, reduction= "umap", label = TRUE)
UMAP
#save UMAP
ggsave("UMAP_Plt.png",units = "cm", width = 20, height = 10, scale = 3)
#create tSNE map 
tSNE <- RunTSNE(cell_data2,dims = 1:10)
tSNE_map<-DimPlot(tSNE,label = TRUE)
tSNE_map
#save tSNE map
ggsave("tSNE_map.png",units="cm",width = 20, height = 10, scale = 3)
#create vector of cluster names 
cluster_names <- c("Alpha","Beta","Delta","Gamma","Epsilon","Acinar","Ductal","Quiescent steliate","Activated steilate","Endothelial",
                   "Macrophage","Mast","Cytotoxic T","Schwann")
names(cluster_names) <- levels(cell_data2)
cell_data2 <-RenameIdents(cell_data2,cluster_names)


#pulling top 5 from each cluster  
top_5 <- cell_markers2 %>% 
  group_by(cluster) %>%
  slice_max(n=5,order_by = avg_log2FC)

#pulling top 3 from each cluster  
top_3 <- cell_markers2 %>% 
  group_by(cluster) %>%
  slice_max(n=3,order_by = avg_log2FC)

top_10 <- cell_markers2 %>% 
  group_by(cluster) %>%
  slice_max(n=10,order_by = avg_log2FC)

#saving top 3 data as png 
png("top3.png", height = 25*nrow(top_3), width =85*ncol(top_3))
grid.table(top_3)
dev.off()
#saving top 5 data as png
png("top5.png", height = 25*nrow(top_5), width =85*ncol(top_5))
grid.table(top_5)
dev.off()
#saving top 10 data as png ** too big dont include 
png("top10.png", height = 25*nrow(top_10), width =85*ncol(top_10))
grid.table(top_10)
dev.off()

#create heatmaps 
Top5_HeatMap<- DoHeatmap(cell_data2, features = top_5$gene) + NoLegend()
Top5_HeatMap
Top3_HeatMap<-DoHeatmap(cell_data2,features = top_3$gene) + NoLegend()
Top3_HeatMap
ggsave("Top 5 HeatMap.png", units = "cm",width = 25, height = 12.5,scale = 2)
DoHeatmap(cell_data2, features = top_3$gene) + NoLegend() 

genes_VlnPlt<- VlnPlot(cells, features = c("GCG","INS","SST","PPY","GHRL",
                                           "CPA1","KRT19","RGS5","PDGFRA","VWF","SDS","TPSAB1","TRAC"), log = TRUE)
genes_VlnPlt
VlnPlot(cell_data2,features = c("GCG","INS","SST","PPY","GHRL","CPA1","KRT19",
                                "RGS5","PDGFRA","VWF","SDS","TPSAB1","TRAC"), log = TRUE)
Vln_plt_1<-VlnPlot(cell_data2,features = c("GCG","INS","SST","PPY","GHRL","CPA1","KRT19"),log = TRUE)
Vln_plt_1
ggsave("Vln_plts.png", units = "cm",width = 20, height = 10,scale = 2)
Vln_plt_2<-VlnPlot(cell_data2,features = c("RGS5","PDGFRA","VWF","SDS","TPSAB1","TRAC"), log = TRUE)
Vln_plt_2
#cowplot::plot_grid(Vln_plt_1,Vln_plt_2)
ggsave("Vln_plts_2.png", units = "cm",width = 20, height = 10,scale = 2)

Gene_UMAPS_1 <- FeaturePlot(UMAP_markers,features = c("GCG","INS","SST","PPY")) 
Gene_UMAPS_1
Gene_UMAPS_2<-FeaturePlot(UMAP_markers,features = c("GHRL","CPA1","KRT19","RGS5"))
Gene_UMAPS_2
Gene_UMAPS_3<-FeaturePlot(UMAP_markers,features = c("PDGFRA","VWF","SDS","TPSAB1"))
Gene_UMAPS_3
Gene_UMAPS_4<-FeaturePlot(UMAP_markers,features = c("TRAC"))
cowplot::plot_grid(Gene_UMAPS_1,Gene_UMAPS_2,Gene_UMAPS_3,Gene_UMAPS_4)
ggsave("Gene_UMAPS.png", units = "cm", width = 20, height = 10, scale = 2)
