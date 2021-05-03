if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("fishpond")
BiocManager::install("EnsDb.Hsapiens.v79")

library(Seurat)
library(tximport)
library(SeqGSEA)
#library(biomaRt)
library(fishpond)
#library(dplyr)
library(EnsDb.Hsapiens.v79)

files <- file.path("/projectnb/bf528/users/frazzled/project_4/Data_curator/salmon/salmon_output/alevin/quants_mat.gz")
file.exists(files)
txi <- tximport(files, type="alevin")
t <- txi$counts



# get ensembl ids, and remove extensions
ensembl <-  rownames(t) #get ensmbl ids
ensembl <- as.character(genes)
ensembl <- sub("[.][0-9]*","",ensembl) #remove ends
rownames(t) <- ensembl

# convert gene names of ensembl ids 
symbols <- select(EnsDb.Hsapiens.v79, keys= ensembl, keytype = "GENEID", columns = c("SYMBOL","GENEID")) 

# subset to rows that have the matches, and rename matrix with genenames
t <- t[rownames(t) %in% symbols$GENEID,]
rownames(t) <- symbols$SYMBOL

pbmc <- CreateSeuratObject(counts = t , min.cells = 3, min.features = 200, project = "10X_PBMC")
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
head(pbmc@meta.data, 20)

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

#filtering 
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 35)

pbmc

#normalize data 
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc[["RNA"]]@meta.features <- data.frame(row.names = rownames(pbmc[["RNA"]]))

#find variable features 
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)


#scale data 
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features=all.genes)

#perform PCA on scaled data 
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

#view PCA results
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
DimPlot(pbmc, reduction = "pca")
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#determine dimensionality of sample set
#resampling test
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:15)
ElbowPlot(pbmc)

#cluster the cells 
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#view cluster ids first 5 cells 
head(Idents(pbmc), 5)
pbmc




saveRDS(pbmc, file = "/projectnb/bf528/users/frazzled/project_4/programmer/pbmc_seurat2.rds")


counts <- as.vector(table(Idents(pbmc)))
names <- names(table(Idents(pbmc)))
png('pie.png')
slices <- as.vector(prop.table(table(Idents(pbmc))))
lbls <- names(prop.table(table(Idents(pbmc))))
pie(slices, labels = lbls, main="relative proportions")
dev.off()

png('barplot.png')
bar <- barplot(height = table(Idents(pbmc)), names.arg = (names), xlab = 'Cluster', ylab = 'total cells')
text(x = bar, y = table(Idents(pbmc))-60, label = table(Idents(pbmc)), pos = 3, cex = 0.8, col = "blue")
dev.off()
