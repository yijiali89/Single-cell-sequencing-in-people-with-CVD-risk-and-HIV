#Optional if you want to use future package if you have limited RAM
library(future)
availableCores()
parallel::detectCores()
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(dplyr)
library(ggplot2)
library(future)
plan("multisession", workers = 12)  
options(future.globals.maxSize = 300 * 1024^3) 

#Read four batches of files----
s1 <- readRDS('CVD1_singlets_Sept10.rds')
s2 <- readRDS('CVD2_singlets_Sept10.rds')
s3 <- readRDS('CVD3_singlets_Sept10.rds')
s4 <- readRDS('CVD4_singlets_Sept10.rds')

# Remove sampletag02 which belongs to neg at the time of joining the cohort
# Of note this participant, at time of sampling, was seronegative
s3 <- subset(s3, subset = Sample_Tag != "SampleTag02_hs")

table(s1$Sample_Tag)
table(s2$Sample_Tag)
table(s3$Sample_Tag)
table(s4$Sample_Tag)
gc()

#Merging 4 batches----
#Of note the files have already gone through initial QC including doublet removal, hemoglobin removal from 
#BD software.
CVD <- merge(s1, y = c(s2, s3, s4), add.cell.ids = c("Sample1", "Sample2", "Sample3", "Sample4"))

#Reassign Group column using Sample_Tag in merged object
CVD$Group <- case_when(
  CVD$Sample_Tag %in% c("SampleTag07_hs", "SampleTag08_hs", "SampleTag05_hs", "SampleTag01_hs", "SampleTag10_hs") ~ "Low",
  CVD$Sample_Tag %in% c("SampleTag09_hs", "SampleTag04_hs", "SampleTag06_hs", "SampleTag03_hs", "SampleTag11_hs", "SampleTag12_hs") ~ "High",
  TRUE ~ "Not_allocated"
)

table(CVD$Group)

CVD@meta.data <- CVD@meta.data[, !grepl("percent.mt", colnames(CVD@meta.data))]
CVD@meta.data <- CVD@meta.data[, !grepl("percent_mt", colnames(CVD@meta.data))]
CVD@meta.data <- CVD@meta.data[, !grepl("percent.ribo", colnames(CVD@meta.data))]
CVD@meta.data <- CVD@meta.data[, !grepl("Cell_Type_Main_Final", colnames(CVD@meta.data))]

DefaultAssay(CVD)<-'RNA'
CVD[["percent.mt"]] <- PercentageFeatureSet(object = CVD, pattern = "^MT-")
CVD <- PercentageFeatureSet(CVD, "^RP[SL]", col.name = "percent.ribo")

VlnPlot(CVD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = .01, 
        raster=T, ncol = 4)

FeatureScatter(CVD, group.by ='orig.ident' , feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(CVD, group.by ='orig.ident' , feature1 = "nFeature_RNA", feature2 = "percent.ribo")
FeatureScatter(CVD, group.by ='orig.ident' ,feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#Remove low quality "cell/barcode" using the following criteria
CVD<-subset(CVD, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20 & percent.ribo>5) 
rm(s1, s2, s3, s4); gc() #Save up some RAM space

VlnPlot(CVD, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = .01, 
        raster=T, ncol = 4)
 
gc()
 
table(CVD$Sample_Name, CVD$Group)

CVD[["RNA"]] <- split(CVD[["RNA"]], f = CVD$orig.ident)


#SCT transformation v2 from Seurat 5 that remove batch effect using IntegratedLayers----
CVD <- SCTransform(CVD)
CVD <- RunPCA(CVD, npcs = 30, verbose = T)

CVD <- IntegrateLayers(
  object = CVD,
  method = RPCAIntegration,
  normalization.method = "SCT",
  verbose = T
)

DimPlot(CVD, reduction = 'integrated.dr')
CVD <- FindNeighbors(CVD, dims = 1:30, reduction = "integrated.dr")
CVD <- FindClusters(CVD, resolution = 1, algorithm=2)

CVD <- RunUMAP(CVD, dims = 1:30, reduction = "integrated.dr", return.model = T)

DimPlot(CVD)

CVD<- JoinLayers(CVD, assay = 'RNA')
DimPlot(CVD, group.by = 'Cell_Type_Experimental')
DimPlot(CVD, split.by = 'orig.ident')

DefaultAssay(CVD)

saveRDS(CVD, file = "CVD_total.rds")
