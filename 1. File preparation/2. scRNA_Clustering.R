library(tidyverse)
library(Seurat)
library(SeuratObject)
library(readxl)
library(openxlsx)
library(conflicted)
conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)
library(future)
library(SingleR)
library(Seurat.utils)
library(reticulate)
library(patchwork)

future::plan("multisession", workers = 8)
options(future.globals.maxSize= 25*1024^3)

scdata<-read_rds('CVD_Total.rds') # Read file from the post-integration

scdata<-UpdateSeuratObject(scdata)

DefaultAssay(scdata)<-'SCT'
FeatureScatter(scdata, feature1 = "percent.mt", feature2 = "percent.ribo")

#QC----
VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), pt.size = 0, 
        raster=T, ncol = 3)

table(scdata@meta.data$Sample_Name, scdata@meta.data$Group)

#Recluster----
features <- c('CD3D', 'CD3E', #T cell
              "CD4", #CD4 T cell and monocytes
              "CD8A", 'CD8B', #CD8 T cell and some NK
              'CCR7', #Naive T cells
              'CD14', "LYZ", #Classic monocytes
              "FCGR3A", #CD16, nonclassical monocytes and NK
              'NCAM1', 'PRF1', #NK
              'CD1C', #Conventional dendritic cells
              'IL3RA', #Plasmocytoid dendritic cells
              'MS4A1', #B cells CD20
              'CD38', #Plasma cells
              'PPBP', #Platlet
              'CD34', #Progenitor
              'MKI67' #Proliferation marker 
              )

FeaturePlot(scdata, features, raster = T, pt.size = 2, ncol = 4, order = F)

DefaultAssay(scdata)<-'SCT'

(DotPlot(scdata, features, assay='SCT', scale = T, group.by = 'seurat_clusters', cluster.idents = T) |
  DimPlot(scdata, group.by = 'seurat_clusters', label = T, raster = T)) +
  plot_layout(widths = c(2,1))

label1<-table(scdata@meta.data$seurat_clusters, scdata@meta.data$Cell_Type_Experimental) %>% as.data.frame()%>%
  as_tibble()%>%filter(Freq>10)%>%arrange(Var1, -Freq)

#Final cell identity----
cell_ident<-read_xlsx('Cell_Type.xlsx')%>%mutate(seurat_clusters=factor(seurat_clusters))

FeaturePlot(subset(scdata, subset= seurat_clusters%in% c('20', '26')), 
        c('CD3D', 'CD4', 'CD8A', 'FCGR3A', 'LINGO2'), label = F, raster = F)

VlnPlot(scdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), raster = T, flip=T, stack = T)

table(cell_ident$Cell_Type_Main_Final)

#Assign cell identity----
metadata<-scdata@meta.data%>% select(-Cell_Type_Main_Final)%>%
  rownames_to_column('Barcode')%>%
  left_join(cell_ident%>%select(seurat_clusters, Cell_Type_Main_Final), by='seurat_clusters')%>%
  select(Barcode, Cell_Type_Main_Final)%>%
  column_to_rownames('Barcode')%>%
  mutate(Cell_Type_Main_Final=factor(Cell_Type_Main_Final, 
                                     levels = c('CD4+ T', 'CD8+ T', 'T', 'NK', "T Prolif", 'NK and T LINGO2+', 'B', 'Plasmablast',
                                                'CD14+ Monocytes', 'CD16+ Monocytes', 'cDC', 'pDC',   
                                                'Progenitor', 'Platelet')))
base::setdiff(cell_ident$Cell_Type_Main_Final, unique(metadata$Cell_Type_Main_Final))

scdata<-AddMetaData(scdata, metadata = metadata)

DimPlot(scdata, group.by = 'Cell_Type_Main_Final', label = T, raster=T)  

FeaturePlot(scdata, features = c('CD3D', 'CD4', 'CD8A', 'FCGR3A', 'MKI67', 'GNLY', 'PPBP', 'LINGO2'), raster = T)

DefaultAssay(scdata)<-'SCT'
VlnPlot(scdata, features = features, pt.size = 0, stack = T, flip = T,
        raster=F, group.by = 'Cell_Type_Main_Final') & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none')

DotPlot(scdata, features = features,  
        group.by = 'Cell_Type_Main_Final', scale = F, dot.min = .1)

table(scdata@meta.data$Cell_Type_Main_Final, scdata@meta.data$Cell_Type_Experimental)

table(scdata@meta.data$Cell_Type_Main_Final, scdata@meta.data$monaco.main)

DimPlot(scdata, label = T, raster = T, pt.size = 2, raster.dpi = c(600,600))+
  DimPlot(scdata, label=T, raster=T, pt.size=2, group.by = 'Cell_Type_Main_Final')

FeaturePlot(scdata, c('CD8A', 'GNLY'), raster = T, blend = T)

#Decipher cluster 26----
future::plan("sequential")
scdata<-PrepSCTFindMarkers(scdata)

DefaultAssay(scdata)
marker_26<-FindMarkers(scdata, ident.1 = '26', ident.2 = '4', 
                       assay = 'SCT', logfc.threshold = .5, min.pct = 0.1)

FeaturePlot(subset(scdata, subset= seurat_clusters %in% c('26')), 
            features = c('CD8A', 'NKG7'), blend = T,
            raster = T, pt.size =2)

DotPlot(subset(scdata, subset= seurat_clusters %in% c('4', '16', '26')), scale = T, col.min = 0, dot.min = .25,
        features = c('CD4', 'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CCR7', 'FCGR3A', 'IL23R', 'SELL', 'TCF7', 'LEF1', 'FAS'))

VlnPlot(subset(scdata, subset= seurat_clusters %in% c('4', '26')), stack = T, flip = T, 
        features = c('CD4', 'CD8A', 'CD8B', 'CD3D', 'CD3E', 'CCR7', 'FCGR3A', 'IL23R', 'SELL', 'FAS'))+
  stat_summary(fun=median, geom="point", shape=18, size=3, color="grey")

#Cluster 26 is a mix of CD4 and CD8 naive 

#Decipher cluster 20----
marker_20<-FindMarkers(scdata, ident.1 = '20', ident.2 = '25', 
                       assay = 'SCT', logfc.threshold = .5, min.pct = 0.1)

DotPlot(subset(scdata, subset= seurat_clusters %in% c('20', '25')), scale = F, col.min = 0, dot.min = .25,
        features = c('CD4', 'CD8A', 'NKG2A', 'NKG7', 'FCGR3A', 'LINGO2'))

VlnPlot(subset(scdata, subset= seurat_clusters %in% c('20', '25')), pt.size = 0, stack = T, flip = T,
        features = c('CD4', 'CD8A', 'NKG2A', 'NKG7', 'FCGR3A', 'LINGO2'))

#Cluster 20 should be NK

#Add metadata----
group<-read_xlsx('Grouping.xlsx')%>%select(Sample_Name, MACSID, US_Risk, Batch)

metadata<-scdata@meta.data%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name, Group, orig.ident)%>%
  left_join(group, 'Sample_Name')%>%
  column_to_rownames('Barcode')

table(metadata$Group, metadata$US_Risk)

scdata<-AddMetaData(scdata, metadata)

#Final cell identity----
DefaultAssay(scdata)<-'SCT'
Idents(scdata)<-'Cell_Type_Main_Final'

DimPlot(scdata, label = T, group.by = 'seurat_clusters')

(FeaturePlot(subset(scdata, subset= seurat_clusters %in% c('11', '4', '24', '0', '10', '6', '20', '25')), features =
             c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD8B', 'FCGR3A', 'GNLY', 'NCAM1', 'PRF1', 'IL32', 'KLRD1'), 
             min.cutoff = 'q5', max.cutoff = 'q75', raster = T, pt.size = 2, raster.dpi = c(600,600),
             order = F) | 
  DimPlot(scdata, label=T, group.by = 'seurat_clusters'))+plot_layout(widths = c(2,1))

VlnPlot(subset(scdata, Cell_Type_Main_Final=='CD4+ T'), group.by = 'seurat_clusters', features = c('CD3D', 'CD4', 'CD8A','CD8B',  'FCGR3A', 'MKI67', 'GNLY', 'PPBP'), pt.size = 0, stack = T, flip = T,
        raster=F) & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none')

VlnPlot(subset(scdata, Cell_Type_Main_Final=='CD8+ T'), group.by = 'seurat_clusters', 
        features = c('CD3D', 'CD4', 'CD8A', 'CD8B', 'FCGR3A', 'NCAM1', 'KLRD1', 'GNLY', 'IL32'), pt.size = 0, stack = T, flip = T,
        raster=F) & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none') |

VlnPlot(subset(scdata, Cell_Type_Main_Final=='NK'), group.by = 'seurat_clusters', 
        features = c('CD3D', 'CD4', 'CD8A', 'CD8B', 'FCGR3A', 'NCAM1', 'KLRD1', 'GNLY', 'IL32'), pt.size = 0, stack = T, flip = T,
        raster=F) & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none')

VlnPlot(subset(scdata, Cell_Type_Main_Final=='B'), group.by = 'seurat_clusters', features = c('CD3D', 'MS4A1'), pt.size = 0, stack = T, flip = T,
        raster=F) & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none')

VlnPlot(subset(scdata, Cell_Type_Main_Final=='B'), group.by = 'seurat_clusters', features = c('CD3D', 'MS4A1'), pt.size = 0, stack = T, flip = T,
        raster=F) & theme(axis.text.y = element_text(size=8, angle = 0)) &
  guides(fill='none')

#Save data----
save.image('CVD_scdata.Rdata')