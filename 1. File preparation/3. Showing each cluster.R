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

load('CVD_scdata.Rdata')

DefaultAssay(scdata)<-'SCT'
Idents(scdata)<-'Cell_Type_Main_Final'

DimPlot(scdata, label = T, group.by = 'seurat_clusters')
 
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
 
#Save data----
save.image('CVD_scdata.Rdata')
saveRDS(scdata, 'CVD_total.rds')