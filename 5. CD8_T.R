library(tidyverse)
library(Seurat)
library(SeuratObject)
library(readxl)
library(openxlsx)
library(conflicted)
library(Seurat.utils)
library(glmGamPoi)
library(DESeq2)
library(fgsea)
library(patchwork)
library(scplotter)
library(monocle3)
library(stringr)
library(clusterProfiler)

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

#load('CD8_T.rdata')
#t<-readRDS('CVD_T.rds')

#Recluster----
DimPlot(t)

DefaultAssay(t)<-'SCT' #Rerun clustering for t
t<-RunPCA(t)
ElbowPlot(t, ndims = 30)
DimHeatmap(t, dims = 1:20, cells = 500, balanced = TRUE)

t<-RunUMAP(t, dims = 1:20)

t<-FindNeighbors(t, dims = 1:20)
t<-FindClusters(t, algorithm = 2, random.seed = 121, resolution = 1)

Idents(t)<-'seurat_clusters'
DefaultAssay(t)<-'SCT'
DimPlot(t, label = TRUE, label.size = 3)
DimPlot(subset(t, subset= Cell_Type_Main_Final=='CD4+ T'), label = TRUE, label.size = 3)+FeaturePlot(t, c('MKI67', 'CD4'))

#Identification----
DotPlot(object =t, scale = T, assay = 'SCT',  cluster.idents = T, scale.by = 'size', 
        dot.min = .05, col.min = -1, col.max = 3, cols = c('grey99', 'firebrick'),
        features = c('CD3D', 'CD4', 'CD8A', 'CCR7', 'FCGR3A', 'SELL', 'CD27', 'FAS', 'CD28', 'IL7R', 'CCL5','GZMH', 'CXCR3', 'ICOS', 'MAF', 'PDCD1',
                     'S100A11', 'TIMP1', 'CST7', 'GZMK', 'GZMA', 'FOXP3', 'CTLA4', 'IL2RA', 'IFNG', 'HLA-DRB1', 'percent.mt')) /
  (FeaturePlot(t, features = c('CD4', 'CD8A', 'FAS', 'CCR7'))|DimPlot(t, label = TRUE, label.size = 3))

label<-table(t@meta.data$Cell_Type_Main_Final, t@meta.data$seurat_clusters)%>%
  as.data.frame()%>%mutate(Var2=as.character(Var2))%>%
  group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq==Freqmax)

VlnPlot(t, c('CD3D', 'CD3E', 'CD4', 'CD8A', 'FCGR3A'), pt.size = 0, stack = T, flip = T)|DimPlot(t, label = TRUE, label.size = 5)

cd8_t<-subset(t, subset=seurat_clusters %in% c(label$Var2[label$Var1=='CD8+ T'], '20', '26'))

#CD8----
DimPlot(cd8_t, group.by = 'Cell_Type_Main_Final')
cd8_t<-subset(cd8_t, Cell_Type_Main_Final %in% c('CD4+ T', 'NK', 'NK and T LINGO2+')==FALSE)
DimPlot(cd8_t, group.by = 'Cell_Type_Main_Final')

FeaturePlot(cd8_t, c('CD4', 'CD8A', 'CD3D'))| DimPlot(cd8_t)

hist(unlist(FetchData(cd8_t, c('CD8A'))), breaks = 100)
hist(unlist(FetchData(cd8_t, c('CD8B'))), breaks = 100)
hist(unlist(FetchData(cd8_t, c('CD4'))), breaks = 100)
hist(unlist(FetchData(cd8_t, c('FCGR3A'))), breaks = 100)

hist(unlist(FetchData(cd8_t, c('CD3D'))), breaks = 100)
hist(unlist(FetchData(cd8_t, c('CD3E'))), breaks = 100)

cd8_t<-subset(cd8_t, subset= CD8A>=CD4 & CD8B>=CD4) #Filter those CD4 high cells
hist(unlist(FetchData(cd8_t, c('CD8A'))), breaks = 1000)
hist(unlist(FetchData(cd8_t, c('FCGR3A'))), breaks = 1000)

DefaultAssay(cd8_t)<-'SCT' #Rerun clustering for t 
cd8_t<-RunPCA(cd8_t)
ElbowPlot(cd8_t)
DimHeatmap(cd8_t, dims = 1:10, cells = 1000, balanced = TRUE)

cd8_t<-RunUMAP(cd8_t, dims = 1:10)

cd8_t<-FindNeighbors(cd8_t, dims = 1:10)

cd8_t<-FindClusters(cd8_t, algorithm = 2, random.seed = 121, resolution = 0.25)

DimPlot(cd8_t, label = TRUE, label.size = 5)

DefaultAssay(cd8_t)<-'SCT'

VlnPlot(cd8_t, features = c('CD3E', 'CD8A', 'FCGR3A', 'NCAM1', 'GNLY', 'CD4', 'CD8A', 'CD8B', 'FCGR3A', 'MKI67', 'FAS', 'CCR7', 'CD27', 'SELL'), 
        pt.size = 0, stack=T, flip=T) |
FeaturePlot(cd8_t, features = c('CD4', 'CD8A', 'CD8B', 'IL7R', 'GZMA', 'GNLY', 'GZMB', 'TCF7', 'MKI67', 'TRDC'), order=T)+
  DimPlot(cd8_t, label = TRUE, label.size = 5)

Idents(cd8_t)<-'seurat_clusters'
VlnPlot(cd8_t, features = c('CD3E', 'CD8A','CD8B', 'GNLY', 'KIR2DL1', 'KIR2DL3', 
                            'KIR3DL1', 'KIR3DL2', 'KLRC3', 'KLRC2', 'KLRC1',
                            'NKG7', 'GNLY', 'GZMB', 'PRF1', 'TCF7', 'SELL', 'BACH2'), 
        pt.size = 0, stack=T, flip=T) +
  DimPlot(cd8_t, label = TRUE, label.size = 5)

FeaturePlot(cd8_t,features = c('CD3E', 'CD8A','CD8B', 'GNLY', 'KIR2DL1', 'KIR2DL3', 
                               'KIR3DL1', 'KIR3DL2', 'KLRC3', 'KLRC2', 'KLRC1',
                               'NKG7', 'GNLY', 'GZMB', 'PRF1', 'TCF7', 'SELL', 'BACH2'))

FeaturePlot(cd8_t, features = c('CD3E', 'CD8A','CD8B', 'GNLY', 'KIR2DL1', 'KIR2DL3', 
                                'KIR3DL1', 'KIR3DL2', 'KLRC3', 'KLRC2', 'KLRC1',
                                'NKG7', 'GNLY', 'GZMB', 'PRF1', 'TCF7', 'SELL', 'BACH2'), 
            order=T)+
  DimPlot(cd8_t, label = TRUE, label.size = 5)

DefaultAssay(cd8_t)<-'SCT'

cd8_t@assays$SCT@SCTModel.list$model1.1<-NULL
cd8_t@assays$SCT@SCTModel.list$model1.2<-NULL
cd8_t@assays$SCT@SCTModel.list$model1.3<-NULL

cd8_t<-PrepSCTFindMarkers(cd8_t)

cd8_t.markers <- FindAllMarkers(cd8_t, only.pos = F, logfc.threshold = .5, min.pct = .1)

write.xlsx(cd8_t.markers, 'Fig 4A.xlsx')

cd8_t.markers<- cd8_t.markers%>%filter(p_val_adj<0.05)

VlnPlot(cd8_t, c('nFeature_RNA', 'nCount_RNA'), pt.size = 0, stack = T, flip = T, split.by = 'orig.ident')

cd8_features<- c('S100B', 'CCR7', 'LEF1', 'TCF7', 'SELL', 
                 'IL7R', 'CD27', 'CD28', 'MKI67', 'FAS', 'CCL5', 'GZMH', 'GZMB', 'GZMA', 
                 'GNLY', 'TRDC', 'KLRB1', 'PDCD1')

Idents(cd8_t)<-'seurat_clusters'
(VlnPlot(cd8_t, features = cd8_features, 
        pt.size = 0, stack = T, flip = T) & guides(fill='none')) + DimPlot(cd8_t, label = TRUE, label.size = 3) /
  DotPlot(object =cd8_t, scale = F, cluster.idents = T,
          features = cd8_features) + plot_layout(widths = c(1, 2))

(VlnPlot(subset(cd8_t, seurat_clusters %in% c('8', '17','16', '11', '12', '15')==FALSE), features = cd8_features, 
         pt.size = 0, stack = T, flip = T) & guides(fill='none')) + DimPlot(cd8_t, label = TRUE, label.size = 3) /
  DotPlot(object =subset(cd8_t, seurat_clusters %in% c('8', '17', '16', '11', '12', '15')==FALSE), scale = F, 
          cluster.idents = T,
          features = cd8_features)

#Rename cluster----
metadata_cd8t<-cd8_t@meta.data
library(SingleR)
library(SingleCellExperiment)
library(celldex)
sce <- as.SingleCellExperiment(cd8_t)
ref.set <- celldex::MonacoImmuneData()
pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.fine)
lbls.keep <- table(pred.cnts$labels)>10

cd8_t$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')

DimPlot(cd8_t, label=T, label.size = 4, group.by='SingleR.labels')+
  DimPlot(cd8_t, label=T, label.size = 4)

label_cd8<-table(cd8_t@meta.data$SingleR.labels, cd8_t@meta.data$seurat_clusters)%>%
  as.data.frame()%>%group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq>=Freqmax/2)

Idents(cd8_t)<-'seurat_clusters'

FeaturePlot(cd8_t, 'PDCD1', order=T)
VlnPlot(cd8_t, cd8_features, flip=T, stack=T)

cd8_features<- c('CD8A', 'CD8B', 'S100B', 'CCR7', 'LEF1', 'TCF7', 'SELL', 
                 'IL7R', 'CD27', 'CD28', 'FAS', 'CCL5', 
                 'GZMH', 'GZMB', 'GZMA', 'PRF1',
                 'GNLY', 'TRDC', 'KLRB1', 'MKI67')

Idents(cd8_t)<-'seurat_clusters'
(VlnPlot(cd8_t, features = cd8_features, 
         pt.size = 0, stack = T, flip = T) & guides(fill='none')) + 
  (DimPlot(cd8_t, label = TRUE, label.size = 3) /
  (FeaturePlot(object =cd8_t, features = cd8_features)  & NoLegend() & NoAxes())) + 
  plot_layout(widths = c(1, 2))

Clusters_Final_CD8<-cd8_t@meta.data%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters %in% c('3') ~ 'CD8_TN', #CD8 TN 
                                  seurat_clusters %in% c('2') ~ 'CD8_TCM', #CCR7+,  IL7R+
                                  seurat_clusters %in% c('4') ~ 'CD8_MAIT', #KLRB1-high SLC4A10-high
                                  seurat_clusters %in% c('5') ~ 'gdT', #TRDC
                                  seurat_clusters %in% c('6') ~ 'CD8_Prolif', #high KI67
                                  seurat_clusters %in% c('0') ~ 'CD8_TTEM', #Lowest IL7R and TCF7, highest GZMH 
                                  seurat_clusters %in% c('1') ~ 'CD8_TEM' 
  ))%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('CD8_TN', 'CD8_TCM','CD8_TEM', 'CD8_TTEM', 'CD8_MAIT', 'gdT', 'CD8_Prolif')))%>%
  select(Barcode, Clusters_Final, seurat_clusters)%>%
  column_to_rownames('Barcode')

cd8_t<-AddMetaData(cd8_t, Clusters_Final_CD8)

DimPlot(cd8_t, group.by = 'seurat_clusters', label=T)
DimPlot(cd8_t, group.by = 'Clusters_Final', label=T)

FeaturePlot(cd8_t, c('CD3D', 'CD3E', 'CD8A', 'CD8B', 'TRDC', 'GNLY', 'FCER1G'), order=T)

Idents(cd8_t)<-'Clusters_Final'

scplotter::CellDimPlot(cd8_t, label = T, theme = "theme_blank", palette='locuszoom', dot_size=1, 
                       label_insitu = T, pt_alpha =0.5)

Idents(cd8_t)<-'Clusters_Final'
DotPlot(object =cd8_t, scale = T, assay = 'SCT',  
        cluster.idents = T, dot.scale = 5, scale.by='size',  col.min = -1.5,
        cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
        features = c(cd8_features)) +theme_bw()+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab('')

VlnPlot(cd8_t, stack = T, 
        features = cd8_features, 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) 

DefaultAssay(cd8_t)<-'SCT'
(FeaturePlot(cd8_t, features = c('CD8A', 'CD8B', 'CCR7', 'IL7R', 'MKI67', 'CCL5', 'GZMH'), raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, #min.cutoff = 0.1, max.cutoff = 3, 
            cols = c('grey95', 'firebrick'), slot = 'data', #split.by = 'Group', 
            reduction = 'umap') & NoAxes())| DimPlot(cd8_t)

(FeaturePlot(cd8_t, features = c('CCR7', 'IL7R','CCL5', 'GZMH', 'TRDC', 'MKI67'), raster = T, order = T,
             raster.dpi = c(600, 600), pt.size = 3, #min.cutoff = 0.1, max.cutoff = 3, 
             cols = c('grey95', 'firebrick'), slot = 'data',  split.by = 'Group', 
             reduction = 'umap') & NoAxes())

DotPlot(cd8_t, features = c('LEF1', 'CCR7', 'FHIT', 'TCF7', 'SELL', 'FAS', 'GZMA', 'FOXP3', 'IL2RA', 'CTLA4', 'PDCD1'), scale =T, cluster.idents = T,
        dot.min = .05, dot.scale = 6, scale.by = 'size', col.max = 1, col.min = -1)+ theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = 'grey90'), panel.grid.minor = element_line(color = 'grey90'))

(FeaturePlot(subset(cd8_t, Group=='High'), features = c('CD4', 'CCR7', 'FAS', 'GZMH', 'FOXP3', 'MKI67'), raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#BC3C29CC'), slot = 'data', 
            reduction = 'umap') & NoAxes()) |
  
  (FeaturePlot(subset(cd8_t, Group=='Low'), features = c('CD4', 'CCR7', 'FAS', 'GZMH', 'FOXP3', 'MKI67'), raster = T, order = T,
               raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#0072B5CC'), slot = 'data', 
               reduction = 'umap') & NoAxes())

#Compare component of cluster and disease----
metadata_cd8<-cd8_t@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(cd8_t, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar") 

scCustomize::Proportion_Plot(seurat_object = subset(cd8_t, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar")

cellpercent<-metadata_cd8%>% #summarize percentage of each cluster_final between control and disease group
  group_by(Group, Sample_Name, Clusters_Final, .drop=FALSE)%>%
  summarise(Count=n())%>%
  mutate(Percentage=Count/sum(Count)*100)%>%
  ungroup()

cellpercent%>%
  pivot_wider(id_cols=Sample_Name, names_from = c(Clusters_Final, Group), values_from = Percentage)%>%
  write.xlsx('tmp.xlsx')

cellpercent%>%
  ggplot(aes(x=Group, y=Percentage, color=Group))+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  geom_jitter(aes(color=Group), width=0.2)+
  scale_color_manual(values = c("#BC3C29CC", "#0072B5CC"))+
  facet_wrap(~Clusters_Final) +theme_bw()+stat_pwc() + scale_y_continuous(limits = c(0, 100))

cvd_metadata<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Data_CVD.xlsx')%>%select(MACSID, VL, US_Risk, USLESION:CD4NADIR), by=c('MACSID', 'US_Risk'))

#Add CVD data to metadata----
metadata_cvd_add<-metadata_cd8%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

cd8_t<-AddMetaData(cd8_t, metadata_cvd_add)

#Calculate CVD correlation----
cvd<-cvd_metadata%>%
  right_join(cellpercent, by='Sample_Name')%>%
  select(Clusters_Final, Sample_Name, Group, Percentage, USLESION:CD4NADIR)

cor_res1<-cvd%>%
  select(Clusters_Final, Sample_Name, Percentage, IMT, contains('SCORE'), contains('COUNT'), MMP2:CD4NADIR)%>%
  pivot_longer(cols = c(IMT, contains('SCORE'), contains('COUNT'), MMP2:CD4NADIR), 
               names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  group_by(Clusters_Final, Var)%>%
  group_modify(.f=~broom::tidy(cor.test(.x$Percentage, .x$Value, method = 'spearman', data = .x))) %>%
  ungroup()%>%
  mutate(Padj=p.adjust(p.value, 'BH')) 

cor_res2<-cvd%>%
  select(Clusters_Final, Sample_Name, Percentage, Group, AS_DI:MP_DI, STEN4, STEN34)

cor_res2[, 4:11]<-map_df(cor_res2[, 4:11], as.factor)

cor_res2<-cor_res2%>%
  pivot_longer(cols = c(Group, contains('DI')), 
               names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  group_by(Clusters_Final, Var)%>%
  group_modify(.f=~broom::tidy(wilcox.test(.x$Percentage~.x$Value, data = .x))) %>%
  ungroup()%>%
  mutate(Padj=p.adjust(p.value, 'BH'))

#DEG----
DefaultAssay(cd8_t)<-'SCT'

##Total CD8 ----
cd8_total_for_deseq<-subset(cd8_t, Clusters_Final %in% c('CD8_TN', 'CD8_TCM', 'CD8_TEM', 'CD8_TTEM', 'CD8_MAIT'))

DimPlot(cd8_total_for_deseq)

pseudo_total <- AggregateExpression(cd8_total_for_deseq, assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(cd8_total_for_deseq)[2]

pseudo_total<-pseudo_total$RNA

pseudo_total<-pseudo_total%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_total)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_total) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_total,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("GMT/h.all.v2025.1.Hs.symbols.gmt")
pathway.reactome<- gmtPathways("GMT/c2.cp.reactome.v2025.1.Hs.symbols.gmt")
pathway.gobp<-gmtPathways("GMT/c5.go.bp.v2025.1.Hs.symbols.gmt")

set.seed(111)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

pdf('Fig 5E.pdf', height=2.5, width = 6)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) + 
  scale_size_continuous(breaks = c(50, 100,150, 200), limits = c(50, 200))
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_MYC_TARGETS_V1']], stats = ranks)+
  ggtitle('HALLMARK_MYC_TARGETS_V1')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

###MYC gene modules----
enrich_df <- data.frame(ID = 'GO:0000001',
                        Description = 'HALLMARK_MYC_TARGETS_V1',
                        setSize = fgseaRes$size[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1'],
                        enrichmentScore = fgseaRes$ES[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1'],
                        NES = fgseaRes$NES[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1'],
                        pvalue = fgseaRes$pval[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1'],
                        p.adjust = fgseaRes$padj[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1'],
                        core_enrichment = paste0(unlist(fgseaRes$leadingEdge[fgseaRes$pathway=='HALLMARK_MYC_TARGETS_V1']), collapse = '/'))

rownames(enrich_df) <- 'GO:0000001'

en <- GseaVis::dfGO2gseaResult(enrich.df = enrich_df,keytype = 'SYMBOL',
                               geneList = sort(ranks, decreasing = T),
                               OrgDb = 'org.Hs.eg.db')

library(org.Hs.eg.db)
en <- DOSE::setReadable(x = en, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")

en@geneSets$`GO:0000001`<- pathways.hallmark$HALLMARK_MYC_TARGETS_V1

enrichplot::gseaplot2(x = en, geneSetID = enrich_df$ID[1], color = 'grey20', pvalue_table = TRUE,
                      title = enrich_df$Description[1], subplots = 1:2) 

myc_genes<-pathways.hallmark$HALLMARK_MYC_TARGETS_V1
cd8_t$MYC1<-NULL

cd8_t<-AddModuleScore(cd8_t, features = list(myc_genes), name = 'MYC')

###Viral infection gene modules----
vi_genes<-pathway.reactome$REACTOME_VIRAL_INFECTION_PATHWAYS

cd8_t<-AddModuleScore(cd8_t, features = list(vi_genes), name = 'VI')

FeaturePlot(cd8_t, 'VI1', split.by = 'Group')

###IFN gene modules----
IFNA_genes<-pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE
IFNG_genes<-pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE

cd8_t<-AddModuleScore(cd8_t, features = list(IFNA_genes), name = 'IFNA_Mod')
cd8_t<-AddModuleScore(cd8_t, features = list(IFNG_genes), name = 'IFNG_Mod')

###Estrogen gene modules----
ES_early_genes<-pathways.hallmark$HALLMARK_ESTROGEN_RESPONSE_EARLY
cd8_t<-AddModuleScore(cd8_t, features = list(ES_early_genes), name = 'ES_Early')

###Apoptosis gene modules----
Apop_genes<-pathways.hallmark$HALLMARK_APOPTOSIS
cd8_t$Apop1<-NULL

cd8_t<-AddModuleScore(cd8_t, features = list(Apop_genes), name = 'Apop', search = T)

FeaturePlot(cd8_t, 'Apop1', split.by = 'Group')
FeaturePlot(cd8_t, 'MYC1', split.by = 'Group')

VlnPlot(cd8_t, c('MYC1', 'Apop1'), split.by = 'Group', pt.size = 0, stack=T, flip=T) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))+
  scale_y_continuous(expand = c(0.25,0))

(VlnPlot(cd8_t, 'IL12_1', split.by = 'Group', pt.size = 0) +
    stat_compare_means(label = 'p.format')+
    scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
    stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
    stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9)))/

(VlnPlot(cd8_t, 'VI1', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))) /

(VlnPlot(cd8_t, 'HIV1', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9)))

VlnPlot(cd8_t, 'BCL3', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))
 

##CD8 TN----
pseudo_tn <- AggregateExpression(subset(cd8_t, Clusters_Final=='CD8_TN'), assays = "RNA", 
                                    slot = 'counts', return.seurat = F, 
                                    group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd8_t, Clusters_Final=='CD8_TN'))[2]

pseudo_tn<-pseudo_tn$RNA

pseudo_tn<-pseudo_tn%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_tn)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_tn) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_tn,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '2. Padj 0.05-0.1'= '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.1)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

deg_list<-list()

deg_list$TN<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


##CD8 TCM----
pseudo_tcm <- AggregateExpression(subset(cd8_t, Clusters_Final=='CD8_TCM'), assays = "RNA", 
                                 slot = 'counts', return.seurat = F, 
                                 group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd8_t, Clusters_Final=='CD8_TCM'))[2]

pseudo_tcm<-pseudo_tcm$RNA

pseudo_tcm<-pseudo_tcm%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_tcm)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_tcm) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_tcm,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '2. Padj 0.05-0.1'= '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

deg_list$TCM<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


##CD8 TEM----

pseudo_tem <- AggregateExpression(subset(cd8_t, Clusters_Final =='CD8_TEM'), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd8_t, Clusters_Final =='CD8_TEM'))[2]

pseudo_tem<-pseudo_tem$RNA

pseudo_tem<-pseudo_tem%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_tem)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_tem) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_tem,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '2. Padj 0.05-0.1'= '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

deg_list$TEM<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


##CD8 MAIT----
pseudo_mait<- AggregateExpression(subset(cd8_t, Clusters_Final=='CD8_MAIT'), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd8_t, Clusters_Final=='CD8_MAIT'))[2]

pseudo_mait<-pseudo_mait$RNA

pseudo_mait<-pseudo_mait%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_mait)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_mait) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_mait,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '2. Padj 0.05-0.1'= '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

deg_list$MAIT<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)



##CD8 TTEM----

pseudo_ttem <- AggregateExpression(subset(cd8_t, Clusters_Final =='CD8_TTEM'), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd8_t, Clusters_Final =='CD8_TTEM'))[2]

pseudo_ttem<-pseudo_ttem$RNA

pseudo_ttem<-pseudo_ttem%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_ttem)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_ttem) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_ttem,
                            colData = samples,
                            design = ~ Group+Batch) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)

res_tidy%>%
  mutate(padj=p.adjust(pvalue, 'BH'))%>%
  mutate(labels=ifelse(padj<0.1&abs(log2FoldChange)>0.5, row, NA))%>%
  mutate(color=case_when(
    padj<0.05  ~ '1. Padj<0.05',
    padj<0.1 & padj>=0.05 ~ '2. Padj 0.05-0.1',
    TRUE ~ 'NS'
  ))%>%
  mutate(size=ifelse(padj<0.1, 'a', 'b'))%>%
  ggplot(aes(x=log2FoldChange, y=-log10(pvalue)))+
  geom_hline(yintercept = 0)+geom_vline(xintercept = 0)+
  geom_point(aes(color=color, size=size))+
  ggrepel::geom_text_repel(aes(label=labels))+scale_size_manual(values = c(1, 0.5))+
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '2. Padj 0.05-0.1'= '#46B8DACC', 'NS'='grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

deg_list$TTEM<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


#Summarize gsea results----
gsea_summary_hallmark<-deg_list$TN$Hallmark%>%as_tibble()%>%
  filter(padj<0.05)%>%
  select(Pathway=pathway, NES, pval, padj, size)%>%
  mutate(Celltype='CD8_TN')%>%
  bind_rows(deg_list$TCM$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TCM'))%>%
  bind_rows(deg_list$TEM$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TEM'))%>%
  bind_rows(deg_list$TTEM$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TTEM'))%>%
  bind_rows(deg_list$MAIT$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_MAIT'))

gsea_summary_reactome<-deg_list$TN$Reactome%>%as_tibble()%>%
  filter(padj<0.05)%>%
  select(Pathway=pathway, NES, pval, padj, size)%>%
  mutate(Celltype='CD8_TN')%>%
  bind_rows(deg_list$TCM$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TCM'))%>%
  bind_rows(deg_list$TEM$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TEM'))%>%
  bind_rows(deg_list$TTEM$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_TTEM'))%>%
  bind_rows(deg_list$MAIT$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD8_MAIT'))%>%
  filter(str_detect(Pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))

gsea_summary<- add_row(gsea_summary_reactome, gsea_summary_hallmark)%>%
  arrange(Celltype, Pathway)

#Summary
pdf('Fig 5F.pdf', width = 6, height = 2.5)
gsea_summary%>%
  mutate(Geneset=ifelse(str_detect(Pathway, 'REACTOME'), 'Reactome', 'Hallmark'))%>%
  mutate(Pathway=str_remove_all(Pathway, 'REACTOME_|HALLMARK_')%>%str_replace_all('_', ' ')%>%str_to_title())%>% 
  mutate(Celltype=factor(Celltype, levels=c('CD8_TN', 'CD8_TCM', 'CD8_TEM', 'CD8_MAIT')))%>%
  mutate(Pathway=paste0(Pathway, ' (', Celltype, ')'))%>% 
  ggplot(aes(x=NES, y=reorder(Pathway, NES))) +
  geom_point(aes(color=-log10(pval), size=size, shape=Geneset))+
  theme_minimal()+scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'F')+
  geom_vline(aes(xintercept = 0)) +
  ylab('')
dev.off()

#Save----
save.image('CD8_T.rdata')
saveRDS(cd8_t, 'Cell interaction/CD8T.rds')
 