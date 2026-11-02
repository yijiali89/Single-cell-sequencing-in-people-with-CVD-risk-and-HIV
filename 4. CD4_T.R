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

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

#load('CD4_T.rdata')
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

cd4_t<-subset(t, subset=seurat_clusters %in% c(label$Var2[label$Var1=='CD4+ T'], '20', '26'))

#CD4----
DimPlot(cd4_t, group.by = 'Cell_Type_Main_Final')

cd4_t<-subset(cd4_t, Cell_Type_Main_Final %in% c('CD4+ T', 'T', 'T Prolif'))

FeaturePlot(cd4_t, c('CD4', 'CD8A'))| DimPlot(cd4_t)

hist(unlist(FetchData(cd4_t, c('CD8A'))), breaks = 100)
hist(unlist(FetchData(cd4_t, c('CD8B'))), breaks = 100)
hist(unlist(FetchData(cd4_t, c('CD4'))), breaks = 100)
hist(unlist(FetchData(cd4_t, c('FCGR3A'))), breaks = 100)

cd4_t<-subset(cd4_t, subset= CD8A<=CD4 & CD8B<=CD4 & FCGR3A<=CD4) #Filter those CD8A and CD8B and FCGR3A high cells
hist(unlist(FetchData(cd4_t, c('CD8A'))), breaks = 1000)
hist(unlist(FetchData(cd4_t, c('FCGR3A'))), breaks = 1000)

DefaultAssay(cd4_t)<-'SCT' #Rerun clustering for t 
cd4_t<-RunPCA(cd4_t)
ElbowPlot(cd4_t)
DimHeatmap(cd4_t, dims = 1:10, cells = 1000, balanced = TRUE)

cd4_t<-RunUMAP(cd4_t, dims = 1:10)

cd4_t<-FindNeighbors(cd4_t, dims = 1:10)

cd4_t<-FindClusters(cd4_t, algorithm = 2, random.seed = 121, resolution = 1)

DimPlot(cd4_t, label = TRUE, label.size = 5)

DefaultAssay(cd4_t)<-'SCT'

VlnPlot(cd4_t, features = c('CD3D', 'FCGR3A', 'CD4', 'CD8A', 'MKI67', 'FAS', 'CCR7'), pt.size = 0, stack=T, flip=T)

FeaturePlot(cd4_t, features = c('FCGR3A', 'CD4', 'CD8A', 'IL7R', 'GZMA', 'GNLY', 'MKI67'))+DimPlot(cd4_t, label = TRUE, label.size = 5)

DefaultAssay(cd4_t)<-'SCT'

cd4_t@assays$SCT@SCTModel.list$model1.1<-NULL
cd4_t@assays$SCT@SCTModel.list$model1.2<-NULL
cd4_t@assays$SCT@SCTModel.list$model1.3<-NULL

cd4_t<-PrepSCTFindMarkers(cd4_t)

cd4_t.markers <- FindAllMarkers(cd4_t, only.pos = F, logfc.threshold = .5, min.pct = .1)

write.xlsx(cd4_t.markers, 'Fig 4A.xlsx')

cd4_t.markers<- cd4_t.markers%>%filter(p_val_adj<0.05)

VlnPlot(cd4_t, c('nFeature_RNA', 'nCount_RNA'), pt.size = 0, stack = T, flip = T, split.by = 'orig.ident')

VlnPlot(cd4_t, features = c('CD4', 'CCR7', 'FAS', 'FOXP3', 'IL7R'), 
        pt.size = 0, stack = T, flip = T) + DimPlot(cd4_t, label = TRUE, label.size = 3) /
  DotPlot(object =cd4_t, scale = F, 
          features = c('CCR7', 'FAS', 'IL7R'))

cd4_features<- c('IL7R', 'SELL', 'CCR7', 'LEF1', 'TCF7', 'FHIT', 'BACH2', 'AQP3', 'FAS', 'KLRB1',  'GZMH', 'GZMB', 'GZMK', 'GNLY', 'CCL5', 'GATA3','ICOS', 'CXCR5', 'PDCD1', 'FOXP3', 'IL2RA', 'HLA-DRB1', 'HLA-DRA', 'MKI67')

DotPlot(object =cd4_t, scale = T, assay = 'SCT',  cluster.idents = T, dot.scale = 8, scale.by='size', col.min = -1,
        cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
        features = cd4_features)/(
  DimPlot(cd4_t, label = T, label.size = 5, label.color = 'black') |
  VlnPlot(object = cd4_t, features = cd4_features, pt.size = 0, stack=T, flip=T))

DotPlot(object =subset(cd4_t, subset= seurat_clusters %in% c('6', '3', '5', '0', '14', '7', '8', '9', '12')==FALSE), 
        scale = T, assay = 'RNA',  cluster.idents = T, scale.by = 'size',
        dot.min = .05, cols = c('grey99', 'firebrick'),
        features = cd4_features) / (FeaturePlot(subset(cd4_t, subset= seurat_clusters %in% c('6', '3', '5', '0', '14', '7', '8', '9', '12')==FALSE), c('CD4', 'FAS', 'CCR7'))+ DimPlot(cd4_t, label = T))

FeaturePlot(cd4_t, c('GZMH', 'GZMK', 'GZMB', 'GNLY', 'CD3D', 'HLA-DRB1'), order=T)
FeaturePlot(cd4_t, c('BACH2'), order=T)

#Rename cluster----
metadata_cd4t<-cd4_t@meta.data
library(SingleR)
library(SingleCellExperiment)
library(celldex)
sce <- as.SingleCellExperiment(cd4_t)
ref.set <- celldex::MonacoImmuneData()
pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.fine)
lbls.keep <- table(pred.cnts$labels)>10

cd4_t$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')

DimPlot(cd4_t, label=T, label.size = 4, group.by='SingleR.labels')+
  DimPlot(cd4_t, label=T, label.size = 4)

label_cd4<-table(cd4_t@meta.data$SingleR.labels, cd4_t@meta.data$seurat_clusters)%>%
  as.data.frame()%>%group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq==Freqmax)

Idents(cd4_t)<-'seurat_clusters'
marker_0v10<-FindMarkers(cd4_t, ident.1 = '0', ident.2 = '10')%>%filter(p_val_adj<0.05, avg_log2FC>0.25)
marker_6v1<-FindMarkers(cd4_t, ident.1 = '6', ident.2 = '1')%>%filter(p_val_adj<0.05, abs(avg_log2FC)>0.25)

FeaturePlot(cd4_t, 'GZMK', order=T)
VlnPlot(cd4_t, cd4_features, flip=T, stack=T)

Clusters_Final_CD4<-metadata_cd4t%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters %in% c('2', '3', '5', '9') ~ 'CD4_TN', #CD4 TN, high FHIT
                                  seurat_clusters =='4'~ 'CD4_Treg', #high FOXP3 
                                  seurat_clusters %in% c('11', '13', '12', '14') ~ 'CD4_TCM', #CCR7+, FAS+, CD28+, IL7R+
                                  seurat_clusters %in% c('1', '8', '6') ~ 'CD4_TEM_TFH', #CXCR5+, ICOS+
                                  seurat_clusters %in% c('0', '10', '7') ~ 'CD4_TEM_Th1/17',
                                  seurat_clusters %in% c('15')~'CD4_TEM_Th2', #GATA3 high
                                  seurat_clusters %in% c('16') ~ 'CD4_Prolif'
  ))%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('CD4_TN', 'CD4_TCM','CD4_TEM_TFH', 'CD4_TEM_Th1/17', 'CD4_TEM_Th2',
                                                        'CD4_Treg', 'CD4_Prolif')))%>%
  select(Barcode, Clusters_Final, seurat_clusters)%>%
  column_to_rownames('Barcode')

cd4_t<-AddMetaData(cd4_t, Clusters_Final_CD4)
DimPlot(cd4_t, group.by = 'Clusters_Final', label=T)

Idents(cd4_t)<-'Clusters_Final'

scplotter::CellDimPlot(cd4_t, label = T, theme = "theme_blank", palette='locuszoom', dot_size=2, 
                       label_insitu = T, pt_alpha =0.5)

scplotter::CellDimPlot(subset(cd4_t, Clusters_Final%in% c('CD4_TEM_TFH', 'CD4_TEM_Th1/17')), label = T, theme = "theme_blank", palette='locuszoom', dot_size=2, 
                       label_insitu = T, pt_alpha =0.5, group_by='seurat_clusters')

DotPlot(object =cd4_t, scale = T, assay = 'SCT',  cluster.idents = F, dot.scale = 5, scale.by='size', col.min = -1,
        cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100, 
        features = cd4_features) +theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab('')

VlnPlot(cd4_t, stack = T, features = c('CD4', 'CD8A', 'CCR7', 'TCF7', 'SELL', 'IL7R', 'CD27', 'FAS', 'AQP3', 'KLRB1', 'S100A4', 'FOXP3', 'IL2RA'), 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) 

DefaultAssay(cd4_t)<-'SCT'
FeaturePlot(cd4_t, features = c('CD4', 'CCR7', 'FHIT', 'FAS', 'PDCD1', 'GZMA', 'ICOS', 'CXCR5', 'FOXP3'), raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 3, cols = c('grey95', 'firebrick'), slot = 'data', 
            reduction = 'umap', split.by = 'Group') & NoAxes()

DotPlot(cd4_t, features = c('LEF1', 'CCR7', 'FHIT', 'TCF7', 'SELL', 'FAS', 'GZMA', 'FOXP3', 'IL2RA', 'CTLA4', 'PDCD1'), scale =T, cluster.idents = T,
        dot.min = .05, dot.scale = 6, scale.by = 'size', col.max = 1, col.min = -1)+ theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(color = 'grey90'), panel.grid.minor = element_line(color = 'grey90'))

(FeaturePlot(subset(cd4_t, Group=='High'), features = c('CD4', 'CCR7', 'FAS', 'GATA3', 'FOXP3', 'MKI67'), raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#BC3C29CC'), slot = 'data', 
            reduction = 'umap') & NoAxes()) |
  
  (FeaturePlot(subset(cd4_t, Group=='Low'), features = c('CD4', 'CCR7', 'FAS', 'GATA3', 'FOXP3', 'MKI67'), raster = T, order = T,
               raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#0072B5CC'), slot = 'data', 
               reduction = 'umap') & NoAxes())

#Compare component of cluster and disease----
metadata_cd4<-cd4_t@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(cd4_t, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar") 

scCustomize::Proportion_Plot(seurat_object = subset(cd4_t, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar")

cellpercent<-metadata_cd4%>% #summarize percentage of each cluster_final between control and disease group
  group_by(Group, Sample_Name, Clusters_Final, .drop=FALSE)%>%
  summarise(Count=n())%>%
  mutate(Percentage=Count/sum(Count)*100)%>%
  ungroup()

#Save a file for it
cellpercent%>%
  pivot_wider(id_cols=Sample_Name, names_from = c(Clusters_Final, Group), values_from = Percentage)%>%
  write.xlsx('tmp.xlsx')

cellpercent%>%
  ggplot(aes(x=Group, y=Percentage, color=Group))+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  geom_jitter(aes(color=Group), width=0.2)+
  scale_color_manual(values = c("#BC3C29CC", "#0072B5CC"))+
  facet_wrap(~Clusters_Final) +theme_bw()+stat_pwc()+scale_y_continuous(limits = c(0, 65))

cvd_metadata<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Data_CVD.xlsx')%>%select(MACSID, VL, US_Risk, USLESION:CD4NADIR), by=c('MACSID', 'US_Risk'))

#Add CVD data to metadata----
metadata_cvd_add<-metadata_cd4%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

cd4_t<-AddMetaData(cd4_t, metadata_cvd_add)

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
DefaultAssay(cd4_t)<-'SCT'
##Total CD4 ----
pseudo_total <- AggregateExpression(cd4_t, assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(cd4_t)[2]

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
                            design = ~ Group+Batch+VL_Blip) #Control batch

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

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_MYC_TARGETS_V1']], stats = ranks)+
  ggtitle('HALLMARK_MYC_TARGETS_V1')

VlnPlot(cd4_t, features = c('NFKBIA'), group.by = 'Group', pt.size = 0)+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20") +
  stat_summary(fun.y=median, geom='point', size=5, shape=20) +stat_compare_means()

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

cd4_t<-AddModuleScore(cd4_t, features = list(myc_genes), name = 'MYC')

###Viral infection gene modules----
vi_genes<-pathway.reactome$REACTOME_VIRAL_INFECTION_PATHWAYS

cd4_t<-AddModuleScore(cd4_t, features = list(vi_genes), name = 'VI')

FeaturePlot(cd4_t, 'MYC1', split.by = 'Group')

###IFN gene modules----
IFNA_genes<-pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE
IFNG_genes<-pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE

cd4_t<-AddModuleScore(cd4_t, features = list(IFNA_genes), name = 'IFNA_Mod')
cd4_t<-AddModuleScore(cd4_t, features = list(IFNG_genes), name = 'IFNG_Mod')

###Estrogen gene modules----
ES_early_genes<-pathways.hallmark$HALLMARK_ESTROGEN_RESPONSE_EARLY
cd4_t<-AddModuleScore(cd4_t, features = list(ES_early_genes), name = 'ES_Early')


VlnPlot(cd4_t, c('IFNA_Mod1', 'ES_Early1'), split.by = 'Group', pt.size = 0, stack=T, flip=T) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))+
  scale_y_continuous(expand = c(0.25,0))

(VlnPlot(cd4_t, 'IL12_1', split.by = 'Group', pt.size = 0) +
    stat_compare_means(label = 'p.format')+
    scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
    stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
    stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9)))/

(VlnPlot(cd4_t, 'VI1', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))) /

(VlnPlot(cd4_t, 'HIV1', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+theme(axis.text.x = element_text(angle = 0))+xlab('')+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9)))

VlnPlot(cd4_t, 'BCL3', split.by = 'Group', pt.size = 0) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))
 

##CD4 TN----
pseudo_tn <- AggregateExpression(subset(cd4_t, Clusters_Final=='CD4_TN'), assays = "RNA", 
                                    slot = 'counts', return.seurat = F, 
                                    group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd4_t, Clusters_Final=='CD4_TN'))[2]

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
                            design = ~ Group+Batch+VL_Blip) #Control batch

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

deg_list<-list()

deg_list$TN<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


##CD4 TCM----
pseudo_tcm <- AggregateExpression(subset(cd4_t, Clusters_Final=='CD4_TCM'), assays = "RNA", 
                                 slot = 'counts', return.seurat = F, 
                                 group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd4_t, Clusters_Final=='CD4_TCM'))[2]

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
                            design = ~ Group+Batch+VL_Blip) #Control batch

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


##CD4 TEM----
DimPlot(subset(cd4_t, Clusters_Final %in% c('CD4_TEM_TFH','CD4_TEM_Th1/17', 
                                            'CD4_TEM_Th2')))

pseudo_tem <- AggregateExpression(subset(cd4_t, Clusters_Final %in% c('CD4_TEM_TFH','CD4_TEM_Th1/17', 
                                                                      'CD4_TEM_Th2')), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd4_t, Clusters_Final %in% c('CD4_TEM_TFH','CD4_TEM_Th1/17', 
                                                'CD4_TEM_Th2')))[2]

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
                            design = ~ Group+Batch+VL_Blip) #Control batch

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

##CD4 Treg----
pseudo_treg <- AggregateExpression(subset(cd4_t, Clusters_Final=='CD4_Treg'), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(cd4_t, Clusters_Final=='CD4_Treg'))[2]

pseudo_treg<-pseudo_treg$RNA

pseudo_treg<-pseudo_treg%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_treg)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_treg) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_treg,
                            colData = samples,
                            design = ~ Group+Batch+VL_Blip) #Control batch

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

deg_list$Treg<-list(DEG=res_tidy, Hallmark=fgseaRes, Reactome=fgseaRes_reactome)


#Summarize gsea results----
gsea_summary_hallmark<-deg_list$TN$Hallmark%>%as_tibble()%>%
  filter(padj<0.05)%>%
  select(Pathway=pathway, NES, pval, padj, size)%>%
  mutate(Celltype='CD4_TN')%>%
  bind_rows(deg_list$TCM$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_TCM'))%>%
  bind_rows(deg_list$TEM$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_TEM'))%>%
  bind_rows(deg_list$Treg$Hallmark%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_Treg'))

gsea_summary_reactome<-deg_list$TN$Reactome%>%as_tibble()%>%
  filter(padj<0.05)%>%
  select(Pathway=pathway, NES, pval, padj, size)%>%
  mutate(Celltype='CD4_TN')%>%
  bind_rows(deg_list$TCM$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_TCM'))%>%
  bind_rows(deg_list$TEM$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_TEM'))%>%
  bind_rows(deg_list$Treg$Reactome%>%as_tibble()%>%
              filter(padj<0.05)%>%
              select(Pathway=pathway, NES, pval, padj, size)%>%
              mutate(Celltype='CD4_Treg'))%>%
  filter(str_detect(Pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))

gsea_summary<- add_row(gsea_summary_reactome, gsea_summary_hallmark)%>%
  arrange(Celltype, Pathway)

#Plot a polar plot showing positive NES pathway from each cell type; each cell type occupy a quarter of the round polar plot-
gsea_summary%>%
  mutate(Geneset=ifelse(str_detect(Pathway, 'REACTOME'), 'Reactome', 'Hallmark'))%>%
  mutate(Pathway=str_remove_all(Pathway, 'REACTOME_|HALLMARK_')%>%str_replace_all('_', ' ')%>%str_to_title())%>% 
  mutate(Celltype=factor(Celltype, levels=c('CD4_TN', 'CD4_TCM', 'CD4_TEM', 'CD4_Treg')))%>%
  ggplot(aes(x=NES, y=-log10(pval)))+
  geom_point(aes(color=Celltype, size=size, shape=Geneset))+
  theme_minimal()+
  theme(axis.text.x = element_text(size=8)) +
  scale_color_manual(values = c('CD4_TN'="#d43f3a", 'CD4_TCM'="#ea9336", 'CD4_TEM'="#4fb8a3", 'CD4_Treg'="#9a45b8"))+
  ggrepel::geom_text_repel(aes(label=ifelse(padj<0.01& abs(NES)>1.5, Pathway, '')), size=3)+ylim(0, 10)

#Save----
save.image('CD4_T.rdata')

#ISR signal----
isr<- read.xlsx('ISR.xlsx')%>%
  as.vector()%>%unlist()
DimPlot(cd4_t)

cd4_t<- AddModuleScore(cd4_t, features = list(isr), name = 'ISR')
FeaturePlot(cd4_t, features = 'ISR1', min.cutoff = 'q50', max.cutoff = 'q75', split.by = 'Group')
