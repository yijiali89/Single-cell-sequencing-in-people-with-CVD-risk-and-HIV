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
library(ggpubr)
library(ggsci)
source('0. Spearman_boot.R')

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

#load('NK.rdata')
#t<-readRDS('CVD_T.rds')

#Recluster----
DimPlot(nk)

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
DimPlot(subset(t, subset= Cell_Type_Main_Final %in% c('NK', 'NK and T LINGO2+')), label = TRUE, label.size = 3)+FeaturePlot(t, c('FCGR3A', 'CD8A'))

#Identification----
DotPlot(object =t, scale = T, assay = 'SCT',  cluster.idents = T, scale.by = 'size', 
        dot.min = .05, col.min = -1, col.max = 3, cols = c('grey99', 'firebrick'),
        features = c('CD3D', 'CD4', 'CD8A', 'CCR7', 'FCGR3A', 'SELL', 'CD27', 'FAS', 'CD28', 'IL7R', 'CCL5','GZMH', 'CXCR3', 'ICOS', 'MAF', 'PDCD1',
                     'S100A11', 'TIMP1', 'CST7', 'GZMK', 'GZMA', 'FOXP3', 'CTLA4', 'IL2RA', 'IFNG', 'HLA-DRB1', 'percent.mt')) /
  (FeaturePlot(t, features = c('CD4', 'CD8A', 'FCGR3A', 'NCAM1'))|DimPlot(t, label = TRUE, label.size = 3))

label<-table(t@meta.data$Cell_Type_Main_Final, t@meta.data$seurat_clusters)%>%
  as.data.frame()%>%mutate(Var2=as.character(Var2))%>%
  group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq==Freqmax)

VlnPlot(t, c('CD3D', 'CD3E', 'CD4', 'CD8A', 'FCGR3A', 'NCAM1', 'GNLY', 'NKG7'), pt.size = 0, stack = T, flip = T)|DimPlot(t, label = TRUE, label.size = 5)

DotPlot(object =t, scale = F, assay = 'SCT',  cluster.idents = T, scale.by = 'size', 
        dot.min = .05, col.min = -1, col.max = 3, cols = c('grey99', 'firebrick'),
        features = c('CD3D', 'CD3E', 'CD4', 'CD8A', 'FCGR3A', 'NCAM1', 'GNLY', 'NKG7')) |DimPlot(t, label = TRUE, label.size = 5)
nk<-subset(t, subset=seurat_clusters %in% c(label$Var2[label$Var1 %in% c('NK', 'NK and T LINGO2+')]))

#NK----
DimPlot(nk, group.by = 'Cell_Type_Main_Final')

nk<-subset(nk, Cell_Type_Main_Final %in% c('NK', 'NK and T LINGO2+'))

DimPlot(nk, group.by = 'Cell_Type_Main_Final')

FeaturePlot(nk, c('CD4', 'CD8A', 'CD3D', 'FCGR3A'))| DimPlot(nk) #The T+NK LINGO+ is a mixed group with CD4 and NK, remove those with high CD4
nk<-subset(nk, subset= umap_1<0) #Remove those CD4 high cells 

hist(unlist(FetchData(nk, c('CD8A'))), breaks = 100)
hist(unlist(FetchData(nk, c('CD8B'))), breaks = 100)
hist(unlist(FetchData(nk, c('CD4'))), breaks = 100)
hist(unlist(FetchData(nk, c('FCGR3A'))), breaks = 100)

hist(unlist(FetchData(nk, c('CD3D'))), breaks = 100)

DefaultAssay(nk)<-'SCT' #Rerun clustering for t 
nk<-RunPCA(nk)
ElbowPlot(nk)
DimHeatmap(nk, dims = 1:10, cells = 2000, balanced = TRUE)

nk<-RunUMAP(nk, dims = 1:10)

nk<-FindNeighbors(nk, dims = 1:10)

nk<-FindClusters(nk, algorithm = 2, random.seed = 111, resolution = 0.25)

DimPlot(nk, label = TRUE, label.size = 5)

DefaultAssay(nk)<-'SCT'

VlnPlot(nk, features = c('CD3E', 'FCGR3A', 'NCAM1', 'GNLY', 'CD4', 'CD8A', 'CD8B', 'FCGR3A', 'MKI67', 'FAS', 'CCR7', 'CD27', 'SELL'), pt.size = 0, stack=T, flip=T) |

FeaturePlot(nk, features = c('CD4', 'CD8A', 'CD8B', 'IL7R', 'GZMA', 'GNLY', 'GZMB', 'TCF7'))+DimPlot(nk, label = TRUE, label.size = 5)

DefaultAssay(nk)<-'SCT'

nk@assays$SCT@SCTModel.list$model1.1<-NULL
nk@assays$SCT@SCTModel.list$model1.2<-NULL
nk@assays$SCT@SCTModel.list$model1.3<-NULL

nk<-PrepSCTFindMarkers(nk)

nk.markers <- FindAllMarkers(nk, only.pos = F, logfc.threshold = .5, min.pct = .1)

write.xlsx(nk.markers, 'Fig 5A.xlsx')

nk.markers<- nk.markers%>%filter(p_val_adj<0.05)

VlnPlot(nk, c('nFeature_RNA', 'nCount_RNA'), pt.size = 0, stack = T, flip = T, split.by = 'orig.ident')

DimPlot(nk, label = TRUE, label.size = 6)

nk_features<- c('CD3G', 'KLRG1', 'FCGR3A', 'NCAM1', 'GNLY', 'FCER1G', 'SPON2', 'CD160', 'ADGRG1', 'GZMB', 'PRF1', 'IL32', 'GZMH', 'GZMK', 'IFNG', 'TNF', 'LINGO2')
                
(VlnPlot(nk, features = nk_features, 
        pt.size = 0, stack = T, flip = T) & guides(fill='none')) + DimPlot(nk, label = TRUE, label.size = 3) /
  DotPlot(object =nk, scale = T, cluster.idents = T,
          features = nk_features) + plot_layout(widths = c(1, 2)) 

#Rename cluster----
metadata_nk<-nk@meta.data
library(SingleR)
library(SingleCellExperiment)
library(celldex)
sce <- as.SingleCellExperiment(nk)
ref.set <- celldex::MonacoImmuneData()
pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.fine)
lbls.keep <- table(pred.cnts$labels)>10

nk$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')

DimPlot(nk, label=T, label.size = 4, group.by='SingleR.labels')+
  DimPlot(nk, label=T, label.size = 4)

label_nk<-table(nk@meta.data$SingleR.labels, nk@meta.data$seurat_clusters)%>%
  as.data.frame()%>%group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq>=Freqmax/2)

Idents(nk)<-'seurat_clusters'

FeaturePlot(nk, 'FCGR3A', order=T)

VlnPlot(nk, nk_features, flip=T, stack=T)

Clusters_Final_NK<-metadata_nk%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters %in% c('0')~'NK0_CD56loIL32hi_2DL3_3DL2hi', #inhibitory KIR
                                      seurat_clusters %in% c('1')~ 'NK1_GZMHhi',
                                  seurat_clusters %in% c('2')~ 'NK2_Mature_KIR2DL4', #high CD160, ADGRG1, FCER1G, SPON2
                                      seurat_clusters %in% c('3')~ 'NK3_Activated', #High IFNG, CD69
                                  seurat_clusters %in% c('4')~ 'NK4_LINGO2hi',
                                  seurat_clusters %in% c('5')~ 'NK5_GZMHhi_KIR2DS4hi',
  ))%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('NK0_CD56loIL32hi_2DL3_3DL2hi', 'NK1_GZMHhi','NK2_Mature_KIR2DL4', 
                                                        'NK3_Activated', 'NK4_LINGO2hi', 'NK5_GZMHhi_KIR2DS4hi')))%>%
  select(Barcode, Clusters_Final, seurat_clusters)%>%
  column_to_rownames('Barcode')

nk<-AddMetaData(nk, Clusters_Final_NK)

DimPlot(nk, group.by = 'seurat_clusters', label=T)
DimPlot(nk, group.by = 'Clusters_Final', label=T)

FeaturePlot(nk, c('CD3G', 'KLRG1', 'FCGR3A', 'NCAM1', 'GNLY', 'GZMK', 'FCER1G', 'IL32', 'GZMH', 'CD160', 'IFNG', 'TNF', 'LINGO2', 'KIR2DS4', 'GZMB', 'PRF1'), order=F, 
            cols = c('grey95', 'firebrick'), slot = 'data', #max.cutoff = 'q75', #min.cutoff = 'q10', 
            reduction = 'umap')  | (DimPlot(nk, label = T, group.by = 'Clusters_Final', label.size = 5) /
  DotPlot(object =nk, scale = T, assay = 'SCT',  cluster.idents = T, dot.scale = 5, scale.by='size', col.min = -1,
          cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
          features = c(nk_features)) +theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab(''))

Idents(nk)<-'Clusters_Final'

scplotter::CellDimPlot(nk, label = T, theme = "theme_blank", palette='locuszoom', dot_size=1, 
                       label_insitu = T, pt_alpha =0.5)

nk_features<- c('FCGR3A', 'NCAM1','CD3G', 'KLRG1', 'GNLY', 'FCER1G', 'SPON2', 'CD160', 'ADGRG1', 'GZMB', 'PRF1', 'IL32', 'GZMH', 'GZMK', 'IFNG', 'TNF', 'CCL3', 'LINGO2')

DotPlot(object =nk, scale = T, assay = 'SCT',  cluster.idents = T, dot.scale = 5, scale.by='size', col.min = -1,
        cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
        features = c(nk_features)) +theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab('')

VlnPlot(nk, stack = T, 
        features = cd8_features, 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) 

DefaultAssay(nk)<-'SCT'
(FeaturePlot(nk, features = c('CD8A', 'CD8B', 'CCR7', 'IL7R', 'MKI67', 'CCL5', 'GZMH'), raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, #min.cutoff = 0.1, max.cutoff = 3, 
            cols = c('grey95', 'firebrick'), slot = 'data', #split.by = 'Group', 
            reduction = 'umap') & NoAxes())| DimPlot(nk)

(FeaturePlot(nk, features = c('CCR7', 'IL7R','CCL5', 'GZMH', 'TRDC', 'MKI67'), raster = T, order = T,
             raster.dpi = c(600, 600), pt.size = 3, #min.cutoff = 0.1, max.cutoff = 3, 
             cols = c('grey95', 'firebrick'), slot = 'data',  split.by = 'Group', 
             reduction = 'umap') & NoAxes())

(FeaturePlot(subset(nk, Group=='High'), features = nk_features, raster = T, order = T,
            raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#BC3C29CC'), slot = 'data', 
            reduction = 'umap') & NoAxes()) |
  
  (FeaturePlot(subset(nk, Group=='Low'), features = nk_features, raster = T, order = T,
               raster.dpi = c(600, 600), pt.size = 3, min.cutoff = 0.1, max.cutoff = 2, cols = c('grey95', '#0072B5CC'), slot = 'data', 
               reduction = 'umap') & NoAxes())

#Compare component of cluster and disease----
metadata_nk<-nk@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(nk, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar") 

scCustomize::Proportion_Plot(seurat_object = subset(nk, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar")

cellpercent<-metadata_nk%>% #summarize percentage of each cluster_final between control and disease group
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
  facet_wrap(~Clusters_Final) +theme_bw()+stat_pwc() + scale_y_continuous(limits = c(0, 75))

cvd_metadata<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Data_CVD.xlsx')%>%select(MACSID, VL, US_Risk, USLESION:CD4NADIR), by=c('MACSID', 'US_Risk'))

#Add CVD data to metadata----
metadata_cvd_add<-metadata_nk%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

nk<-AddMetaData(nk, metadata_cvd_add)

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
DefaultAssay(nk)<-'SCT'

##Total NK ----
nkotal_for_deseq<-nk

DimPlot(nkotal_for_deseq)

pseudo_total <- AggregateExpression(nkotal_for_deseq, assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nkotal_for_deseq)[2]

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

set.seed(51)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

pdf('Fig 6E.pdf', width = 6, height = 2)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) +xlim(1, 2)
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_IL2_STAT5_SIGNALING']], stats = ranks)+
  ggtitle('HALLMARK_IL2_STAT5_SIGNALING')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

###IL2 infection gene modules----
il2_genes<-pathways.hallmark$HALLMARK_IL2_STAT5_SIGNALING

nk<-AddModuleScore(nk, features = list(il2_genes), nbin = 23, name = 'IL2_Sig')

FeaturePlot(nk, 'IL2_Sig1', split.by = 'Group')

###Oxphos gene modules----
Ox_genes<-pathways.hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION

nk<-AddModuleScore(nk, features = list(Ox_genes), nbin = 23, name = 'Ox')

###TNF gene modules----
TNF_genes<-pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB

nk<-AddModuleScore(nk, features = list(TNF_genes), nbin = 23, name = 'TNF')

###TNFA gene----
nk<-AddModuleScore(nk, features = 'TNF', nbin = 23, name = 'TNFA')

###IFNG gene modules----
IFNG_genes<-pathways.hallmark$HALLMARK_INTERFERON_GAMMA_RESPONSE

nk<-AddModuleScore(nk, features = list(IFNG_genes), nbin = 23, name = 'IFNG_Sig')

VlnPlot(nk, c('IL2_Sig1', 'Ox1', 'TNF1', 'TNFA1'), split.by = 'Group', pt.size = 0, stack=T, flip=T) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))+
  scale_y_continuous(expand = c(0.25,0))

nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, IL2_Sig1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(IL2_Sig1=mean(IL2_Sig1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=IL2_Sig1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, Ox1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(Ox1=mean(Ox1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=Ox1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNF1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(TNF1=mean(TNF1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=TNF1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNFA1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(TNFA1=mean(TNFA1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=TNFA1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  plot_layout(guides = "collect") & theme(legend.position = 'bottom')

##NK0 ----
nk0_for_deseq<-subset(nk, subset=Clusters_Final=='NK0_CD56loIL32hi_2DL3_3DL2hi')

DimPlot(nk0_for_deseq)

pseudo_nk0 <- AggregateExpression(nk0_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk0_for_deseq)[2]

pseudo_nk0<-pseudo_nk0$RNA

pseudo_nk0<-pseudo_nk0%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk0)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk0) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk0,
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

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)
 
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##NK1 ----
nk1_for_deseq<-subset(nk, subset=Clusters_Final=='NK1_GZMHhi')

DimPlot(nk1_for_deseq)

pseudo_nk1 <- AggregateExpression(nk1_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk1_for_deseq)[2]

pseudo_nk1<-pseudo_nk1$RNA

pseudo_nk1<-pseudo_nk1%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk1)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk1) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk1,
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

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##NK2 ----
nk2_for_deseq<-subset(nk, subset=Clusters_Final=='NK2_Mature_KIR2DL4')

DimPlot(nk2_for_deseq)

pseudo_nk2 <- AggregateExpression(nk2_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk2_for_deseq)[2]

pseudo_nk2<-pseudo_nk2$RNA

pseudo_nk2<-pseudo_nk2%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk2)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk2) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk2,
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
 
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

pdf('Fig 6G.pdf', width = 6, height = 2)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##NK3 ----
nk3_for_deseq<-subset(nk, subset=Clusters_Final=='NK3_Activated')

DimPlot(nk3_for_deseq)

pseudo_nk3 <- AggregateExpression(nk3_for_deseq, assays = "RNA", 
                                    slot = 'counts', return.seurat = F, 
                                    group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk3_for_deseq)[2]

pseudo_nk3<-pseudo_nk3$RNA

pseudo_nk3<-pseudo_nk3%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk3)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk3) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk3,
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

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

pdf('Fig 6G.pdf', width = 6, height = 2)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##NK4 ----
nk4_for_deseq<-subset(nk, subset=Clusters_Final=='NK4_LINGO2hi')

DimPlot(nk4_for_deseq)

pseudo_nk4 <- AggregateExpression(nk4_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk4_for_deseq)[2]

pseudo_nk4<-pseudo_nk4$RNA

pseudo_nk4<-pseudo_nk4%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk4)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk4) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk4,
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
 
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)
 
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 



##NK5 ----
nk5_for_deseq<-subset(nk, subset=Clusters_Final=='NK5_GZMHhi_KIR2DS4hi')

DimPlot(nk5_for_deseq)

pseudo_nk5 <- AggregateExpression(nk5_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(nk5_for_deseq)[2]

pseudo_nk5<-pseudo_nk5$RNA

pseudo_nk5<-pseudo_nk5%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_nk5)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_nk5) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_nk5,
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

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

pdf('Fig 6G.pdf', width = 6, height = 2)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0)) 
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB, stats = ranks)+
  ggtitle('HALLMARK_TNFA_SIGNALING_VIA_NFKB')

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 


#Modules and CVD correlations----
nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNF1, IMT:SEG_SCORE, MMP2:CD4NADIR)%>%
  #filter(Clusters_Final %in% c('pDC', 'cDC')==FALSE)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(TNF1=mean(TNF1), across(IMT:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  pivot_longer(cols=IMT:CD4NADIR, names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  filter(Var %in% c('Intact', 'CRP', 'TNFARI', 'OSTEOP'), Clusters_Final=='NK3_Activated')%>%
  mutate(Var=factor(Var, levels = c('Intact', 'CRP', 'TNFARI', 'OSTEOP')))%>%
  ggplot(aes(x=TNF1, y=log10(Value)))+
  geom_smooth(method = 'lm', se=F, linetype='dashed', alpha=0.1, color='grey')+
  geom_point(aes(color=Group))+
  scale_color_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  theme_bw()+facet_wrap(~Var, scales = 'free')+
  stat_cor(method = 'spearman', label.y.npc = 'bottom', label.x.npc = 'middle')

nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNF1, IMT:SEG_SCORE, MMP2:CD4NADIR)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(TNF1=mean(TNF1), across(IMT:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  pivot_longer(cols=IMT:CD4NADIR, names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  filter(Var %in% c('Intact', 'CRP', 'TNFARI', 'OSTEOP'), Clusters_Final=='NK3_Activated')%>%
  mutate(Var=factor(Var, levels = c('Intact', 'CRP', 'TNFARI', 'OSTEOP')))%>%
  ggplot(aes(x=TNF1, y=log10(Value)))+
  geom_point(aes(color=Group))+
  geom_smooth(method = 'lm', se=F)+
  scale_color_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  theme_bw()+facet_wrap(~Var, scales = 'free')+
  stat_cor(method = 'spearman', label.y.npc = 'bottom', label.x.npc = 'middle')

#Final results CI bootstrap
set.seed(111)
cor_boot <- nk@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNF1, Intact, CRP, TNFARI, OSTEOP)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(TNF1=mean(TNF1), across(Intact:OSTEOP))%>%
  ungroup()%>%unique() %>%
  filter(Clusters_Final=='NK3_Activated') %>%
  pivot_longer(cols=Intact:OSTEOP, names_to = 'Item', values_to ='Value')%>%
  filter(is.na(Value)==FALSE)%>%
  group_by(Item) %>%
  group_modify(~ spearman_boot_ci(data = .x, x = 'TNF1', y = 'Value', R = 2000)) %>%
  ungroup()

cor_boot

#Save----
save.image('NK.rdata')
 