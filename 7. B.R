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

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

#load('B.rdata')
#b<-readRDS('CVD_B_Plasma.rds')

#Recluster----
DimPlot(b)

DefaultAssay(b)<-'SCT' #Rerun clustering for B
b<-RunPCA(b)
ElbowPlot(b, ndims = 30)
DimHeatmap(b, dims = 1:20, cells = 500, balanced = TRUE)

b<-RunUMAP(b, dims = 1:20)

b<-FindNeighbors(b, dims = 1:20)
b<-FindClusters(b, algorithm = 2, random.seed = 121, resolution = 0.25)

Idents(b)<-'seurat_clusters'
DefaultAssay(b)<-'SCT'
DimPlot(b, label = TRUE, label.size = 3)

DimPlot(b, group.by = 'Cell_Type_Main_Final')

#Identification----
label<-table(b@meta.data$Cell_Type_Main_Final, b@meta.data$seurat_clusters)%>%
  as.data.frame()%>%mutate(Var2=as.character(Var2))%>%
  group_by(Var2)%>%mutate(Freqmax=max(Freq))%>%ungroup()%>%filter(Freq==Freqmax)
 
FeaturePlot(b, c('MS4A1', 'CD38'))| DimPlot(b)  
 
b@assays$SCT@SCTModel.list$model1.1<-NULL
b@assays$SCT@SCTModel.list$model1.2<-NULL
b@assays$SCT@SCTModel.list$model1.3<-NULL

b<-PrepSCTFindMarkers(b)

b.markers <- FindAllMarkers(b, only.pos = F, logfc.threshold = .5, min.pct = .1)

write.xlsx(b.markers, 'Fig 7A.xlsx')

b.markers<- b.markers%>%filter(p_val_adj<0.05)

VlnPlot(b, c('nFeature_RNA', 'nCount_RNA'), pt.size = 0, stack = T, flip = T, split.by = 'orig.ident')

DimPlot(b, label = TRUE, label.size = 6)

b_features<-c('MS4A1', 'CD19', 'CD38', 'CD79A', 'IGHD', 'IGHM', 'IGHG1', 'IGHA1',  'IL4R', 'TCL1A', 'IGKC', 'IGLC2', 'BCL6')
#TCL1A, IGHD and IL4R are markers for naive B, which is consistent with the annotation of cluster 0, also IGKC+. 
#TCL1A, IGHD and IL4R are markers for naive B, which is consistent with the annotation of cluster 0, also IGLC+. 
#Cluster 2 is annotated as memory B AIM2+
#Cluster 3 is annotated as plasma cell because of high CD38 and low MS4A1.
#Cluster 4 has high CD3E and will be removed because of potential T cell contamination.

b<-subset(b, subset = seurat_clusters != '4')

(VlnPlot(b, features = b_features, 
        pt.size = 0, stack = T, flip = T) & guides(fill='none')) + DimPlot(b, label = TRUE, label.size = 3) /
  DotPlot(object =b, scale = T, cluster.idents = T,
          features = b_features) + plot_layout(widths = c(1, 2)) 

marker0v1<-FindMarkers(b, ident.1 = '0', ident.2 = '1', logfc.threshold = 0.5, min.pct = 0.1)%>%filter(p_val_adj<0.05)
marker2v3<-FindMarkers(b, ident.1 = '2', ident.2 = '3', logfc.threshold = 0.5, min.pct = 0.1)%>%filter(p_val_adj<0.05)
marker0v3<-FindMarkers(b, ident.1 = '0', ident.2 = '3', logfc.threshold = 0.5, min.pct = 0.1)%>%filter(p_val_adj<0.05)

#Rename cluster----
metadata_b<-b@meta.data

Clusters_Final_b<-metadata_b%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters %in% c('0')~'B0_Naive_IGKC+', #inhibitory KIR
                                      seurat_clusters %in% c('1')~ 'B1_Naive_IGLC+',
                                  seurat_clusters %in% c('2')~ 'B3_Memory', #high CD160, ADGRG1, FCER1G, SPON2
                                      seurat_clusters %in% c('3')~ 'B2_Intermediate', #High IFNG, CD69
                                  seurat_clusters %in% c('5')~ 'Plasma cells',
  ))%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('B0_Naive_IGKC+', 'B1_Naive_IGLC+', 'B2_Intermediate', 
                                                        'B3_Memory', 'Plasma cells')))%>%
  select(Barcode, Clusters_Final, seurat_clusters)%>%
  column_to_rownames('Barcode')

b<-AddMetaData(b, Clusters_Final_b)

DimPlot(b, group.by = 'seurat_clusters', label=T)
DimPlot(b, group.by = 'Clusters_Final', label=T)

b_features<-c('MS4A1', 'CD19', 'CD38', 'CD27', 'IGHD', 'IGHM', 'IGHG1', 'IGHA1',  
              'IL4R', 'TCL1A', 'IGKC', 'IGLC2', 'BCL6', 'VPREB1')

FeaturePlot(b, b_features, order=F, 
            cols = c('grey95', 'firebrick'), slot = 'data', #max.cutoff = 'q75', #min.cutoff = 'q10', 
            reduction = 'umap')  | (DimPlot(b, label = T, group.by = 'Clusters_Final', label.size = 5) /
  DotPlot(object =b, scale = T, assay = 'SCT',  cluster.idents = T, dot.scale = 5, scale.by='size', col.min = -1,
          cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
          features = c(b_features)) +theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab(''))

Idents(b)<-'Clusters_Final'

scplotter::CellDimPlot(b, label = T, theme = "theme_blank", palette='seurat', dot_size=1, 
                       label_insitu = T, pt_alpha =0.5)

DotPlot(object =b, scale = T, assay = 'SCT',  cluster.idents = T, dot.scale = 5, scale.by='size', col.min = -1,
        cols = c('grey99', 'firebrick4'), scale.min = 1, scale.max = 100,
        features = c(b_features)) +theme_bw()+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+xlab('')+ylab('')

VlnPlot(b, stack = T, 
        features = b_features, 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) 

b_features<-c('MS4A1', 'CD19', 'CD38', 'IL4R', 'IGHD', 'IGHM', 'IGHG1', 'IGHA1',  
              'IGKC', 'IGLC2')

FeaturePlot(b, features= b_features, 
           min.cutoff = 'q5', max.cutoff = 'q95', cols =c('grey90', "#BC3C29CC"), raster=T, pt.size = 3, raster.dpi = c(600, 600),
            reduction = "umap", ncol = 5)

((FeaturePlot(subset(b, subset= Group=='High'), features = b_features, 
              keep.scale = 'all', min.cutoff = 0, max.cutoff = 4.5, cols =c('grey90', "#BC3C29CC"), raster=T, pt.size = 3, raster.dpi = c(600, 600),
              reduction = "umap", ncol = 4)+plot_annotation('US high-risk')) /
    
    (FeaturePlot(subset(b, subset= Group=='Low'), features = b_features, 
                 keep.scale = 'all',  min.cutoff = 0, max.cutoff = 4.5, cols =c('grey90',"#0072B5CC"), raster=T, pt.size = 3, raster.dpi = c(600, 600),
                 reduction = "umap", ncol = 4) +plot_annotation('US low-risk'))) & NoAxes()

FeaturePlot(b, 'ATP6V0E1', min.cutoff = 'q10', max.cutoff = 'q90', split.by = 'Group')

#Compare component of cluster and disease----
metadata_b<-b@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(b, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar") 

scCustomize::Proportion_Plot(seurat_object = subset(b, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar")

cellpercent<-metadata_b%>% #summarize percentage of each cluster_final between control and disease group
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
metadata_cvd_add<-metadata_b%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

b<-AddMetaData(b, metadata_cvd_add)

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
DefaultAssay(b)<-'SCT'

##Total b ----
btotal_for_deseq<-subset(b, subset= Clusters_Final != 'Plasma cells')

DimPlot(btotal_for_deseq)

pseudo_total <- AggregateExpression(btotal_for_deseq, assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(btotal_for_deseq)[2]

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

pdf('Fig 6E.pdf', width = 6, height = 2)
fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 
dev.off()

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_OXIDATIVE_PHOSPHORYLATION']], stats = ranks)+
  ggtitle('HALLMARK_OXIDATIVE_PHOSPHORYLATION')

fgsea_all<-fgseaRes%>%as_tibble()%>%filter(padj<0.05)

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

###Oxphos gene modules----
Oxphos<-pathways.hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION

b$Oxphos1<-NULL

b<-AddModuleScore(b, features = list(Oxphos), name = 'Oxphos')

###TNF gene modules----
TNF_genes<-pathways.hallmark$HALLMARK_TNFA_SIGNALING_VIA_NFKB

b$TNF1<-NULL

b<-AddModuleScore(b, features = list(TNF_genes), name = 'TNF')

###Glycolysis gene modules----
Glyco_genes<-pathways.hallmark$HALLMARK_GLYCOLYSIS

b<-AddModuleScore(b, features = list(Glyco_genes), name = 'Glyco')

VlnPlot(b, c('Oxphos1', 'TNF1'), split.by = 'Group', pt.size = 0, stack=T, flip=T) +
  stat_compare_means(label = 'p.format')+
  scale_fill_manual(values=c( "#BC3C29CC", "#0072B5CC"))+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(.9))+
  scale_y_continuous(expand = c(0.25,0))

##B_naive ----
bnaive_for_deseq<-subset(b, subset=Clusters_Final %in% c('B0_Naive_IGKC+', 'B1_Naive_IGLC+'))

DimPlot(bnaive_for_deseq)

pseudo_bnaive0 <- AggregateExpression(bnaive_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(bnaive_for_deseq)[2]

pseudo_bnaive0<-pseudo_bnaive0$RNA

pseudo_bnaive0<-pseudo_bnaive0%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_bnaive0)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_bnaive0) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_bnaive0,
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

fgsea_naive<-fgseaRes%>%as_tibble()%>%filter(padj<0.05)

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##B_inter ----
binter_for_deseq<-subset(b, subset=Clusters_Final=='B2_Intermediate')

DimPlot(binter_for_deseq)

pseudo_binter <- AggregateExpression(binter_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(binter_for_deseq)[2]

pseudo_binter<-pseudo_binter$RNA

pseudo_binter<-pseudo_binter%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_binter)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_binter) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_binter,
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

fgsea_binter<-fgseaRes%>%as_tibble()%>%filter(padj<0.05)

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  filter(pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(10, 50, 100, 500, 1000), limits = c(10, 800))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(0, 3)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) 

##B_memory ----
bmemory_for_deseq<-subset(b, subset=Clusters_Final=='B3_Memory')

DimPlot(bmemory_for_deseq)

pseudo_bmem <- AggregateExpression(bmemory_for_deseq, assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(bmemory_for_deseq)[2]

pseudo_bmem<-pseudo_bmem$RNA

pseudo_bmem<-pseudo_bmem%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_bmem)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_bmem) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_bmem,
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

##Plasma cells ----
plasma_for_deseq<-subset(b, subset=Clusters_Final=='Plasma cells')

DimPlot(plasma_for_deseq)

pseudo_plasma <- AggregateExpression(plasma_for_deseq, assays = "RNA", 
                                    slot = 'counts', return.seurat = F, 
                                    group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(plasma_for_deseq)[2]

pseudo_plasma<-pseudo_plasma$RNA

pseudo_plasma<-pseudo_plasma%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_plasma)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL)%>%mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

rownames(samples)<-colnames(pseudo_plasma) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_plasma,
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

fgsea_plasma<-fgseaRes%>%as_tibble()%>%filter(padj<0.05)

#Combine the three B cell subsets----
fgsea_bcell<-bind_rows(fgsea_all%>%mutate(Subtype='All B cells'),
                        fgsea_naive%>%mutate(Subtype='Naive B cells'), 
                       fgsea_binter%>%mutate(Subtype='Intermediate B cells'), 
                       fgsea_plasma%>%mutate(Subtype='Plasma cells')) %>%
  filter(str_detect(pathway, pattern = 'TNF|MYC|OXI'))%>% 
  mutate(Subtype=factor(Subtype, levels = c('All B cells', 'Naive B cells', 'Intermediate B cells', 'Plasma cells')))
  
fgsea_bcell%>%
  mutate(pathway=str_remove_all(pathway, 'HALLMARK_')%>%
  str_replace_all('_', ' ')%>%
  str_to_title(.))%>%
  mutate(pathway=paste0(pathway, '(', Subtype, ')'))%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size), position = position_dodge(0.8))+
  theme_minimal()+scale_size_continuous(breaks = c(1, 50, 100, 150, 200), limits = c(10, 250))+
  scale_color_viridis_c(begin = 0.3, end=0.9, direction = -1, option = 'B')+geom_vline(aes(xintercept = 0))+ylab('')

#Modules and CVD correlations----
pal_b<-c('B0_Naive_IGKC+'='#fecf8c', 'B1_Naive_IGLC+'='#87ab00', 'B2_Intermediate'='#00c0af', 
         'B3_Memory'='#948eff', 'Plasma cells'='#ff68a1')
  
b@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNF1, MMP2:CD4NADIR)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(TNF1=mean(TNF1), across(MMP2:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  pivot_longer(cols=MMP2:CD4NADIR, names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  group_by(Clusters_Final, Var)%>%
  group_modify(.f=~broom::tidy(cor.test(.x$TNF1, .x$Value, method = 'spearman', data = .x))) %>%
  ungroup()%>%
  mutate(padj=p.adjust(p.value, 'BH'))%>%
  arrange(padj)%>%mutate(text=ifelse(padj<0.25& p.value<0.05, Var, NA))%>%
  mutate(shape=ifelse(is.na(text), 'a', 'b'))%>%
  mutate(text=ifelse(str_detect(text, 'CRP|Intact|TNF'), text, NA))%>%
  ggplot(aes(x=estimate, y=-log10(p.value)))+
  geom_point(aes(color=Clusters_Final, shape=shape), position = position_jitter(.007, .007, seed = 111))+ 
  ggrepel::geom_text_repel(aes(label=text), size=4, min.segment.length = 0.05, position = position_jitter(.007, .007, seed = 111))+scale_shape_manual(values = c(46, 19))+
  theme_bw() +scale_color_manual(values=pal_b)+xlim(0.5, 1) +guides(shape='none')+xlab('rho')

#Save----
save.image('B.rdata')
SaveSeuratRds(b, 'Cell interaction/B.rds') 
