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
source('0. Spearman_boot.R')

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

load('DC.rdata')

#Recluster----
DimPlot(dc)

DefaultAssay(dc)<-'SCT' #Rerun clustering for dc 
dc<-RunPCA(dc)
ElbowPlot(dc)
dc<-RunUMAP(dc, dims = 1:15)

dc<-FindNeighbors(dc, dims = 1:15)
dc<-FindClusters(dc, algorithm = 2, random.seed = 111, resolution = 0.3)
DimPlot(dc, label = TRUE, label.size = 3)

DefaultAssay(dc)<-'SCT'
FeaturePlot(dc, features = c('CD1C', 'sct_CD14', 'IL3RA', 'THBD'))

pal_dc<-c('cDC0'='#66c2a5', 'cDC1'='#fc8d62', 'cDC2'='#a6d854', 'pDC0'='#e78ac3', 'pDC1'='#8da0cb')
 
DefaultAssay(dc)<-'SCT'
dc<-PrepSCTFindMarkers(dc)
dc.markers <- FindAllMarkers(dc, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 1)

dc.markers%>%
  filter(p_val_adj<0.05)%>%write.xlsx('DC.xlsx')

VlnPlot(dc, stack = T, features = c('CD1C', 'CLEC9A', 'FCER1A','CADM1', 'VCAN', 'S100A9', 'THBD', 'sct_CD14', 'IL3RA', 'CD209'), flip = T)+
  DimPlot(dc, label = T)

#Rename cluster----
metadata_dc<-dc@meta.data

Clusters_Final<-metadata_dc%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters =='0' ~ 'cDC0',
                                         seurat_clusters == '1' ~ 'cDC1',
                                         seurat_clusters == '2' ~ 'pDC0',
                                         seurat_clusters == '3' ~ 'cDC0', 
                                         seurat_clusters == '4' ~ 'cDC2', 
                                  seurat_clusters == '5' ~ 'pDC1'
                                  ))%>%
  select(Barcode, Clusters_Final)%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('cDC0', 'cDC1', 'cDC2', 'pDC0', 'pDC1')))%>%
  column_to_rownames('Barcode')

dc<-AddMetaData(dc, Clusters_Final)
DimPlot(dc, group.by = 'Clusters_Final', label=T)

Idents(dc)<-'Clusters_Final'

scplotter::CellDimPlot(dc, label = T, theme = "theme_blank", palcol=pal_dc, dot_size=2, 
                       label_insitu = T)

VlnPlot(dc, stack = T, features = c('CD1C', 'CLEC9A', 'FCER1A','CADM1', 'VCAN', 'S100A9', 'THBD', 'IL3RA', 'GZMB', 'AXL'), flip = T)+
  theme(axis.text.x = element_text(angle = 0))& xlab('Cluster') & guides(fill='none')

VlnPlot(dc, stack = T, features = c('CD1C', 'FCER1A','VCAN', 'S100A9', 'CADM1', 'CLEC9A', 'THBD', 'IL3RA'), 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) 

dc.markers_final <- FindAllMarkers(dc, only.pos = TRUE, min.pct = 0.1)

#Compare component of cluster and disease----
metadata_dc<-dc@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(dc, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar", colors_use = pal_dc) 

scCustomize::Proportion_Plot(seurat_object = subset(dc, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar", colors_use = pal_dc)

cellpercent<-metadata_dc%>% #summarize percentage of each cluster_final between control and disease group
  group_by(Group, Sample_Name, Clusters_Final, .drop = FALSE)%>%
  summarise(Count=n())%>%
  mutate(Percentage=Count/sum(Count)*100)%>%
  ungroup()

cellpercent%>%
  ggplot(aes(x=Group, y=Percentage, color=Group))+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  geom_jitter(aes(color=Group), width=0.2)+
  scale_color_manual(values = c("#BC3C29CC", "#0072B5CC"))+
  facet_wrap(~Clusters_Final) +theme_bw()#+stat_pwc()+scale_y_continuous(limits = c(0, 100))

#Save original data
cellpercent%>%
  pivot_wider(id_cols = Sample_Name, names_from = c(Clusters_Final, Group), values_from = Percentage)%>%
  write.xlsx('tmp.xlsx')

cvd_metadata<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Data_CVD.xlsx')%>%select(MACSID, Age, US_Risk, IMT_residual, USLESION:CD4NADIR), by=c('MACSID', 'US_Risk'))

#Add CVD data to metadata----
metadata_cvd_add<-metadata_dc%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

#Calculate cvd correlation----
cvd<-cvd_metadata%>%
  right_join(cellpercent, by='Sample_Name')%>%
  select(Clusters_Final, Sample_Name, Group, Age, Percentage, USLESION:CD4NADIR, IMT_residual)

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

#Plotting----
cvd%>%select(Clusters_Final, Percentage, Def5, Total, MMP9, MMP2, SERPINE1, ICAM1, FIBR, TNFARII)%>%
  filter(Clusters_Final=='pDC0')%>%
  ggplot(aes(x=log10(Percentage+1), y=log10(MMP9)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F, color='grey50', linetype='dashed')+
  theme_bw()+scale_color_manual(values=pal_dc)+guides(color='none')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.75, label.y.npc = 0.2)+
  xlab('log10(pDC0%+1)') |
  
  cvd%>%select(Clusters_Final, Percentage, Def5, Total, MMP9, MMP2, SERPINE1, ICAM1, FIBR, TNFARII)%>%
  filter(Clusters_Final=='pDC0')%>%
  ggplot(aes(x=log10(Percentage+1), y=log10(MMP2)))+
  geom_point()+
  geom_smooth(method = 'lm', se = F, color='grey50', linetype='dashed')+
  theme_bw()+scale_color_manual(values=pal_dc)+guides(color='none')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.75, label.y.npc = 0.2)+
  xlab('log10(pDC0%+1)')  
  

#Bootstrap these two
set.seed(111)
cor_boot <- cvd%>%select(Clusters_Final, Percentage, Def5, Total, MMP9, MMP2, SERPINE1, ICAM1, FIBR, TNFARII)%>%
  filter(Clusters_Final=='pDC0') %>%
  select(Clusters_Final, Percentage, MMP9, MMP2)%>%
  pivot_longer(cols= MMP9:MMP2, names_to = 'Item', values_to ='Value')%>%
  filter(is.na(Value)==FALSE)%>%
  group_by(Item) %>%
  group_modify(~ spearman_boot_ci(data = .x, x = 'Percentage', y = 'Value', R = 2000)) %>%
  ungroup()

cor_boot

#DEG----
##Classic DC ----
pseudo_cdc <- AggregateExpression(subset(dc, subset= Clusters_Final%in% c('cDC0', 'cDC1', 'cDC2')), assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-ncol(subset(dc, subset= Clusters_Final%in% c('cDC0', 'cDC1', 'cDC2')))

pseudo_cdc<-pseudo_cdc$RNA

pseudo_cdc<-pseudo_cdc%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_cdc)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) 

rownames(samples)<-colnames(pseudo_cdc) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_cdc,
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

VlnPlot(dc, features = c('HLA-DQB2', 'PLCL1'), group.by = 'Group', pt.size = 0, stack = T, flip=T)&
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20") &
  stat_summary(fun.y=median, geom='point', size=5, shape=20) & stat_compare_means()

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("GMT/h.all.v2025.1.Hs.symbols.gmt")
pathway.reactome<- gmtPathways("GMT/c2.cp.reactome.v2025.1.Hs.symbols.gmt")

fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))+scale_size_continuous(limits=c(50, 200))

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']], stats = ranks)+
  ggtitle('HALLMARK TNFA Signaling via NFKB')
 
#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) |

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))

write.xlsx(list(res_tidy, fgseaRes, fgseaRes_reactome), 'tmp.xlsx')

##Non-Classic dc ----
pseudo_pdc <- AggregateExpression(subset(dc, subset= Clusters_Final%in% c('pDC0', 'pDC1')), assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-ncol(subset(dc, subset= Clusters_Final%in% c('pDC0', 'pDC1')))

pseudo_pdc<-pseudo_pdc$RNA

pseudo_pdc<-pseudo_pdc%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_pdc)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) 

rownames(samples)<-colnames(pseudo_pdc) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_pdc,
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

VlnPlot(dc, features = c('HLA-DQB2', 'PTGDS'), group.by = 'Group', pt.size = 0, stack = T, flip=T)&
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20") &
  stat_summary(fun.y=median, geom='point', size=5, shape=20) & stat_compare_means()

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("GMT/h.all.v2025.1.Hs.symbols.gmt")
pathway.reactome<- gmtPathways("GMT/c2.cp.reactome.v2025.1.Hs.symbols.gmt")

set.seed(111)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))+scale_size_continuous(limits=c(50, 200))

fgsea::plotEnrichment(pathway = pathways.hallmark[['HALLMARK_TNFA_SIGNALING_VIA_NFKB']], stats = ranks)+
  ggtitle('HALLMARK TNFA Signaling via NFKB')

#Reactome
set.seed(132)
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0)) |
  
  fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))

VlnPlot(dc, features = c('PTGDS', 'HLA-DQA2', 'HLA-DQB2'), split.by = 'Group', pt.size = 0, stack = T, flip=T)+
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#Show aggregate comparison
dc<- AddModuleScore(dc, features = c('PTGDS'), name = 'PTGDS1', nbin=20)
dc<- AddModuleScore(dc, features = c('HLA-DQA2'), name = 'HLADQA2_', nbin=20)

dc@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, PTGDS11, HLADQA2_1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(PTGDS11=mean(PTGDS11), HLADQA2_1=mean(HLADQA2_1))%>%
  pivot_longer(cols = c(PTGDS11, HLADQA2_1), names_to = 'Gene', values_to = 'Score')%>%
  ggplot(aes(x=Clusters_Final, y=Score, color=Group))+
  geom_boxplot(outlier.shape = NA, fill=NA, width=0.5)+geom_jitter(height = 0, width = .1)+
  facet_wrap(~Gene, scales = 'free_y')+theme_bw()+
  scale_color_manual(values = c( "#BC3C29CC", "#0072B5CC"))+stat_pwc()+
  ylab('Module Score')+xlab('')

AverageExpression(dc, features = c('PTGDS', 'HLA-DQA2'), assays = 'SCT', slot = 'counts', return.seurat = F, 
                  group.by = c('Clusters_Final', 'Sample_Name'))%>%
  as.data.frame()%>%t()%>%as.data.frame()%>%
  rownames_to_column('Sample_Name')%>%
  as_tibble()%>% 
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>% #split sample_name to Clusters_Final and Sample_Name
  separate(Sample_Name, into=c('Clusters_Final', 'Sample_Name'), sep='_', extra = "merge", remove = FALSE)%>%  
  left_join(cvd%>%select(Sample_Name, Group)%>%unique(), by='Sample_Name')%>%
  pivot_longer(cols = c(PTGDS, `HLA-DQA2`), names_to = 'Gene', values_to = 'Expression')%>%
  ggplot(aes(x=Clusters_Final, y=Expression, color=Group))+
  geom_boxplot(outlier.shape = NA, fill=NA, width=0.6, position = position_dodge(width = 0.8))+
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = .8))+
  facet_wrap(~Gene, scales = 'free_y')+theme_bw()+
  scale_color_manual(values = c( "#BC3C29CC", "#0072B5CC"))+stat_pwc(label = 'p.adj.format', p.adjust.by = 'panel', p.adjust.method = 'BH')  

FeaturePlot(dc, features = c('PTGDS', 'HLA-DQA2'), split.by = 'Group', order=T)

AverageExpression(subset(dc, subset= Clusters_Final=='pDC0'), features = 'PTGDS', assays = 'SCT', slot = 'counts', return.seurat = F, 
                  group.by = c('Sample_Name'))%>%
  as.data.frame()%>%t()%>%as.data.frame()%>%
  rownames_to_column('Sample_Name')%>%
  as_tibble()%>%dplyr::rename('Gene'='V1')%>%
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>%
  left_join(cvd%>%select(Sample_Name, Group, MMP9, Intact, Total, Def5, Def3, MMP9, IMT)%>%unique(), by='Sample_Name')%>%
  with(., cor.test(Gene, Age, method='spearman'))

AverageExpression(subset(dc, subset= Clusters_Final=='pDC0'), features = 'PTGDS', assays = 'SCT', slot = 'count', return.seurat = F, 
                  group.by = c('Sample_Name'))%>%
  as.data.frame()%>%t()%>%as.data.frame()%>%
  rownames_to_column('Sample_Name')%>%
  as_tibble()%>%dplyr::rename('Gene'='V1')%>%
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>%
  left_join(cvd%>%select(Sample_Name, Group, MMP9, Intact, Total, Def5, Def3, MMP9, IMT)%>%unique(), by='Sample_Name')%>%
  ggplot(aes(x=Group, y=Gene, color=Group))+
  geom_jitter(height = 0, width = .1)+
  geom_boxplot(outlier.shape = NA, fill=NA, width=0.5)+ylab('PTGDS Expression (SCT)')+
  theme_bw()+scale_color_manual(values = c( "#BC3C29CC", "#0072B5CC"))+stat_pwc() +

AverageExpression(subset(dc, subset= Clusters_Final=='pDC0'), features = 'PTGDS', assays = 'SCT', slot = 'counts', return.seurat = F, 
                  group.by = c('Sample_Name'))%>%
  as.data.frame()%>%t()%>%as.data.frame()%>%
  rownames_to_column('Sample_Name')%>%
  as_tibble()%>%dplyr::rename('Gene'='V1')%>%
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>%
  left_join(cvd%>%select(Sample_Name, Group, Age, MMP9, Intact, Total, Def5, Def3, MMP9, IMT_residual)%>%unique(), by='Sample_Name')%>%
  ggplot(aes(x=Gene, y=IMT_residual))+
  geom_point()+
  geom_smooth(method = 'lm', se = F, color='grey50', linetype='dashed')+
  theme_bw()+scale_color_manual(values=pal_dc)+guides(color='none')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.75, label.y.npc = 0.2)+
  xlab('PTGDS Expression (SCT)')+ylab('IMT Residual')

set.seed(111)
cor_boot1 <- AverageExpression(subset(dc, subset= Clusters_Final=='pDC0'), features = 'PTGDS', assays = 'SCT', slot = 'counts', return.seurat = F, 
                               group.by = c('Sample_Name'))%>%
  as.data.frame()%>%t()%>%as.data.frame()%>%
  rownames_to_column('Sample_Name')%>%
  as_tibble()%>%dplyr::rename('Gene'='V1')%>%
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>%
  left_join(cvd%>%select(Sample_Name, Group, IMT_residual)%>%unique(), by='Sample_Name') %>%
  with(., spearman_boot_ci(data=., x = 'Gene', y = 'IMT_residual', R = 2000))  

cor_boot1
 
write.xlsx(list(res_tidy, fgseaRes, fgseaRes_reactome), 'tmp.xlsx')

FeaturePlot(dc, c('PTGDS', 'HLA-DQA2', 'HLA-DQB2', 'HLA-DRB6'), order=T, pt.size = 0.5, keep.scale = 'feature',
            split.by = 'Group', min.cutoff = 0.5, max.cutoff = 10)

#Save----
save.image('DC.rdata')
 