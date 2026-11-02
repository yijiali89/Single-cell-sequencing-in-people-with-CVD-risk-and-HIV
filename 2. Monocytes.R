library(tidyverse)
library(Seurat)
library(SeuratObject)
library(readxl)
library(openxlsx)
library(conflicted)
library(SingleR)
library(Seurat.utils)
library(DESeq2)
library(fgsea)
library(patchwork)
library(scplotter)
library(gee)
source('0. Spearman_boot.R')

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::count)

load('Monocyte_DC.rdata')

#Recluster----
DimPlot(mono)

DefaultAssay(mono)<-'SCT' #Rerun clustering for monocytes
mono<-RunPCA(mono)
ElbowPlot(mono)
mono<-RunUMAP(mono, dims = 1:20)

mono<-FindNeighbors(mono, dims = 1:20)
mono<-FindClusters(mono, algorithm = 2, random.seed = 121, resolution = 0.3)

DimPlot(mono, label = TRUE, label.size = 3, group.by = 'seurat_clusters') /

DotPlot(mono,  features = c('CD14', 'FCGR3A', 'CD1C', 'IL3RA', 'VCAN', 'S100A8', 
                            'THBS1', 'THBD', 'IL1B', 'IL1RN', 'APOBEC3A', 'EGR1', 'BST2',
                            'PDE4B', 'CXCL2', 'CCR2'), 
        cluster.idents = T, scale=F)

DefaultAssay(mono)<-'SCT'
FeaturePlot(mono, features = c('CD14', 'FCGR3A', 'LYZ', 'CD1C', 'IL3RA'))

DefaultAssay(mono)<-'SCT'
Idents(mono)<-'seurat_clusters'
mono<-PrepSCTFindMarkers(mono)
monocyte.markers <- FindAllMarkers(mono, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.5)

monocyte.markers%>%
  filter(pct.1>0.5, pct.2<0.5, avg_log2FC>.5, p_val_adj<0.01)%>%write.xlsx('Mono.xlsx')

mono_5vs7<-FindMarkers(mono, ident.1 = '5', ident.2 = '7') #APOBEC3A difference otherwise not much diff
mono_0vs8<-FindMarkers(mono, ident.1 = '0', ident.2 = '8') #Not much difference
mono_0vs1<-FindMarkers(mono, ident.1 = '0', ident.2 = '1') #Not much difference
mono_0vs3<-FindMarkers(mono, ident.1 = '0', ident.2 = '3') #Not much difference

mono_4vs9<-FindMarkers(mono, ident.1 = '4', ident.2 = '9') #Not much difference

DimPlot(subset(mono, seurat_clusters=='1'))+DimPlot(subset(mono, seurat_clusters=='3'))

FeaturePlot(mono, c('AREG'), order=T)

pal_mono<-c('CM0'='#A6CEE3', 'CM1'='#FB9A99', 'CM2'= '#FDBF6F', 
            'NCM0'="#D43F3ACC" , 'NCM1'= "#46B8DACC", 
            'cDC'="#357EBDCC", 'pDC'="#9632B8CC")

VlnPlot(mono, stack = T, features = c('CD14', 'THBS1', 'FCGR3A', 'IL1B', 'IL1RN', 'VEGFA'), flip = T)
 
#Rename cluster----
metadata_mono<-mono@meta.data

Clusters_Final<-metadata_mono%>%rownames_to_column('Barcode')%>%
  mutate(Clusters_Final=case_when(seurat_clusters %in% c('0', '1', '3', '8') ~ 'CM0', #Not much difference except THBS1 but Cluster 3 and Cluster0 basically superimpose on UMAP
                                         seurat_clusters %in% c('2') ~ 'CM1',
                                         seurat_clusters %in% c('5', '7') ~ 'CM2', 
                                         seurat_clusters == '4' ~ 'NCM0',
                                         seurat_clusters == '9' ~ 'NCM0', 
                                         seurat_clusters == '11' ~ 'NCM1',
                                      seurat_clusters == '6' ~ 'cDC', 
                                      seurat_clusters == '10' ~ 'pDC'))%>%
  select(Barcode, Clusters_Final)%>%
  mutate(Clusters_Final=factor(Clusters_Final, levels=c('CM0', 'CM1', 'CM2', 'CM3', 'CM4','NCM0',  'NCM1', 
                                                        'NCM2', 'cDC', 'pDC')))%>%
  column_to_rownames('Barcode')

mono<-AddMetaData(mono, Clusters_Final)

DimPlot(mono, group.by = 'Clusters_Final', label=T)+ DimPlot(mono, label=T)

Idents(mono)<-'Clusters_Final'
mono<-PrepSCTFindMarkers(mono)
monocyte.markers_final <- FindAllMarkers(mono, only.pos = T, min.pct = 0.1, logfc.threshold = 1)

write.xlsx(monocyte.markers_final, 'Mono.xlsx')

scplotter::CellDimPlot(mono, label = T, theme = "theme_blank", palcolor = pal_mono, dot_size=5, 
                       label_insitu = T)

VlnPlot(mono, stack = T, features = c('CD14', 'FCGR3A', 'CD1C', 'IL3RA', 'IL1RN', 'THBS1', 
                                      'SERPINB2','IL1B', 'CXCL2', 'C1QA'), 
        flip = T) & stat_summary(fun = median, geom = "point", color = 'black', 
                             position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) & theme(axis.text.x = element_text(angle = 0)) & xlab('')
 
DefaultAssay(mono)<-'SCT'
((FeaturePlot(subset(mono, subset= Group=='High'), features = c('sct_CD14', 'FCGR3A', 'CD1C', 'IL3RA', 'THBS1', 'IL1B', 'SERPINB2', 'C1QA'), 
            keep.scale = 'all', min.cutoff = 0, max.cutoff = 4.5, cols =c('grey90', "#BC3C29CC"), raster=T, pt.size = 3, raster.dpi = c(600, 600),
                reduction = "umap", ncol = 4)+plot_annotation('US high-risk')) /

(FeaturePlot(subset(mono, subset= Group=='Low'), features = c('sct_CD14', 'FCGR3A', 'CD1C', 'IL3RA', 'THBS1', 'IL1B','SERPINB2',  'C1QA'), 
            keep.scale = 'all',  min.cutoff = 0, max.cutoff = 4.5, cols =c('grey90',"#0072B5CC"), raster=T, pt.size = 3, raster.dpi = c(600, 600),
            reduction = "umap", ncol = 4) +plot_annotation('US low-risk'))) & NoAxes()


ms1sig<-read_xlsx('MS1.xlsx')%>%select(Gene)%>%unlist()

mono <- AddModuleScore(mono, features = list(ms1sig[1:100]), name = 'MS', assay = 'SCT')

VlnPlot(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), 
        cols = c('High'= "#BC3C29CC", 'Low'="#0072B5CC"),  
        features=c('MS1', 'IL1R2', 'RETN', 'ALOX5AP', 'SERPINB1'), stack = T, flip = T, split.by = 'Group')& 
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9))  &
  stat_compare_means(label = 'p.format', method = 'wilcox') &
  scale_y_continuous(expand = expansion(c(.1, .1)))

#Generate individual level MS1 module score comparison
MS1_CM0<-mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, MS1)%>%
  filter(Clusters_Final %in% c('CM0'))%>%
  arrange(Sample_Name)%>%
  mutate(Group=factor(Group))%>%
  mutate(ID=as.numeric(as.factor(Sample_Name)))

mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, MS1)%>%
  filter(Clusters_Final %in% c('CM0', 'CM1', 'CM2'))%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(MS1=mean(MS1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=MS1, fill=Group))+
  geom_violin(position = position_dodge(width = 0.9), alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) +
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) +
  stat_compare_means(label = 'p.format', method = 'wilcox')


gee_MS1_CM0<-gee(MS1 ~ Group, id = ID, data = MS1_CM0, family = gaussian, corstr = "exchangeable", silent=FALSE)

VlnPlot(mono, features = c('THBS1', 'IL1B', 'SERPINB1', 'AREG', 'CXCL8', 'TNF'), stack=T, flip=T, split.by = 'Group', pt.size = 0, cols = c("#0072B5CC", "#BC3C29CC")) &
  theme(axis.text.x = element_text(angle = 0))& xlab('') & 
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) & stat_compare_means(label = 'p.format')

VlnPlot(mono, features = c('MS1'), split.by = 'Group', pt.size = 0, cols = c("#0072B5CC", "#BC3C29CC")) &
  theme(axis.text.x = element_text(angle = 0))& xlab('') & 
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) & stat_compare_means(label = 'p.format')
 
#Compare component of cluster and disease----
metadata_mono<-mono@meta.data

scCustomize::Proportion_Plot(seurat_object = subset(mono, subset= Group=='High'), 
                             split.by = 'Sample_Name', plot_type = "bar", colors_use = pal_mono) 

scCustomize::Proportion_Plot(seurat_object = subset(mono, subset= Group=='Low'), 
                             split.by = 'Sample_Name', plot_type = "bar", colors_use = pal_mono)

cellpercent<-metadata_mono%>% #summarize percentage of each cluster_final between control and disease group
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(Count=n())%>%
  filter(Clusters_Final %in% c('pDC', 'cDC', 'NCM0', 'NCM1')==FALSE)%>%
  mutate(Percentage=Count/sum(Count)*100)%>%
  ungroup()

cellpercent%>%
  ggplot(aes(x=Group, y=Percentage, color=Group))+
  geom_boxplot(outlier.shape = NA, alpha=0)+
  geom_jitter(aes(color=Group), width=0.2)+
  scale_color_manual(values = c("#BC3C29CC", "#0072B5CC"))+
  facet_wrap(~Clusters_Final, scales = 'free') +theme_bw()+stat_pwc()

cvd_metadata<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Data_CVD.xlsx')%>%select(MACSID, VL, US_Risk, USLESION:CD4NADIR), by=c('MACSID', 'US_Risk'))%>%
  mutate(VL_Blip=ifelse(VL>50, 'Yes', 'No'))

#Senescence score----
pathway.sen<-gmtPathways('GMT/SAUL_SEN_MAYO.v2025.1.Hs.gmt')

mono<-AddModuleScore(mono, features = list(pathway.sen$SAUL_SEN_MAYO), name = 'SASP', assay = 'SCT')

FeaturePlot(mono, c('SASP1', 'MS1'), split.by = 'Group', #cols = c('grey90', "#0072B5CC", "#BC3C29CC"), keep.scale = 'all',
            raster=T, pt.size = 3, raster.dpi = c(600, 600), min.cutoff = 'q5', max.cutoff = 'q95',
            reduction = "umap", ncol = 1 ) +plot_annotation('SASP score')
mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, MS1)%>%
  filter(Clusters_Final %in% c('CM0', 'CM1', 'CM2'))%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(MS1=mean(MS1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=MS1, fill=Group))+
  geom_violin(position = position_dodge(width = 0.9), alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) +
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) +
  stat_compare_means(label = 'p.format', method = 'wilcox') +

mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, SASP1)%>%
  filter(Clusters_Final %in% c('CM0', 'CM1', 'CM2'))%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(SASP1=mean(SASP1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=SASP1, fill=Group))+
  geom_violin(position = position_dodge(width = 0.9), alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = .9), alpha=0.5)+
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) +
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9)) +
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  plot_layout(ncol = 2, guides = 'collect') & theme(legend.position = 'bottom')

VlnPlot(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), 
        cols = c('High'= "#BC3C29CC", 'Low'="#0072B5CC"),  
        features=c('MS1', 'SASP1'), stack = T, flip = T, split.by = 'Group')& 
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) &
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9))  &
  stat_compare_means(label = 'p.format', method = 'wilcox') &
  scale_y_continuous(expand = expansion(c(.1, .1)))

#Add CVD data to metadata----
metadata_cvd_add<-metadata_mono%>%
  rownames_to_column('Barcode')%>%
  select(Barcode, Sample_Name)%>%
  left_join(cvd_metadata%>%select(Sample_Name, MACSID, US_Risk, USLESION:CD4NADIR), by=c('Sample_Name'))%>%
  column_to_rownames('Barcode')

mono<-AddMetaData(mono, metadata_cvd_add)

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
  pivot_longer(cols = c(Group, contains('DI'), STEN4, STEN34), 
               names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  group_by(Clusters_Final, Var)%>%
  group_modify(.f=~broom::tidy(wilcox.test(.x$Percentage~.x$Value, data = .x))) %>%
  ungroup()%>%
  mutate(Padj=p.adjust(p.value, 'BH'))

#Plotting----
cvd%>%select(Clusters_Final, Percentage, MMP9, MMP2, CRP, TNFARI)%>%
  filter(Clusters_Final %in% c('CM1', 'CM2'))%>%
  pivot_longer(cols = c(MMP9, MMP2, CRP, TNFARI), 
               names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  filter((Clusters_Final=='CM1') | (Clusters_Final=='CM2' & Var %in% c('MMP2', 'MMP9')))%>%
  mutate(Item=paste0(Clusters_Final, ' percent vs. ', Var))%>%
  #mutate(Item=factor(Item, levels=c(paste0('M1 percent vs. ', c('Total', 'Def5', 'MMP9', 'MMP2', 'SERPINE1', 'ICAM1', 'FIBR')), 
  #                                  paste0('M3 percent vs. ', c('ICAM1', 'TNFARII')))))%>%
  ggplot(aes(x=log10(Percentage+1), y=log10(Value)))+
  geom_point(aes(color=Clusters_Final), shape=19)+
  geom_smooth(method = 'lm', se = F, color='grey50', linetype='dashed')+
  facet_wrap(~Item, scales = 'free')+
  theme_bw()+scale_color_manual(values=pal_mono)+guides(color='none')+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.5, label.y.npc = 0.5)

#Final results CI bootstrap
set.seed(121)
cor_boot <- cvd %>%
  select(Clusters_Final, Percentage, MMP9, MMP2, CRP, TNFARI) %>%
  filter(Clusters_Final %in% c('CM1', 'CM2')) %>%
  pivot_longer(cols = c(MMP9, MMP2, CRP, TNFARI),
               names_to = 'Var', values_to = 'Value') %>%
  filter(!is.na(Value), !is.na(Percentage)) %>%
  filter((Clusters_Final == 'CM1') | (Clusters_Final == 'CM2' & Var %in% c('MMP2', 'MMP9'))) %>%
  mutate(Item = paste0(Clusters_Final, ' percent vs. ', Var)) %>%
  group_by(Item) %>%
  group_modify(~ spearman_boot_ci(data = .x, x = 'Percentage', y = 'Value', R = 2000)) %>%
  ungroup()

cor_boot
   
#DEG----
##Classic monocytes----
DimPlot(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')))
pseudo_cm <- AggregateExpression(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), assays = "RNA", 
                                       slot = 'counts', return.seurat = F, 
                                       group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<- dim(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')))[2]

pseudo_cm<-pseudo_cm$RNA

pseudo_cm<-pseudo_cm%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_cm)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL_Blip)%>%
              mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')
  
rownames(samples)<-colnames(pseudo_cm) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_cm,
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
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC','2. Padj 0.05-0.1'= '#46B8DACC', 'NS'= 'grey50'))

hla_df<-AggregateExpression(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), assays = "SCT", 
                    slot = 'counts', return.seurat = F, 
                    group.by = c('Sample_Name')) %>% 
  as.data.frame()%>%t()%>%as.data.frame()%>% 
  rownames_to_column('Sample_Name')%>%
  as_tibble() %>%
  select(Sample_Name, 'HLA-DQA2', 'HLA-DQB2')%>%
  mutate(Sample_Name=str_remove(Sample_Name, 'SCT.')%>%str_replace(pattern='.hs', replacement = '_hs'))%>%
  left_join(cvd_metadata%>%select(Sample_Name, Group=US_Risk, Batch)) 
  
hla_df%>%
  mutate(Group=case_when(
    Group==1 ~ 'High-risk',
    Group==0 ~ 'Low-risk'
  ))%>%
  mutate(Group=factor(Group, levels=c('Low-risk', 'High-risk')), Batch=factor(Batch))%>%
  pivot_longer(cols = c('HLA-DQA2', 'HLA-DQB2'), names_to = 'Gene', values_to = 'Expression')%>%
  ggplot(aes(x=Group, y=log10(Expression+1), fill=Group))+
  geom_boxplot(outlier.shape = NA, alpha=0.1)+geom_jitter(aes(color=Batch), width = 0.2, alpha = 0.5)+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+
  theme(axis.text.x = element_text(angle = 0))+
  stat_compare_means(label = 'p.format', method = 'wilcox')+
  facet_grid(~Gene, scales = 'free')+
  scale_y_continuous(expand = expansion(c(.1, .1)))+
  theme_bw()

VlnPlot(subset(mono, subset= Clusters_Final %in% c( 'CM0' )),
        features = c('ZFP36', 'NFKBIA'), group.by = 'Group', pt.size = 0, stack=T, flip=T)+
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20") +
  stat_summary(fun.y=median, geom='point', size=5, shape=20) +stat_compare_means()

res2 <- res_tidy %>%
  dplyr::select(row, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(row) %>%
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("GMT/h.all.v2025.1.Hs.symbols.gmt")
pathway.reactome<- gmtPathways("GMT/c2.cp.reactome.v2025.1.Hs.symbols.gmt")

set.seed(12)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))

gene_mtor<-fgseaRes$leadingEdge[31]%>%unlist()

mono<-AddModuleScore(mono, features = gene_mtor, name = 'HM_MTORC')

VlnPlot(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), cols = c("#BC3C29CC", "#0072B5CC"),
        c('HM_MTORC1', 'TNF'), 
        split.by  = 'Group', stack=T, flip=T)+
    theme_minimal()+
  stat_compare_means()+
  scale_y_continuous(expand = expansion(c(.1, .1)))+
  stat_summary(fun = median, geom = "point", color = 'black', 
               position = position_dodge(width = 0.9), size = 1.5) +
  stat_summary(fun= function(x) quantile(x, c(.25, .75)), geom='line',
               position = position_dodge(width = 0.9))
  

VlnPlot(subset(mono, subset= Group=='Low'), c('HM_MTORC1', 'TNF'), group.by  = 'orig.ident', stack=T, flip=T)+
  VlnPlot(subset(mono, subset= Group=='High'), c('HM_MTORC1', 'TNF'), group.by  = 'orig.ident', stack=T, flip=T)

AggregateExpression(subset(mono, subset= Clusters_Final!='M3'), assays = 'SCT', features = 'ZFP36', 
                    slot = 'counts', return.seurat = F, group.by = c('Group', 'orig.ident', 'Sample_Name'))%>%
  as.data.frame()%>%rownames_to_column('Gene')%>%
  pivot_longer(cols = -Gene, names_to = 'Sample', values_to = 'Expression')%>%
  mutate(Group=str_split(Sample, '_', n=3, simplify = T)[,1],
         Sample_Name=str_split(Sample, '_', n=3, simplify = T)[,3])%>%
  ggplot(aes(x=Group, y=log10(Expression+1), fill=Group))+
  geom_boxplot(outlier.shape = NA)+geom_jitter(width = 0.2, alpha = 0.5)+
  scale_fill_manual(values=c("#0072B5CC", "#BC3C29CC"))+
  theme(axis.text.x = element_text(angle = 90))+stat_compare_means()

#Reactome
set.seed(111)
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

((fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(50, 100, 150, 200), limits = c(10, 200))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(-2.5, 2.5)+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))) / plot_spacer()) |

(fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))+
   scale_size_continuous(breaks = c(50, 100, 200, 400, 800), limits = c(10, 800))+
   scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 30)))

#For CROI only
fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05, NES>1.5)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL'))%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))+
  scale_size_continuous(transform = 'sqrt', breaks = c(50, 100, 200, 400, 800), limits = c(20, 1000))+
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 30))+
  xlim(1.4, 2.25)

write.xlsx(list(res_tidy, fgseaRes, fgseaRes_reactome), 'tmp.xlsx')

Oxphos_leadingedge<-fgseaRes$leadingEdge[36]

###MS1----
library(org.Hs.eg.db)
ms1list<-list(ms1sig)
names(ms1list)<-'MS1'
fgseaRes_ms1 <- fgseaMultilevel(pathways=ms1list, stats=ranks, nPermSimple = 10000)%>%
  mutate(leadingEdge=paste0(unlist(leadingEdge),collapse = "/"))

fgsea::plotEnrichment(ms1list$MS1, ranks)

enrich_df <- data.frame(ID = 'GO:0000001',
                        Description = 'MS1',
                        setSize = fgseaRes_ms1$size,
                        enrichmentScore = fgseaRes_ms1$ES,
                        NES = fgseaRes_ms1$NES,
                        pvalue = fgseaRes_ms1$pval,
                        p.adjust = fgseaRes_ms1$padj,
                        core_enrichment = fgseaRes_ms1$leadingEdge)

rownames(enrich_df) <- 'GO:0000001'

en <- GseaVis::dfGO2gseaResult(enrich.df = enrich_df,keytype = 'SYMBOL',
                               geneList = sort(ranks, decreasing = T),
                               OrgDb = 'org.Hs.eg.db')


en <- DOSE::setReadable(x = en, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")

en@geneSets$`GO:0000001`<- c(ms1sig)

enrichplot::gseaplot2(x = en,geneSetID = enrich_df$ID[1], pvalue_table = TRUE,
                      title = enrich_df$Description[1]) 

GseaVis::gseaNb(object = en, geneSetID = 'GO:0000001', subPlot = 2, addPval = T,
                newGsea = F, addGene = c('CCL7', 'IL1RN', 'TNFSF14', 'SERPINB2'),
                addPoint = F)

###SASP----
fgseaRes_sasp <- fgseaMultilevel(pathways=pathway.sen, stats=ranks, nPermSimple = 10000)%>%
  mutate(leadingEdge=paste0(unlist(leadingEdge),collapse = "/"))

fgsea::plotEnrichment(pathway.sen$SAUL_SEN_MAYO, ranks)

enrich_df <- data.frame(ID = 'GO:0000001',
                        Description = 'MS1',
                        setSize = fgseaRes_ms1$size,
                        enrichmentScore = fgseaRes_ms1$ES,
                        NES = fgseaRes_ms1$NES,
                        pvalue = fgseaRes_ms1$pval,
                        p.adjust = fgseaRes_ms1$padj,
                        core_enrichment = fgseaRes_ms1$leadingEdge)

rownames(enrich_df) <- 'GO:0000001'

en <- GseaVis::dfGO2gseaResult(enrich.df = enrich_df,keytype = 'SYMBOL',
                               geneList = sort(ranks, decreasing = T),
                               OrgDb = 'org.Hs.eg.db')


en <- DOSE::setReadable(x = en, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")

en@geneSets$`GO:0000001`<- c(ms1sig)

enrichplot::gseaplot2(x = en,geneSetID = enrich_df$ID[1], pvalue_table = TRUE,
                      title = enrich_df$Description[1]) 

GseaVis::gseaNb(object = en, geneSetID = 'GO:0000001', subPlot = 2, addPval = T,
                newGsea = F, addGene = c('CCL7', 'IL1RN', 'TNFSF14', 'SERPINB2'),
                addPoint = F)

##Classic monocytes including VLBlip----
DimPlot(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')))
pseudo_cm <- AggregateExpression(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')), assays = "RNA", 
                                 slot = 'counts', return.seurat = F, 
                                 group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<- dim(subset(mono, subset= Clusters_Final %in% c('CM0', 'CM1', 'CM2')))[2]

pseudo_cm<-pseudo_cm$RNA

pseudo_cm<-pseudo_cm%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_cm)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL_Blip)%>%
              mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')

rownames(samples)<-colnames(pseudo_cm) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_cm,
                            colData = samples,
                            design = ~ Group+Batch+VL_Blip) #Control batch

dds <- DESeq(dds)
res<-results(dds, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy <- results(dds, tidy = T, contrast = c('Group', 'High', 'Low'), alpha = 0.05, independentFiltering = FALSE)
res_tidy_vlblip <- results(dds, tidy = T, contrast = c('VL_Blip', 'Yes', 'No'), alpha = 0.05, independentFiltering = FALSE)

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
  theme_minimal()+scale_color_manual(values=c('1. Padj<0.05'='#D43F3ACC', '#46B8DACC', 'NS'= 'grey50'))

res_tidy_vlblip%>%
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
  theme_minimal()+scale_color_manual(values=c('#D43F3ACC', '#46B8DACC', 'grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

pathways.hallmark <- gmtPathways("GMT/h.all.v2025.1.Hs.symbols.gmt")
pathway.reactome<- gmtPathways("GMT/c2.cp.reactome.v2025.1.Hs.symbols.gmt")

set.seed(12)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))

#Reactome
set.seed(111)
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

(fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
    mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace('_', ' ')%>%str_to_title())%>%
    ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(50, 100, 150, 200), limits = c(10, 200))+
    geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(-2.5, 2.5)+
    scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))) /
  
  (fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
     filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL|NFKB'))%>%
     mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
     ggplot(aes(x=NES, y=reorder(pathway, NES)))+
     geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
     scale_color_viridis_c(begin = 0.3, end=0.9)+geom_vline(aes(xintercept = 0))+
     scale_size_continuous(breaks = c(50, 100, 400), limits = c(10, 800))+
     scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 30)))+plot_layout(ncol=1, heights = c(2,1))

##Non-Classic monocytes----
DimPlot(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')))
pseudo_ncm <- AggregateExpression(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')) , assays = "RNA", 
                                 slot = 'counts', return.seurat = F, 
                                 group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')))[2]

pseudo_ncm<-pseudo_ncm$RNA

pseudo_ncm<-pseudo_ncm%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_ncm)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) 

rownames(samples)<-colnames(pseudo_ncm) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_ncm,
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
  theme_minimal()+scale_color_manual(values=c('#D43F3ACC', '#46B8DACC', 'grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)
 
set.seed(121)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, option = 'A')+geom_vline(aes(xintercept = 0))
 
#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
  filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|VIR|INTERLEUKIN|INFECT|INFLAMM|RESPIRAT|TNF|NF_KB|VIRAL|OXIDATI|MHC|MACROPHAGE|MONOCYTE'))%>%
  mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_size_continuous(breaks=c(1, 5, 10, 50, 100, 200), limits = c(5, 300))+
  scale_color_viridis_c(begin = 0.3, end=0.9, option='A')+geom_vline(aes(xintercept = 0))+
  scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "),
                                               width = 40))

##Non-Classic monocytes with VLBlip adjustment----
DimPlot(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')))
pseudo_ncm <- AggregateExpression(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')) , assays = "RNA", 
                                  slot = 'counts', return.seurat = F, 
                                  group.by = c('Group', 'orig.ident', 'Sample_Name')) #Aggregate by group and sample
nsample<-dim(subset(mono, subset = Clusters_Final %in% c('NCM0', 'NCM1')))[2]

pseudo_ncm<-pseudo_ncm$RNA

pseudo_ncm<-pseudo_ncm%>%as.data.frame()%>%
  rownames_to_column('Gene')%>%
  rowwise()%>%
  mutate(Genesum=sum(c_across(`High_CVD123-1_SampleTag09-hs`:`Low_CVD4Aprip14_SampleTag10-hs`)))%>%
  ungroup()%>%
  filter(Genesum>nsample*.05)%>%select(-Genesum)%>% #Filter those genes that were expressed more than 5% in sample
  column_to_rownames('Gene')

samples<-colnames(pseudo_cm)%>%
  str_split('_', n = 3, simplify = T)%>%
  as.data.frame()%>%
  select(Group=V1, Sample_Name=V3, Batch=V2) %>%
  left_join(cvd_metadata%>%select(Sample_Name, VL_Blip)%>%
              mutate(Sample_Name=str_replace(Sample_Name, '_hs', '-hs')), by='Sample_Name')

rownames(samples)<-colnames(pseudo_ncm) #Prepare DESeq2 label

dds<-DESeqDataSetFromMatrix(countData = pseudo_ncm,
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
  theme_minimal()+scale_color_manual(values=c('#D43F3ACC', '#46B8DACC', 'grey50'))

res2 <- res_tidy %>% 
  dplyr::select(row, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(row) %>% 
  summarize(stat=mean(stat)) #Prepare for fgsea

ranks <- deframe(res2)

set.seed(121)
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks, nPermSimple =10000)

fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
  ggplot(aes(x=NES, y=reorder(pathway, NES)))+
  geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
  scale_color_viridis_c(begin = 0.3, end=0.9, option = 'A')+geom_vline(aes(xintercept = 0))

#Reactome
fgseaRes_reactome <- fgseaMultilevel(pathways=pathway.reactome, stats=ranks, nPermSimple = 10000)

(fgseaRes%>%as_tibble()%>%filter(padj<0.05)%>%
    mutate(pathway=str_remove(pathway, 'HALLMARK_')%>%str_replace('_', ' ')%>%str_to_title())%>%
    ggplot(aes(x=NES, y=reorder(pathway, NES)))+scale_size_continuous(breaks = c(50, 100, 150, 200), limits = c(10, 200))+
    geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+xlim(-2.5, 2.5)+
    scale_color_viridis_c(begin = 0.3, end=0.9, option = 'A')+geom_vline(aes(xintercept = 0))) /
  
  (fgseaRes_reactome%>%as_tibble()%>%filter(padj<0.05)%>%
     filter(str_detect(pathway, pattern = 'IMMUN|HIV|IFN|INTERLEUKIN|INFEC|RESPIRAT|TNF|NF_KB|VIRAL|NFKB'))%>%
     filter(pathway!='REACTOME_RESPIRATORY_ELECTRON_TRANSPORT'& pathway!='REACTOME_SARS_COV_INFECTIONS')%>%
     mutate(pathway=str_remove(pathway, 'REACTOME_')%>%str_replace_all('_', ' ')%>%str_to_title())%>%
     ggplot(aes(x=NES, y=reorder(pathway, NES)))+
     geom_point(aes(color=-log10(pval), size=size))+ylab(NULL)+theme_minimal()+
     scale_color_viridis_c(begin = 0.3, end=0.9, option = 'A')+geom_vline(aes(xintercept = 0))+
     scale_size_continuous(breaks = c(50, 100, 400), limits = c(10, 800))+
     scale_y_discrete(labels = function(x) str_wrap(str_replace_all(x, "foo" , " "), width = 30)))#+plot_layout(ncol=1, heights = c(2,1))

#MTOR signal----
enrich_df <- data.frame(ID = 'GO:0000001',
                        Description = 'MTORC1 Signaling',
                        setSize = fgseaRes$size[31],
                        enrichmentScore = fgseaRes$ES[31],
                        NES = fgseaRes$NES[31],
                        pvalue = fgseaRes$pval[31],
                        p.adjust = fgseaRes$padj[31],
                        core_enrichment = paste0(unlist(fgseaRes$leadingEdge[31]), collapse = '/'))

rownames(enrich_df) <- 'GO:0000001'

library(org.Hs.eg.db)
en <- GseaVis::dfGO2gseaResult(enrich.df = enrich_df,keytype = 'SYMBOL',
                               geneList = sort(ranks, decreasing = T),
                               OrgDb = 'org.Hs.eg.db')


en <- DOSE::setReadable(x = en, OrgDb = org.Hs.eg.db, keyType = "SYMBOL")

en@geneSets$`GO:0000001`<- pathways.hallmark$HALLMARK_MTORC1_SIGNALING

enrichplot::gseaplot2(x = en,geneSetID = enrich_df$ID[1], pvalue_table = TRUE,
                      title = enrich_df$Description[1]) 

GseaVis::gseaNb(object = en, geneSetID = 'GO:0000001', subPlot = 2, addPval = T,
                newGsea = F,  addGene = 'NAMPT',
                addPoint = F)


gene_mtor<-pathways.hallmark[str_detect(names(pathways.hallmark), 'MTORC')]%>%unlist()%>%as.vector()
names(gene_mtor)<-NULL

gene_mtor<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_mtor)

mono@meta.data <- mono@meta.data[, !grepl("MTOR", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_mtor)), assay = 'SCT', name = 'MTORC')

FeaturePlot(mono, features = 'MTORC1', split.by = 'Group')

#Aggregate MTOR score by participant
mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, MTORC1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(MTORC1=mean(MTORC1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=MTORC1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox') +
  
  mono@meta.data%>%
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
  stat_compare_means(label = 'p.format', method = 'wilcox')+
  
  mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, Complement1)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  summarise(Complement1=mean(Complement1))%>%
  ungroup()%>%
  ggplot(aes(x=Clusters_Final, y=Complement1, fill=Group))+
  geom_boxplot(position = position_dodge(width = 0.9), outlier.shape = NA, alpha=0.5)+
  geom_jitter(aes(color=Group), position = position_jitterdodge(jitter.width = .1, jitter.height = 0, dodge.width = 0.9), alpha=0.5)+ 
  scale_fill_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  scale_color_manual(values = c("#0072B5CC", "#BC3C29CC"))+
  theme_bw()+
  stat_compare_means(label = 'p.format', method = 'wilcox')+
  
  plot_layout(guides = 'collect', ncol = 2, nrow = 2)& theme(legend.position = 'bottom')
 
#Host HIV signal----
gene_hiv<-pathway.reactome[str_detect(names(pathway.reactome), 'HOST_INTERACTIONS_OF_HIV')]%>%unlist()%>%as.vector()

gene_hiv<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_hiv)

mono@meta.data <- mono@meta.data[, !grepl("HIV", colnames(mono@meta.data))]
mono<-AddModuleScore(mono, features=list(c(gene_hiv)), assay = 'SCT', name = 'HIV')


#Innate signal----
gene_inn<-pathway.reactome[str_detect(names(pathway.reactome), 'REACTOME_INNATE_IMMUNE_SYSTEM')]%>%unlist()%>%as.vector()

gene_inn<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_inn)

mono@meta.data <- mono@meta.data[, !grepl("Innate", colnames(mono@meta.data))]
mono<-AddModuleScore(mono, features=list(c(gene_inn)), assay = 'SCT', name = 'Innate')

#TNFA signal----
gene_tnf<-pathways.hallmark[str_detect(names(pathways.hallmark), 'TNFA')]%>%unlist()%>%as.vector()
names(gene_tnf)<-NULL

gene_tnf<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_tnf)

mono@meta.data <- mono@meta.data[, !grepl("TNF", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_tnf)), assay = 'SCT', name = 'TNFA')

FeaturePlot(mono, features = 'TNFA1', split.by = 'Group')

VlnPlot(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==FALSE), c('MTORC1', 'HIV1', 'TNFA1'), split.by = 'Group', pt.size = 0, flip=T, stack = T) + 
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#TNFA noncanonical signal----
gene_tnf2<-pathway.reactome[str_detect(names(pathway.reactome), 'TNFR2')]%>%unlist()%>%as.vector()
names(gene_tnf)<-NULL

gene_tnf2<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_tnf2)

mono@meta.data <- mono@meta.data[, !grepl("TNFA_NC", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_tnf2)), assay = 'SCT', name = 'TNFA_NC')

FeaturePlot(mono, features = 'TNFA_NC1', split.by = 'Group', order = T, max.cutoff = 10)

VlnPlot(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==FALSE), c('MTORC1', 'HIV1', 'TNFA1', 'TNFA_NC1'), split.by = 'Group', pt.size = 0, flip=T, stack = T) + 
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#Complement signal----
gene_complement<-pathway.reactome[str_detect(names(pathway.reactome), 'COMPLEMENT')]%>%unlist()%>%as.vector()
names(gene_complement)<-NULL

gene_complement<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_complement)

mono@meta.data <- mono@meta.data[, !grepl("Complement", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_complement)), assay = 'SCT', name = 'Complement')

FeaturePlot(mono, features = 'Complement1', split.by = 'Group', order = T, max.cutoff = 10)

VlnPlot(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==FALSE), c('MTORC1',  'TNFA1', 'Complement1'), split.by = 'Group', pt.size = 0, flip=T, stack = T) + 
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 1, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=5, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#IFNalpha signal----
gene_ifna<-pathways.hallmark[str_detect(names(pathways.hallmark), 'HALLMARK_INTERFERON_ALPHA_RESPONSE')]%>%unlist()%>%as.vector()
names(gene_ifna)<-NULL

gene_ifna<-base::intersect(mono@assays[["SCT"]]@data@Dimnames[[1]], gene_ifna)

mono@meta.data <- mono@meta.data[, !grepl("IFNA", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_ifna)), assay = 'SCT', name = 'IFNA_Response')

FeaturePlot(mono, features = 'IFNA_Response1', split.by = 'Group', order = T, max.cutoff = 10)

VlnPlot(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==FALSE), c('TNFA_NC1',  'HIV1', 'Innate1'), split.by = 'Group', pt.size = 0, flip=T, stack = T) + 
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 0.5, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=2, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#OXPHOS and Glycolysis signal----
gene_oxphos<-pathways.hallmark$HALLMARK_OXIDATIVE_PHOSPHORYLATION
gene_glyco<-pathways.hallmark$HALLMARK_GLYCOLYSIS

names(gene_oxphos)<-NULL
names(gene_glyco)<-NULL

mono@meta.data <- mono@meta.data[, !grepl("OXPHOS", colnames(mono@meta.data))]

mono<-AddModuleScore(mono, features=list(c(gene_oxphos)), assay = 'SCT', name = 'OXPHOS')
mono<-AddModuleScore(mono, features=list(c(gene_glyco)), assay = 'SCT', name = 'Glyco', search = T)

FeaturePlot(mono, features = 'OXPHOS1', split.by = 'Group', order = T,  keep.scale = 'all')
FeaturePlot(mono, features = 'Glyco1', split.by = 'Group', order = T,  keep.scale = 'all')

VlnPlot(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==FALSE), c('OXPHOS1', 'Glyco1',  'MTORC1', 'TNFA1'), split.by = 'Group', pt.size = 0, flip=T, stack = T) + 
  scale_fill_manual(values = c( "#BC3C29CC", "#0072B5CC")) +
  stat_summary(fun.y = function(x)quantile(x, probs =c(.25, .75)), geom='line', size = 0.5, colour = "grey20", position = position_dodge(0.9)) +
  stat_summary(fun.y=median, geom='point', size=2, shape=20, position = position_dodge(0.9)) + 
  stat_compare_means(label = "p.format", hide.ns = T)+
  theme(axis.text.x = element_text(angle = 0))+xlab('')+ylab('HIV Score')+
  scale_y_continuous(expand = expansion(.1, .1))

#Mtor, TNF signal and reservoir----
ggarrange(
  mono@meta.data%>%
    select(Group, Sample_Name, Clusters_Final, MTORC1, MMP2:CD4NADIR)%>%
    filter(Clusters_Final %in% c('pDC', 'cDC')==FALSE)%>%
    group_by(Group, Sample_Name, Clusters_Final)%>%
    reframe(MTORC1=mean(MTORC1), across(MMP2:CD4NADIR))%>%
    ungroup()%>%unique()%>%
    pivot_longer(cols=MMP2:CD4NADIR, names_to = 'Var', values_to = 'Value')%>%
    filter(!is.na(Value))%>%
    group_by(Clusters_Final, Var)%>%
    group_modify(.f=~broom::tidy(cor.test(.x$MTORC1, .x$Value, method = 'spearman', data = .x))) %>%
    ungroup()%>%
    mutate(padj=p.adjust(p.value, 'BH'))%>%
  arrange(padj)%>%mutate(text=ifelse(padj<0.25& p.value<0.05, Var, NA))%>%
  mutate(shape=ifelse(is.na(text), 'a', 'b'))%>%
  mutate(text=ifelse(str_detect(text, 'CRP|Intact|TNF|CD'), text, NA))%>%
  ggplot(aes(x=estimate, y=-log10(p.value)))+
  geom_point(aes(color=Clusters_Final, shape=shape))+ 
  scale_color_manual(values=pal_mono)+
  ggrepel::geom_text_repel(aes(label=text), size=4, min.segment.length = 0.05)+scale_shape_manual(values = c(46, 19))+
  theme_bw()+scale_y_continuous(breaks = c(0, 1, 2))  , 
  
  mono@meta.data%>%
    select(Group, Sample_Name, Clusters_Final, TNFA1, MMP2:CD4NADIR)%>%
    filter(Clusters_Final %in% c('pDC', 'cDC')==FALSE)%>%
    group_by(Group, Sample_Name, Clusters_Final)%>%
    reframe(TNFA1=mean(TNFA1), across(MMP2:CD4NADIR))%>%
    ungroup()%>%unique()%>%
    pivot_longer(cols=MMP2:CD4NADIR, names_to = 'Var', values_to = 'Value')%>%
    filter(!is.na(Value))%>%
    group_by(Clusters_Final, Var)%>%
    group_modify(.f=~broom::tidy(cor.test(.x$TNFA1, .x$Value, method = 'spearman', data = .x))) %>%
    ungroup()%>%
    mutate(padj=p.adjust(p.value, 'BH'))%>%
    arrange(padj)%>%mutate(text=ifelse(padj<0.25& p.value<0.05, Var, NA))%>%
    mutate(shape=ifelse(is.na(text), 'a', 'b'))%>%
    mutate(text=ifelse(str_detect(text, 'CRP|Intact|TNF|CD'), text, NA))%>%
    ggplot(aes(x=estimate, y=-log10(p.value)))+
    geom_point(aes(color=Clusters_Final, shape=shape))+ 
    scale_color_manual(values=pal_mono)+
    ggrepel::geom_text_repel(aes(label=text), size=4, min.segment.length = 0.05)+scale_shape_manual(values = c(46, 19))+
    theme_bw()  , 
  
  common.legend = T, ncol=1, legend = 'right')

mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, MTORC1, MMP2:CD4NADIR)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(MTORC1=mean(MTORC1), across(MMP2:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  filter(Clusters_Final=='NCM0')%>%
  ggplot(aes(x=MTORC1, y=log10(Intact+1)))+
  geom_point()+
  geom_smooth(method = 'lm', se=F)+
  theme_bw()+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.5, label.y.npc = 0.9)

mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, TNFA1, MMP2:CD4NADIR)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(TNFA1=mean(TNFA1), across(MMP2:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  filter(Clusters_Final=='NCM1')%>%
  ggplot(aes(x=TNFA1, y=log10(Intact+1)))+
  geom_point()+
  geom_smooth(method = 'lm', se=F)+
  theme_bw()+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.5, label.y.npc = 0.9)

mono@meta.data%>%
  select(Group, Sample_Name, Clusters_Final, Complement1, MMP2:CD4NADIR)%>%
  group_by(Group, Sample_Name, Clusters_Final)%>%
  reframe(Complement1=mean(Complement1), across(MMP2:CD4NADIR))%>%
  ungroup()%>%unique()%>%
  filter(Clusters_Final=='NCM1')%>%
  ggplot(aes(x=Complement1, y=log10(Intact+1)))+
  geom_point()+
  geom_smooth(method = 'lm', se=F)+
  theme_bw()+
  stat_cor(method = 'spearman', cor.coef.name = 'rho', label.x.npc = 0.5, label.y.npc = 0.9)

#Save----
saveRDS(mono, 'CVD_Monocyte_DC.rds')
save.image('Monocyte_DC.rdata')
saveRDS(subset(mono, Clusters_Final%in% c('pDC', 'cDC')==T), 'CVD_DC.rds')
 