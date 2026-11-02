library(tidyverse)
library(Seurat)
library(SeuratObject)
library(readxl)
library(openxlsx)
library(conflicted)
library(future)
library(SingleR)
library(Seurat.utils)
library(glmGamPoi)
library(DESeq2)
library(fgsea)
library(scplotter)

plan("multisession", workers = 6)
options(future.globals.maxSize = 50000 * 1024^3)

conflict_prefer('select', winner = 'dplyr')
conflict_prefer('filter', winner = 'dplyr')
conflicts_prefer(dplyr::count)

scdata<-readRDS('CVD_total.rds') 

metadata<-scdata@meta.data
idents<-unique(scdata@meta.data$Cell_Type_Main_Final)

pal<-c('B'='#D43F3A', 'CD8+ T'='#e37a37', 'CD14+ Monocytes'='#d0a63d', 
       'CD4+ T'='#79b354', 'NK'='#53b88e', "CD16+ Monocytes"='#46b8da',
       'Progenitor'='#3b95c8', "T"='#83bf97',
       "cDC"='#8241B9', 'NK and T LINGO2+'='#6F99ADCC',
       'Plasmablast'='#a367b8', "pDC"='#486DBE', 'Platelet'='grey') 

Idents(scdata)<-'Cell_Type_Main_Final'

DefaultAssay(scdata)<-'SCT' #Use SCT for violin plot
features<-c('CD3D', 'CD3E', 'CD4', 'CD8A', 'CD14', 'LYZ', 'FCGR3A', 'GNLY', 'MS4A1', 'CD1C', 'IL3RA', 
            'CD38', 'TNFRSF17', 'MKI67', 'CD34', 'PPBP', 'LINGO2')

VlnPlot(scdata, features = features, stack = T, flip = T)+
  guides(fill='none')+
  xlab('')

DimPlot(subset(scdata, seurat_clusters=='20'), group.by = 'seurat_clusters', raster=T, size=2, label = T)

FeaturePlot(scdata, c('CD4', 'GNLY', 'LINGO2', 'CD3D'), raster=T, pt.size = 3)

scplotter::CellDimPlot(scdata, theme = "theme_blank", 
                       palcolor =pal, label=T, label_size = 3, raster=T, raster_dpi=c(900, 900), pt_size=1.2) 

FeaturePlot(scdata, c('CD4'))

scplotter::CellDimPlot(scdata, theme = "theme_blank", split_by='orig.ident', legend.position = "none",
                       palcolor =pal, label=T, label_size = 3, raster=T, raster_dpi=c(900, 900), pt_size=1.5) 

scplotter::CellDimPlot(scdata, theme = "theme_blank", split_by='Group', legend.position = "none", ncol=1, 
                       palcolor =pal, label=T, label_size = 3, raster=T, raster_dpi=c(900, 900), pt_size=1.5) 

#Calculate proportion of cell types per sample and per group----
p1<-scCustomize::Proportion_Plot(seurat_object = subset(scdata, subset= Group=='High'), plot_type = "bar", colors_use = pal, 
                             split.by  = 'Sample_Name')
p2<-scCustomize::Proportion_Plot(seurat_object = subset(scdata, subset= Group=='Low'), plot_type = "bar", colors_use = pal, 
                                 split.by  = 'Sample_Name')

ggarrange(p1, p2, common.legend = T, legend = 'right', align = 'hv')

cellpercent<-scdata@meta.data%>%
  group_by(Group, Sample_Name, Cell_Type_Main_Final)%>% #Calculate percentage of cell type for each sample
  summarise(n=n())%>%
  ungroup()%>%
  group_by(Group, Sample_Name)%>%
  mutate(Percentage=n/sum(n))%>%
  ungroup()

cellpercent%>%
  ggplot(aes(x=Sample_Name, y=Percentage, fill=Cell_Type_Main_Final))+
  geom_bar(stat = 'identity', position = 'fill')+
  scale_fill_manual(values = pal)+
  theme_classic()+facet_grid(~Group)

cellpercent%>%
  group_by(Cell_Type_Main_Final)%>%
  group_modify(.f= ~as.data.frame(wilcox.test(Percentage~Group, data = .x)$p.value))

cellpercent%>%
  filter(Cell_Type_Main_Final %in% c('Plasmablast', 'Progenitor', 'Platelet', 'T', 'NK and T LINGO2+')==FALSE)%>%
  #mutate(Cell_Type_Main_Final=factor(Cell_Type_Main_Final, levels=c('CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells', 'cDC', 'pDC', 'Monocytes', 'CD16+ Non-classical Monocytes')))%>%
  mutate(Percentage=100*Percentage)%>%
  ggplot(aes(x=Group, y=Percentage))+
  geom_jitter(aes(color=Cell_Type_Main_Final), size = 2, alpha = 0.5, height = 0, width = 0.1)+
  geom_boxplot(alpha=0, linewidth=0.1)+
  facet_wrap(facets=~Cell_Type_Main_Final, nrow=4, scales='free')+
  scale_color_manual(values=pal)+theme_bw()+guides(color='none')+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1)))

cellpercent%>%
  filter(Cell_Type_Main_Final %in% c('Plasmablast', 'Progenitor', 'Platelet', 'T', 'NK and T LINGO2+')==FALSE)%>%
  #mutate(Cell_Type_Main_Final=factor(Cell_Type_Main_Final, levels=c('CD4+ T cells', 'CD8+ T cells', 'B cells', 'NK cells', 'cDC', 'pDC', 'Monocytes', 'CD16+ Non-classical Monocytes')))%>%
  mutate(Percentage=100*Percentage)%>%
  pivot_wider(id_cols = 'Sample_Name', names_from = c(Cell_Type_Main_Final, Group), values_from = Percentage)%>%
  select(contains('CD4+ T'), contains('CD8+ T'), contains('NK'), contains('T Prolif'), contains('B'), 
         contains('CD14+'), contains('CD16+'), contains('cDC'), contains('pDC'))%>% 
  write.xlsx('Fig3D.xlsx')
  
cvd<-read_xlsx('Grouping.xlsx', sheet='Main')%>%
  left_join(read_xlsx('Grouping.xlsx', sheet='Data')%>%select(MACSID, USLESION:STEN34), by=c('MACSID', 'US_Risk'))%>%
  right_join(cellpercent, by='Sample_Name')%>%
  select(Cell_Type_Main_Final, Sample_Name, Group, Percentage, USLESION:STEN34)

cor_res1<-cvd%>%
  select(Cell_Type_Main_Final, Sample_Name, Percentage, IMT, contains('SCORE'), contains('COUNT'))%>%
  pivot_longer(cols = c(IMT, contains('SCORE'), contains('COUNT')), 
               names_to = 'Var', values_to = 'Value')%>%
  filter(!is.na(Value))%>%
  group_by(Cell_Type_Main_Final, Var)%>%
  group_modify(.f=~broom::tidy(cor.test(.x$Percentage, .x$Value, method = 'spearman', data = .x))) %>%
  ungroup()%>%
  mutate(Padj=p.adjust(p.value, 'BH')) 

#Subsets to different cell types----
SaveSeuratRds(subset(scdata, subset= Cell_Type_Main_Final %in% c('CD14+ Monocytes', 'CD16+ Monocytes', 'pDC', 'cDC')),
              'CVD_Monocyte_DC.rds')

SaveSeuratRds(subset(scdata, subset= Cell_Type_Main_Final %in% c('CD8+ T', 'CD4+ T', 'NK', 'T', 'T Prolif', 'NK and T LINGO2+')),
              'CVD_T.rds')

SaveSeuratRds(subset(scdata, subset= Cell_Type_Main_Final=='B'|Cell_Type_Main_Final=='Plasmablast'), 'CVD_B_Plasma.rds')

#Save file----
save.image('CVD_YL_Downstream.RData')
