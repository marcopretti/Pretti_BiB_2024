# master script for Cross-reactivity analysis

pacman::p_load(tibble, RColorBrewer, Seurat, cowplot, corrplot, Hmisc, immunarch, ggseqlogo, qs, pheatmap)
library(hpar) # Human Protein Atlas info
set.seed(123)

setwd('~/Pretti2024')
load("results/workspace202302.RData")
### LOAD DATA----

### Select peptides to model, HLA-A2, 9mers
# source('bin/01a_public_peptides.R') # Bulek, Zhang, etc
# source('bin/01b_outgroup_peptides.R')
# source('bin/01_Peptide_selection/SupData1_selected_peptides.R')
peptide_cdr3_db <- readRDS('results/Supplementary_Data1_lista_annotated.rds') %>%
  dplyr::filter(nchar(Epitope)==9) # remove 10mers, remained from some of the studies and were modeled
peptide_cdr3_db = dplyr::filter(peptide_cdr3_db, origin %in% c('Viral','Human','synthetic'))
#saveRDS(peptide_cdr3_db, 'results/Supplementary_Data1_lista_annotated.rds')

### Read and process spatial, charge and radius information
source('bin/03_read_charge.R')

### Dimension reduction based on tabular summarized data
# peptide sequence matrix, encoded using Atchley factors;
# spatial and charge data using the entire pHLA, exposed pHLA residues, or peptide
source('bin/05_DimReduction.R') # depends on 03_.R

### pHLA images processing
# 1) Read images into RGB(Transparency) arrays and does R-B+T, RData_Final
source('bin/pHLA_mask/read_and_add_pHLA_images.R')
# 2) cbind all matrix into a dataframe and calculate SD regions, produces pHLA_masks_sd_mat_l
##source('bin/pHLA_mask/pHLA_to_array_and_SD.R')
pHLA_masks_sd_mat_l <- qread('results/RData/pHLA_masks_sd_mat_l_qsave.rds')
pHLA_masks_dist <- qread('results/RData/pHLA_masks_dist_simple_qsave.rds') # Images
pHLA_masks_diff_l_cbind_all <- qread('results/RData/pHLA_masks_diff_l_cbind_all.rds')
pHLA_masks_diff_l <- readRDS('results/RData/pHLA_masks_diff_l.rds')

# 3) Select relevant regions from image (manual choise + high SD) and perform umap
# source('bin/pHLA_mask/select_pHLA_masks_and_umap.R')
umap_pHLA_masks_l <- qread('results/RData/umap_pHLA_masks_l_qsave.rds')

#####

# Intro and Methods: Illustrations----
m_ <- dist(matrix(c(4,6,5,4,5,5,7,1,5,3,0),nrow=5))

pheatmap::pheatmap(m_, color = RColorBrewer::brewer.pal(9, 'Oranges'),
                   cluster_rows = F, cluster_cols = F, legend = FALSE)

pheatmap::pheatmap(m_, color = RColorBrewer::brewer.pal(9, 'Blues'),
                   cluster_rows = F, cluster_cols = F, legend = FALSE)

pheatmap::pheatmap(m_, color = RColorBrewer::brewer.pal(9, 'Greens'),
                   cluster_rows = F, cluster_cols = F, legend = FALSE)
#####

# Figuras e tabelas da introdução
# source('bin/Figuras_Tabelas.R')

# ANALYSIS
# 1)
# Statistics of the modeled epitopes, Table S1----
unique(peptide_cdr3_db$Epitope) %>% length()
unique(peptide_cdr3_db$CDR3_beta) %>% length()
unique(peptide_cdr3_db$CDR3_alpha) %>% length()

## Table S1
# Count peptides
peptide_cdr3_db %>%
  group_by(Epitope) %>% summarise(db=toString(sort(unique(db)))) %>%
  group_by(Reference=db) %>% summarise('Modeled Epitopes'=n()) %>%
  arrange(desc(`Modeled Epitopes`)) %>% formattable::formattable()

# Per origin: human, viral, synthetic
peptide_cdr3_db %>%
  dplyr::select(Epitope, origin) %>% unique() %>%
  group_by(origin) %>% summarise('Modeled Epitopes'=n()) %>% formattable::formattable()

# Cross-reactive Epitopes
unique(peptide_cdr3_db) %>%
  group_by(CDR3_beta) %>% dplyr::filter(n() > 1)
#####

# Supplementary Data 1, CDR3 from public databases recognizing multiple epitopes----
suppData <- readRDS('results/Supplementary_Data1_lista_annotated.rds')
suppData_ <- suppData %>%
  # select peptides in the public databases
  dplyr::filter(grepl('VDJdb|McPAS|TBAdb', db)) %>%
  group_by(CDR3_beta, CDR3_alpha) %>%
  summarise(n=n(),
            Epitope=toString(sort(unique(Epitope))), 
            origin=toString(sort(unique(origin)))) %>% na.omit()

suppData_ %>%
  dplyr::filter(grepl(',', origin)) %>% # among different 'species'
  dplyr::filter(n>1) %>% # select cross-reactive
  formattable::formattable()

#####

# Figure S3: pHLA mask with variable regions (SD) and selected regions----
# plot and save variable regions

# >> go to Figures_Manuscript
# cols = colorRampPalette(c('white','grey'))
# 
# # Save images of the pHLA marked for high SD
# lapply(names(pHLA_masks_sd_mat_l), function(x){
#   p=pheatmap::pheatmap(pHLA_masks_sd_mat_l[[x]], cluster_rows = FALSE, cluster_cols = FALSE, color = cols(5), 
#                        silent = TRUE, main = x)
# #  png(paste0('results/Figures/pHLA_masks_sd_', x, '.png'), 640*2, 480*2, res=150)
#   print(p)
#   dev.off()
# })
# 
# # Save images of the selected regions
# lapply(names(pHLA_masks_sd_mat_l), function(perspective_){
#   ind = ind[[perspective_]]
#   obj = pHLA_masks_sd_mat_l[[perspective_]]
#   obj[ind$rowInd,] = .5 # max 480, draw rows
#   obj[,ind$colInd] = .5 # max 640, draw columns
#   
#   cols = colorRampPalette(c('white','grey'))
#   p=pheatmap::pheatmap(obj, cluster_rows = FALSE, cluster_cols = FALSE, color = cols(5), 
#                        silent = TRUE, main = perspective_)
# #  png(paste0('results/Figures/pHLA_masks_sd_', perspective_, '_selected.png'), 640*2, 480*2, res=150)
#   print(p)
#   dev.off()
# })
#####


# Figure 1: Illustration in Pymol and Gimp


# Figure 2, Relationship between aa charge and pHLA image - Bulek2012----
# Heatmap of residues color / polarity / charge
# replaces _pmask_to_array.R


pHLA_masks_diff_l_cbind_bulek <- lapply(names(pHLA_masks_diff_l_cbind_all), function(perspective){
  # Transform array into vector to rbind into a df
  pHLA_masks_diff_l_array = lapply(pHLA_masks_diff_l[[perspective]][rownames(meta_bulek)], as.vector)
  print('Reducing')
  pHLA_masks_diff_l_cbind = Reduce("cbind", pHLA_masks_diff_l_array)
  colnames(pHLA_masks_diff_l_cbind) = names(pHLA_masks_diff_l_array)
  return(pHLA_masks_diff_l_cbind)
}); names(pHLA_masks_diff_l_cbind_bulek) = names(pHLA_masks_diff_l_cbind_all)

aa = c('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y')

# Get color of mutated epitope and compare with the reference PPI
bulek_heatmap_l <- lapply(names(pHLA_masks_diff_l_cbind_all), function(perspective_){
  l_ <- lapply(aa, function(pattern_){
    l <- lapply(c(1:9), function(pos_){
      # WT sequence (PPI)
      ref = pHLA_masks_diff_l_cbind_all[[perspective_]][,'ALWGPDPAA']
      # query: mutated PPI at given position
      query = `substr<-`('ALWGPDPAA',pos_,pos_,pattern_)
      
      if(query %in% colnames(pHLA_masks_diff_l_cbind_all[[perspective_]])){
        query = pHLA_masks_diff_l_cbind_all[[perspective_]][,query]
        query-ref
      }else{NULL}
    }); names(l) = paste0('P', c(1:9))
    return(l)
  }); names(l_) = aa
  l_
}); names(bulek_heatmap_l) = names(pHLA_masks_diff_l_cbind_all)

# Heatmap
pre_plot_bulek_ <- lapply(bulek_heatmap_l, function(perspective_) sapply(perspective_, function(aa_) sapply(aa_, mean)))

# >> go Figures_Manuscript.R
# cols = colorRampPalette(c('blue','white','red'))
# myBreaks <- c(seq(min(pre_plot_bulek_$Top, na.rm=TRUE), 0, length.out=ceiling(100/2) + 1), 
#               seq(max(pre_plot_bulek_$Top, na.rm=TRUE)/100, max(pre_plot_bulek_$Top, na.rm=TRUE), 
#                   length.out=floor(100/2)))
# 
# plots_l <- lapply(names(pre_plot_bulek_), function(x){
#   p=pheatmap::pheatmap(t(pre_plot_bulek_[[x]]), angle_col = 0, cluster_cols = FALSE, cluster_rows = FALSE, color = cols(100),
#                        breaks = myBreaks, main = x)
#   p[[4]]
# })
# plot_grid(plotlist = plots_l)
# 
# # single heatmap
# pheatmap::pheatmap(t(pre_plot_bulek_$Top), angle_col = 0, cluster_cols = FALSE, color = cols(100),
#                    breaks = myBreaks, cluster_rows = FALSE, fontsize = 14)
#####


# Table: Illustrate the problem of aggregating atoms into a single metric----
pqr_ex <- fread('models2input_CNN/GILGFVFTL/GILGFVFTL.pqr', fill=TRUE) %>% head(-2) 
colnames(pqr_ex) = c('Field_name','Atom_number','Atom_name','Residue_name','Residue_number','X', 'Y', 'Z','Charge','Radius')
pqr_ex %>% 
  dplyr::filter(Residue_number==1) %>% formattable::formattable()

pqr_gilgf <- dplyr::filter(pqr_df, Epitope=='GILGFVFTL')
dplyr::filter(pqr_gilgf, Molecule=='Epitope') %>% formattable::formattable()

# View selected atoms (Residue_number <=180) quadrante superior direito
# plot_ly(data = pqr_list[[1]], x=~Y, y=~Z, color=~Residue_number,
#         # Hover text:
#         text = ~paste(Residue_number))
#####

# FigSup - pHLA spatially to select the extracellular atoms----
cols = c('B2M'='red', 'HLA_alfa'='darkred', 'Epitope'='cyan3')
limitesZ = quantile(pqr_list[[1]]$Z)[c(1,5)]
limitesY = quantile(pqr_list[[1]]$Y)[c(1,5)]
tema = list(scale_color_manual(values = cols),
            scale_x_continuous(limits = limitesZ), scale_y_continuous(limits = limitesY),
            xlab('Z-coordinate in xyz space'), ylab('Y-coordinate in xyz space'))

pqr_gilgf = dplyr::filter(pqr_df, Epitope=='GILGFVFTL')

p1=pqr_gilgf %>%
  ggplot(aes(x=Z, y=Y, color=Molecule))+geom_point()+tema+
  geom_line(aes(group=Molecule), data = . %>% dplyr::filter(Molecule=='Epitope'))+
  guides(color=guide_legend(override.aes = list(linetype=0)))

p2=pqr_gilgf %>%
  dplyr::filter(Molecule!='B2M' & Residue_number<=180) %>%
  ggplot(aes(x=Z, y=Y, color=Molecule))+geom_point()+tema+guides(color='none')+
  geom_line(aes(group=Molecule), data = . %>% dplyr::filter(Molecule=='Epitope'))

p3=pqr_gilgf %>%
  dplyr::filter(Molecule=='Epitope') %>%
  ggplot(aes(x=Z, y=Y, color=Molecule, group=Molecule))+geom_point()+tema+guides(color='none')+geom_line()

# Figure paper
pqr_gilgf %>%
  ggplot(aes(x=Z, y=Y, color=Molecule))+geom_point()+
  geom_line(aes(group=Molecule), data = . %>% dplyr::filter(Molecule=='Epitope'))+
  guides(color=guide_legend(override.aes = list(linetype=0)))+
  scale_color_manual(values = c('B2M'='chartreuse2','HLA_alfa'='cyan3','Epitope'='red'))+
  scale_x_continuous(limits = limitesZ)+scale_y_continuous(limits = limitesY)+
  xlab('Z-coordinate')+ylab('Y-coordinate')+
  theme(axis.line = element_line(), panel.background = element_blank(), legend.title = element_blank(),
        axis.title = element_text(size=12, face='bold'), legend.text = element_text(size=11),
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_rect(fill="white"))

# plot_grid(p1, p2, p3, ncol=1, align = 'v', axis = 'lr')
# plot_grid(p1, p2, p3, nrow=1, rel_widths = c(1,.7,.7))

#####


# Figure S2A-B: Charge and spatial distribution across Position----

# >> go to Figures_Manuscript.R
# tema=list(theme(axis.text = element_text(size=11), axis.title = element_text(size=12), axis.line = element_line(),
#                 panel.background = element_blank()),scale_y_continuous(expand = c(0,NA)))
# 
# prePlot_charge_HuVi <- pep_l$epitope_charge[grepl('Charge', rownames(pep_l$epitope_charge)),rownames(meta_data_HuVi)]
# prePlot_charge_bulek <- pep_l$epitope_charge[grepl('Charge', rownames(pep_l$epitope_charge)),rownames(meta_bulek)]
# 
# p1_ <- prePlot_charge_bulek %>%
#   `rownames<-`(gsub('Epitope_(\\d)_Charge','P\\1',rownames(.))) %>% t() %>%
#   pheatmap(cluster_cols = F, color = colorRampPalette(brewer.pal(n = 6, name = "Greens"))(100), 
#            show_rownames = FALSE, border_color = NA, angle_col = 0)
# 
# p2_ <- prePlot_charge_HuVi %>%
#   `rownames<-`(gsub('Epitope_(\\d)_Charge','P\\1',rownames(.))) %>% t() %>%
#   pheatmap(cluster_cols = F, color = colorRampPalette(brewer.pal(n = 6, name = "Blues"))(100), 
#            show_rownames = FALSE, border_color = NA, angle_col = 0)
# 
# 
# prePlot_spatial <- pep_l$epitope_charge[!grepl('Charge', rownames(pep_l$epitope_charge)),] %>%
#   rownames_to_column('ID') %>%
#   reshape2::melt() %>% 
#   mutate(coord=gsub('Epitope_\\d_','',ID),
#          Pos=gsub('Epitope_(\\d)_\\S+', 'P\\1', ID)) 
# 
# 
# prePlot_spatial_bulek = dplyr::filter(prePlot_spatial, variable %in% rownames(meta_bulek)) %>%
#   # absolute positions are not meaningful when comparing different peptides
#   group_by(coord, Pos) %>% summarise(sd=sd(value), q1=quantile(value)[2], q3=quantile(value)[4], inter_q=q3-q1) %>%
#   group_by(Pos) %>% mutate(mean_sd=mean(sd), mean_inter_q=mean(inter_q))
# 
# prePlot_spatial_HuVi = dplyr::filter(prePlot_spatial, variable %in% rownames(meta_data_HuVi)) %>%
#   group_by(coord, Pos) %>% summarise(sd=sd(value), q1=quantile(value)[2], q3=quantile(value)[4], inter_q=q3-q1) %>%
#   group_by(Pos) %>% mutate(mean_sd=mean(sd), mean_inter_q=mean(inter_q)) 


# Figure paper

# p3_ <- prePlot_spatial_bulek %>%
#   ggplot(aes(x=Pos, fill=coord))+
#   geom_col(aes(y=sd+inter_q), position = 'dodge')+
#   geom_col(aes(y=inter_q), position = 'dodge', color='black')+
#   geom_line(aes(y=mean_inter_q+mean_sd, group=coord))+
#   ylab('Dispersion index\n(IQR+SD)')+
#   tema+theme(axis.title.x = element_blank())+guides(fill=guide_legend(title = 'coordinate'))
# 
# 
# p4_ <- prePlot_spatial_HuVi %>%
#   ggplot(aes(x=Pos, fill=coord))+
#   geom_col(aes(y=sd+inter_q), position = 'dodge')+
#   geom_col(aes(y=inter_q), position = 'dodge', color='black')+
#   geom_line(aes(y=mean_inter_q+mean_sd, group=coord))+
#   ylab('Dispersion index\n(IQR+SD)')+
#   tema+theme(axis.title.x = element_blank())+guides(fill='none')
# 
# plot_grid(p3_, p4_)
# 
# plot_grid(
#   p1_$gtable, NULL, p3_,
#   NULL, NULL, NULL,
#   p2_$gtable, NULL, p4_, nrow=3, rel_widths = c(1,.2,1), rel_heights = c(1,.1,1), align = 'hv', axis = 'lr'
# )

# backup
# prePlot_spatial_bulek %>%
#   ggplot(aes(x=Pos, y=inter_q, fill=coord))+geom_col(position = 'dodge')+#xlab('Residue position')+
#   geom_point(aes(y=mean_inter_q))+
#   geom_line(aes(y=mean_inter_q, group=coord))+
#   ylab('Interquartile range')+
#   tema+theme(axis.title.x = element_blank())+guides(fill='none')
# 
# prePlot_spatial_bulek %>%
#   ggplot(aes(x=Pos, y=sd, fill=coord))+geom_col(position = 'dodge')+#xlab('Residue position')+
#   geom_point(aes(y=mean_sd))+
#   geom_line(aes(y=mean_sd, group=coord))+
#   ylab('Interquartile range')+
#   tema+theme(axis.title.x = element_blank())+guides(fill='none')


prePlot_spatial_bulek %>%
  ggplot(aes(x=Pos, y=sd, fill=coord))+geom_col(position = 'dodge')+#xlab('Residue position')+
  geom_segment(aes(y=mean_, yend=mean_, xend=..x..+.1))+
  geom_line(aes(y=mean_, group=coord))+
  ylab('Standard deviation')+
  tema+theme(axis.title.x = element_blank())+guides(fill='none')


plot_grid(plot_grid(p1_2, p1_3, nrow=3), 
          plot_grid(p3_2, p3_3, nrow=3, align = 'v', axis = 'lr'), rel_widths = c(1,1.3))
#####

# Figure S1C-D: Same as above per residue (AA): Spatial----
prePlot_spatial_perResidue <- pep_l$epitope_charge[!grepl('Charge', rownames(pep_l$epitope_charge)),] %>%
  rownames_to_column('ID') %>%
  reshape2::melt() %>% 
  mutate(coord=gsub('Epitope_\\d_','',ID),
         Pos=as.numeric(gsub('Epitope_(\\d)_\\S+', '\\1', ID)),
         Residue=substr(variable, Pos, Pos),
         Pos=paste0('P',Pos)) %>%
  group_by(coord, Pos, Residue) %>% 
  summarise(sd=sd(value), q1=quantile(value)[2], q3=quantile(value)[4], inter_q=q3-q1)


#https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html#refs
aa_prop_imgt <- read.csv('data/Amino_acid_characteristics.csv') %>% as_tibble()
data.frame(Residue=gsub(' ','',aa_prop_imgt$Abbreviation2), Volume=aa_prop_imgt$Volume..A3.) %>%
  inner_join(prePlot_spatial_perResidue) %>%
  group_by(Residue, Volume, coord) %>% #summarise(sd=me)
  ggplot(aes(x=Volume, y=sd))+geom_point()+facet_grid(coord~Pos)+
  ggpubr::stat_cor(method = 'pearson', size=3.5)+
  ggtitle('Correlação da distribuição espacial com o volume dos aminoácidos nos epítopos modelados')+
  ylab('Desvio padrão das coordenadas espaciais')+xlab('Volume (Ångström³)')+tema

#####


# 2) individual scripts of Bulek and Lee


# 3)
# Statistics on Self and viral binding predictions----
# Self peptides (GRCh38)

## netMHCpan
# grep 'SB' GRCh38.pep.netMHCpan.txt > GRCh38.pep.netMHCpan.SB.txt
GRCh38_netMHC <- fread('data/netMHC_predictions_GRCh38/GRCh38.pep.netMHCpan.SB.txt.gz')
colnames(GRCh38_netMHC) = c('Pos', 'HLA', 'Peptide', 'Core', 'Of','Gp','Gl','Ip','Il','Icore','Identity', 'Score', 'Aff(nM)', '%Rank', 'nada','BindLevel')


# Viral proteins (NCBI)
## netMHCpan4
NCBIViral_netMHC <- fread('data/netMHC_predictions_NCBI_viral/Viral.pep.SB.netMHCpan.sed.txt.gz')
colnames(NCBIViral_netMHC) = c('Pos', 'HLA', 'Peptide', 'Core', 'Of','Gp','Gl','Ip','Il','Icore','Identity', 'Score', 'Aff(nM)', '%Rank', 'nada','BindLevel')

# See classes (?) of peptides
strsplit(NCBIViral_netMHC$Identity, '_') %>% lapply('[[', 1) %>% unlist() %>% table

# ANALYSIS
# - overlap between viral / human
# - overrepresentation of epitopes from each species (n=2,3,etc)
# - from which species they come more often
# - overrepresentation disconsidering P2 and P9

filter_ = function(x) dplyr::filter(x, !grepl('X', Peptide)) %>% dplyr::select(-c(HLA, nada)) %>% mutate(Pos=as.numeric(Pos))
HuVi_netMHC = bind_rows(filter_(GRCh38_netMHC) %>% mutate(source='Human'),
                        NCBIViral_netMHC %>% 
                          mutate(source='Viral', Pos=gsub('.*:','',Pos)) %>% filter_())

HuVi_netMHC_summ = HuVi_netMHC %>%
  mutate(Pep_len=nchar(Peptide)) %>%
  group_by(Peptide, Pep_len, source) %>% summarise(n=n())

# Overrepresented epitopes within species
length(unique(HuVi_netMHC_summ$Peptide)) # unique epitopes
table(HuVi_netMHC_summ$source) # how many from each organism

HuVi_netMHC_summ %>%
  group_by(source, n) %>% 
  summarise(nn=n())

# >> go to Figures_Manuscript.R
# tema_ = list(theme(axis.title = element_text(size=12), axis.text = element_text(size=11), 
#                    legend.text = element_text(size=11), legend.title = element_text(size=12)))
# 
# HuVi_netMHC_summ %>%
#   mutate(n=ifelse(n>10,10,n)) %>%
#   group_by(Pep_len, source, n) %>% summarise(overRep=n()) %>%
#   ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
#   facet_grid(~Pep_len)+scale_y_log10(labels = scales::label_comma())+scale_x_continuous(breaks = 1:10)+
#   ggtitle('Overrepresentation of human and viral epitopes')+
#   ylab('Epitopes')+xlab('Overrepresentation')+
#   theme_bw()+tema_

# Human / Viral shared epitopes
shared_pep_HuVi_ = HuVi_netMHC_summ %>% ungroup() %>%
  group_by(Peptide) %>% dplyr::filter(n()>1) %>% pull(Peptide)

# annotate shared peptides
library(AnnotationDbi)
library(EnsDb.Hsapiens.v86)

HuVi_netMHC_shared_pep = dplyr::filter(HuVi_netMHC, Peptide %in% shared_pep_HuVi_)

HuVi_netMHC_shared_pep = HuVi_netMHC_shared_pep %>%
  mutate(SYMBOL=ifelse(grepl('^ENSP', Identity), mapIds(EnsDb.Hsapiens.v86, keys = Identity, column = 'SYMBOL', keytype = 'PROTEINID'),
                       gsub('([AP|NP|WP|YP]_\\d+)_(\\d+)_.*','\\1.\\2', Identity))) 

# Get protein name from fasta file
ViralID = system("zgrep '^>' data/Viral_Proteins_NCBI/Viruses_Protein_clusters_NCBI_NO_HYPOTHETICAL_PTN.fasta.gz", intern = TRUE)
substr(ViralID, 2,3) %>% table # see what are the prefixes in the viral ID [AP, NP, WP, YP]
ViralID_ = data.frame(fastaHeader=ViralID,
                      ViralID=gsub('>(\\S+) .*','\\1',ViralID))

HuVi_netMHC_shared_pep$ViralProtein = ViralID_$fastaHeader[match(HuVi_netMHC_shared_pep$SYMBOL, ViralID_$ViralID)]
HuVi_netMHC_shared_pep$ViralProtein = gsub('>\\S+ ','',HuVi_netMHC_shared_pep$ViralProtein)
HuVi_netMHC_shared_pep$ViralSpecies = gsub('.*\\[(.*)\\]','\\1',HuVi_netMHC_shared_pep$ViralProtein)
HuVi_netMHC_shared_pep$ViralProtein = gsub(' \\[.*\\]','',HuVi_netMHC_shared_pep$ViralProtein)

HuVi_netMHC_shared_pep = HuVi_netMHC_shared_pep %>%
  mutate(HumanProtein=ifelse(grepl('^ENSP', Identity), Identity, '')) %>%
  group_by(Peptide, Score, `Aff(nM)`, `%Rank`, source) %>% 
  summarise(SYMBOL=toString(sort(unique(SYMBOL))), HumanProtein=toString(sort(unique(HumanProtein))),
            ViralSpecies=toString(sort(unique(ViralSpecies))), ViralProtein=toString(sort(unique(ViralProtein)))) %>%
  mutate(HumanSYMBOL=ifelse(source=='Human', SYMBOL, ''),
         ViralNCBI_ID=ifelse(source=='Viral', SYMBOL, ''), SYMBOL=NULL, source=NULL) %>%
  summarise(HumanProtein=toString(HumanProtein), ViralProtein=toString(ViralProtein), ViralSpecies=toString(ViralSpecies),
            HumanSYMBOL=toString(HumanSYMBOL), ViralNCBI_ID=toString(ViralNCBI_ID)) %>%
  mutate(across(.fns = function(x) gsub('^, ','',x))) %>%
  mutate(across(.fns = function(x) gsub(', $','',x)))

# Supplementary Table: shared peptides
# HuVi_netMHC_shared_pep %>% 
#   rename_with(.cols = c(Score, `Aff(nM)`, `%Rank`), .fn = function(x) gsub('^','netMHCpan4_',x)) %>%
#   write.csv(file = 'results/Supplementary_Data_2.csv', row.names = FALSE)


# Table of percentage of shared peptides
HuVi_netMHC_shared = HuVi_netMHC_summ %>%
  mutate(carac=ifelse(Peptide %in% shared_pep_HuVi_, 'shared', 'restricted'))

HuVi_netMHC_shared %>% group_by(Peptide, Pep_len, carac) %>% summarise() %>%
  group_by(Pep_len, carac) %>% summarise(n=n()) %>%
  reshape2::dcast(Pep_len~carac) %>%
  mutate(perc_shared=signif((shared/(restricted+shared)*100),2)) %>% 
  formattable::formattable()

# between species
# >> go to Figures_Manuscript.R
# plot_ = HuVi_netMHC_shared %>%
#   mutate(overRep=ifelse(n>10,10,n)) %>%
#   group_by(overRep, Pep_len, carac) %>% summarise(n=n()) %>%
#   group_by(Pep_len, carac) %>% mutate(prop=prop.table(n)) 
# 
# plot_ %>%
#   ggplot(aes(x=overRep, y=prop, color=carac))+
#   geom_point()+geom_line()+
#   scale_x_continuous(breaks = 1:10)+facet_grid(~Pep_len)+
#   ggtitle('Proportion of shared epitopes between human and virus')+
#   ylab('Proportion of Epitopes')+xlab('Overrepresentation')+
#   scale_color_manual(values = c('black','gray'))+theme_bw()+tema_


# Shared epitopes desconsidering P2 and P9
## within species
# >> go to Figures_Manuscript.R
# HuVi_netMHC_summ_tr = HuVi_netMHC %>%
#   dplyr::filter(nchar(Peptide) %in% c(9,10)) %>%
#   mutate(Peptide_trimm=gsub('(\\S)\\S(.*)\\S','\\1\\2',Peptide)) %>%
#   group_by(Peptide_trimm, source) %>% summarise(n=n())
# 
# HuVi_netMHC_summ_tr %>%
#   mutate(n=ifelse(n>10,10,n),
#          Pep_len=nchar(Peptide_trimm)) %>%
#   group_by(source, n, Pep_len) %>% summarise(overRep=n()) %>%
#   ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
#   scale_y_log10(labels = scales::label_comma())+scale_x_continuous(breaks = 1:10)+
#   ggtitle('Overrepresentation of human and viral trimmed epitopes')+
#   facet_grid(~Pep_len)+ylab('Epitopes')+xlab('Overrepresentation')+theme_bw()+tema_


# between species
shared_pep_HuVi_tr_ = HuVi_netMHC_summ_tr %>% ungroup() %>%
  group_by(Peptide_trimm) %>% dplyr::filter(n()>1) %>% pull(Peptide_trimm)

HuVi_netMHC_summ_tr_shared = HuVi_netMHC_summ_tr %>%
  mutate(carac=ifelse(Peptide_trimm %in% shared_pep_HuVi_tr_, 'shared', 'restricted'),
         Pep_len=nchar(Peptide_trimm))

HuVi_netMHC_summ_tr_shared %>% group_by(Peptide_trimm, Pep_len, carac) %>% summarise() %>%
  group_by(Pep_len, carac) %>% summarise(n=n()) %>%
  reshape2::dcast(Pep_len~carac) %>%
  mutate(perc_shared=signif((shared/(restricted+shared)*100),2)) %>% 
  formattable::formattable()


# Are the shared epitopes overrepresented many times ?
# >> go to Figures_Manuscript.R
# plot_ = HuVi_netMHC_summ_tr_shared %>%
#   mutate(overRep=ifelse(n>10,10,n),
#          Pep_len=nchar(Peptide_trimm)) %>%
#   group_by(overRep, Pep_len, carac) %>% summarise(n=n()) %>%
#   group_by(Pep_len, carac) %>% mutate(prop=prop.table(n)) 
# 
# plot_ %>%
#   ggplot(aes(x=overRep, y=prop, color=carac))+
#   geom_point()+geom_line()+
#   scale_x_continuous(breaks = 1:10)+facet_grid(~Pep_len)+
#   ggtitle('Proportion of shared trimmed epitopes between human and virus')+
#   ylab('Proportion of Epitopes')+xlab('Overrepresentation')+
#   scale_color_manual(values = c('black','gray'))+theme_bw()+tema_


# From where do the shared epitopes come from
library(AnnotationDbi)
detach('IRanges')
.libPaths('/scr/R/Rpackages/4.3.0/')
library(EnsDb.Hsapiens.v86)

HuVi_netMHC_annotated = readRDS('results/RData/HuVi_netMHC_annotated.rds')

mutate(mapIds(EnsDb.Hsapiens.v86, keys = Identity, keytype = 'SYMBOL', column = 'GENEID'))

#
Viruses_Protein_clusters = Biostrings::readAAStringSet('data/Viral_Proteins_NCBI/Viruses_Protein_clusters_NCBI_NO_HYPOTHETICAL_PTN.fasta')
Viruses_Protein_clusters = data.frame(ID=names(Viruses_Protein_clusters),
                                      ProteinID=gsub(' .*','',names(Viruses_Protein_clusters)),
                                      Description=gsub('\\S+\\.1 (.*) \\[.*','\\1',names(Viruses_Protein_clusters)),
                                      Virus_Species=gsub('.*\\[(.*)\\]','\\1',names(Viruses_Protein_clusters))) %>% as_tibble()

unique(Viruses_Protein_clusters$Species_high)

shared_viral_species = dplyr::filter(HuVi_netMHC, Peptide %in% shared_pep_HuVi_ & !grepl('^ENSP', Identity)) %>% pull(Identity) %>%
  gsub('(\\S+)_(\\S+)_1.*','\\1_\\2.1',.) %>% unique()

Viruses_Protein_clusters_shared = dplyr::filter(Viruses_Protein_clusters, ProteinID %in% shared_viral_species)
table(Viruses_Protein_clusters_shared$Species_high) %>% sort()

# Among the viral epitopes, do the shared epitopes associate with some species ?
# Species less likely to be in contact with human / or to infect humans have more epitopes shared with human ptns?

del = inner_join(
  Viruses_Protein_clusters_shared %>%
    group_by(Virus_Species) %>% summarise(n_shared=n()),
  Viruses_Protein_clusters %>%
    group_by(Virus_Species) %>% summarise(n_total=n())
)

#####


# Figure 5: Cross-reactivity between viral and human epitopes, waterfall----
dist_selfViral_pHLA_filt = dist(t(pep_trimm_l$pHLA_charge_filt[,rownames(meta_data_HuVi)]), method = 'euclidean') %>% as.matrix()
dist_selfViral_epitope = dist(t(pep_trimm_l$epitope_charge[,rownames(meta_data_HuVi)]), method = 'euclidean') %>% as.matrix()
dist_selfViral_Images = lapply(pHLA_masks_dist$canberra, as.matrix)
dist_selfViral_Images = lapply(dist_selfViral_Images, function(x) x[rownames(meta_data_HuVi),rownames(meta_data_HuVi)])

melt_dist_selfViral = function(dist_selfViral_){
  diag(dist_selfViral_) = NA
  
  dist_selfViral_ = reshape2::melt(dist_selfViral_) %>%
    na.omit() %>% 
    arrange(value) %>%
### wrong !    dplyr::filter(!duplicated(value)) %>% # nao adequado
    mutate(ID=paste0(Var1,'_',Var2),  # Create unique ID
           Var1_origin=ifelse(Var1 %in% rownames(meta_data)[meta_data$origin=='Human'], 'Human', 'Viral'),
           Var2_origin=ifelse(Var2 %in% rownames(meta_data)[meta_data$origin=='Human'], 'Human', 'Viral'),
           selfViral=ifelse(Var1_origin!=Var2_origin, 'Human/Viral',
                            ifelse(Var1_origin=='Viral', 'Viral/Viral', 'Human/Human')))
  
  table(dist_selfViral_$selfViral) %>% prop.table()*100
  dist_selfViral_ %>% group_by(selfViral) %>% summarise(median(value))
  return(dist_selfViral_)
}

dist_selfViral_pHLA_filt_melted = melt_dist_selfViral(dist_selfViral_pHLA_filt)
dist_selfViral_epitope_melted = melt_dist_selfViral(dist_selfViral_epitope)
dist_selfViral_Images_melted = lapply(dist_selfViral_Images, melt_dist_selfViral)


# >> go to Figures_Manuscript.R
# plot_dist_selfViral_melted = function(dist_selfViral_melted_, subtitle_=NULL){
#   dist_selfViral_melted_ %>% 
#     ggplot(aes(x=selfViral, y=value))+geom_violin()+geom_boxplot()+
#     ggpubr::stat_compare_means()+
#     #    ggtitle('Similarity between human and viral', subtitle_)+
#     ggtitle(subtitle_)+
#     ylab(paste0('Distance between ', subtitle_, '\n(method=camberra)'))+xlab('')+theme_bw()+
#     theme(axis.text = element_text(size=11), axis.title = element_text(size=12))
# }
# p1=plot_dist_selfViral_melted(dist_selfViral_Images_melted$Top, 'pHLA Top image')
# p2=plot_dist_selfViral_melted(dist_selfViral_epitope_melted, 'epitopes')
# p3=plot_dist_selfViral_melted(dist_selfViral_pHLA_filt_melted, 'ec pHLA')
# 
# plot_grid(p1, p2, p3, nrow=1)

# Stats for boxplot
dist_selfViral_Images_melted$Top %>% group_by(selfViral) %>% summarise(median(value))
dist_selfViral_pHLA_filt_melted %>% group_by(selfViral) %>% summarise(median(value))
dist_selfViral_epitope_melted %>% group_by(selfViral) %>% summarise(median(value))

# function to Waterfall dotplot
wat_dotPlot <- function(dist_HumVir_melted_, ggtitle_='Similarity between human and viral epitopes', ylab_=NULL){
  require('ggbreak')  
  
  tmp_ = bind_rows(dist_HumVir_melted_ %>% ungroup() %>%
              # get the top 300, distant epitopes between them
              arrange(value) %>% top_n(300, wt = -value),
            dist_HumVir_melted_ %>% ungroup() %>%
              # get the 300, closer epitopes between them
              arrange(value) %>% top_n(300, wt = value))
  
  y_l = sort(tmp_$value) %>% head(300) %>% tail(1)*1.01
  y_u = sort(tmp_$value) %>% tail(300) %>% head(1)*0.99
  print(c(y_l, y_u))
  
  tmp_ %>%
    ggplot(aes(x=reorder(ID, value), y=value, color=selfViral))+geom_point(size=.7)+
    theme(axis.text.x = element_blank(), panel.background = element_blank(),
          axis.text.y = element_text(size=11),
          axis.line = element_line(), axis.title = element_text(size=12), 
          legend.text = element_text(size=11), plot.title = element_text(hjust=.5, size=12))+
    ylab(paste0('Distance between ',ylab_,'\n(method=camberra)'))+
    ggtitle(ggtitle_)+xlab('Epitope pair')+
    scale_color_manual(values = c('black','red','gray'))+
    scale_y_break(c(y_l,y_u))
}

# >> go to Figures_Manuscript.R
# wat_dotPlot(dist_selfViral_Images_melted$rightTop, 'pHLA rightTop image', 'images')+NoLegend()
# wat_dotPlot(dist_selfViral_Images_melted$leftTop, 'pHLA leftTop image')+NoLegend()
# wat_dotPlot(dist_selfViral_Images_melted$downTop, 'pHLA downTop image')+NoLegend()
# wat_dotPlot(dist_selfViral_Images_melted$upTop, 'pHLA upTop image')+NoLegend()

# p4=wat_dotPlot(dist_selfViral_Images_melted$Top, 'pHLA Top perspective', 'images')+NoLegend()
# p5=wat_dotPlot(dist_selfViral_pHLA_filt_melted, 'ec pHLA', 'pHLA')+NoLegend()
# p6=wat_dotPlot(dist_selfViral_epitope_melted, 'Epitope', 'epitopes')+theme(axis.title.x = element_text(hjust=0.3))
# 
# aplot::plot_list(p4, p5, p6, nrow = 1, widths = c(1,1,1.4))

#####



# Figure 6: CR between viral and cancer epitopes----
# Select epitopes from tumor or tumor associates and get extra info in McPAS
mcpas_sep = fread("data/McPAS/McPAS-TCR_10_Sept_2022.csv") 
mcpas_sep = mcpas_sep %>%
  dplyr::filter(Epitope.peptide %in% peptide_cdr3_db$Epitope)

# Select epitopes from tumor or tumor associates
mcpas_cancer = mcpas_sep %>%
  dplyr::filter(grepl('leukemia|Lymphoma|Melanoma|Tumor|carcinoma|Neoantigen', Pathology)) %>% # Neoantigen (GLYDGMEHL)
  group_by(Epitope=Epitope.peptide) %>% summarise(Pathology=toString(sort(unique(Pathology))))

# ALTPVVVTL, AML
# fvyvfIthl (MUT)	fvyvfTthl (WT), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4993154/#SD1 (Sup6)

# Check distance between them and human/viral epitopes, likelyhood of CR
CancerViralCR_Top = dplyr::rename(dist_selfViral_Images_melted$Top, 'Epitope'=Var1) %>% inner_join(mcpas_cancer)
CancerViralCR_rightTop = dplyr::rename(dist_selfViral_Images_melted$rightTop, 'Epitope'=Var1) %>% inner_join(mcpas_cancer)

CancerViralCR_ecpHLA = dplyr::rename(dist_selfViral_pHLA_filt_melted, 'Epitope'=Var1) %>% inner_join(mcpas_cancer)
CancerViralCR_epitope = dplyr::rename(dist_selfViral_epitope_melted, 'Epitope'=Var1) %>% inner_join(mcpas_cancer)

process_CancerViralCR = function(CancerViralCR_){
  # Find the closer epitope from Cancer epitopes
  CancerViralCR_ = CancerViralCR_ %>%
    mutate(Pathology=gsub('Acute myeloid leukemia','AML', gsub('Tumor associated antigen \\(TAA\\)','TAA',Pathology))) %>%
    group_by(Epitope, selfViral, Pathology) %>%
    summarise(value_min=min(value))
  
  # Remove epitopes found in human and virus (same sequence)
  ind_ = dplyr::filter(CancerViralCR_, selfViral=='Viral/Viral') %>% pull(Epitope) %>% unique()
  CancerViralCR_ %>%
    dplyr::filter(!Epitope %in% ind_)
}

# >> go to Figures_Manuscript.R
plot_CR_dotplot = function(CancerViralCR_, dist_selfViral_, ggtitle_=NULL){
  HuVi_ = dplyr::filter(dist_selfViral_, selfViral=='Human/Viral')
  HuHu_ = dplyr::filter(dist_selfViral_, selfViral=='Human/Human')
  
  limits_ = quantile(CancerViralCR_$value_min)[c(1,5)]*c(.9,1.1)
  width_ = max(limits_)/50
  
  CancerViralCR_ %>% 
    reshape2::dcast(Epitope+Pathology~selfViral) %>%
    mutate(Pathology=factor(Pathology, levels=c('Neoantigen','AML','Lymphoma','Melanoma, Neoantigen','TAA'))) %>% arrange(Pathology) %>%
    ggplot()+geom_point(aes(x=`Human/Human`, y=`Human/Viral`, color=Pathology))+geom_abline()+
    geom_boxplot(data = HuVi_, aes(x=value, y=min(CancerViralCR_$value_min)), outlier.shape = NA, width=width_)+
    geom_boxplot(data = HuHu_, aes(y=value, x=min(CancerViralCR_$value_min)), outlier.shape = NA, width=width_)+
    scale_y_continuous(limits = limits_, labels = scales::comma)+
    scale_x_continuous(limits = limits_, labels = scales::comma)+
    ggtitle(ggtitle_)+
    ylab('Distance to the closest viral epitope')+xlab('Distance to the closest self epitope')+
    theme_bw()+theme(axis.title = element_text(size=12), axis.text = element_text(size=11), plot.title = element_text(hjust = .5),
                     legend.text = element_text(size=10), legend.position = 'bottom', plot.margin = margin(5,15,5,5))+
    guides(color=guide_legend(title = element_blank()))+
    scale_color_manual(values = c('Neoantigen'='gray','AML'='red','Lymphoma'='purple','Melanoma, Neoantigen'='darkgreen','TAA'='orange'))
}

# Figure 6A-C
CancerViralCR_Top_ = process_CancerViralCR(CancerViralCR_Top)
# plot_CR_dotplot(CancerViralCR_Top_, dist_selfViral_ = dist_selfViral_Images_melted$Top, ggtitle_ = 'pHLA Top image')

# Stats
unique(CancerViralCR_Top_$Epitope) %>% length()

# To search for the pHLAs
ind_ = CancerViralCR_Top_ %>% ungroup() %>% dplyr::filter(Pathology=='TAA') %>% dplyr::select(Var1=Epitope, value=value_min)
dist_selfViral_Images_melted$Top %>% inner_join(ind_)

# IMNDMPIYM ALLDWVTSV (Human, 1347) AMFDLIYPI (Viral, 1424)
# FLCMKALLL FLAHVLNPV (Viral, 1761) FLACHLFVV (Human, 1813)

dplyr::filter(mcpas_sep, Epitope.peptide=='IMNDMPIYM') %>% group_by(Epitope.peptide,Antigen.protein) %>% summarise
dplyr::filter(mcpas_sep, Epitope.peptide=='FLCMKALLL') %>% group_by(Epitope.peptide,Antigen.protein) %>% summarise


# Figure 6B
# CancerViralCR_rightTop_ = process_CancerViralCR(CancerViralCR_rightTop)
# plot_CR_dotplot(CancerViralCR_rightTop_, dist_selfViral_ = dist_selfViral_Images_melted$rightTop, ggtitle_ = 'pHLA rightTop image')

# CancerViralCR_ecpHLA_ = process_CancerViralCR(CancerViralCR_ecpHLA)
# plot_CR_dotplot(CancerViralCR_ecpHLA_, dist_selfViral_ = dist_selfViral_pHLA_filt_melted, ggtitle_ = 'ec pHLA')

# To search for the pHLAs
ind_ = CancerViralCR_rightTop_ %>% ungroup() %>% dplyr::filter(Pathology=='TAA') %>% dplyr::select(Var1=Epitope, value=value_min)
dist_selfViral_Images_melted$rightTop %>% inner_join(ind_)
# blastp

ind_ = CancerViralCR_ecpHLA_ %>% ungroup() %>% dplyr::filter(Pathology=='TAA') %>% dplyr::select(Var1=Epitope, value=value_min)
dist_selfViral_pHLA_filt_melted %>% inner_join(ind_)


# Figure 6C
# CancerViralCR_epitope_ = process_CancerViralCR(CancerViralCR_epitope)
# plot_CR_dotplot(CancerViralCR_epitope_, dist_selfViral_ = dist_selfViral_epitope_melted, ggtitle_ = 'Epitope')

# To search for the pHLAs
ind_ = CancerViralCR_epitope_ %>% ungroup() %>% dplyr::filter(Pathology=='TAA') %>% dplyr::select(Var1=Epitope, value=value_min)
dist_selfViral_epitope_melted %>% inner_join(ind_) %>% mutate(value=signif(value, 3))


#####


# Figure 7: CR between T1D (autoimmune), viral and self epitopes----
head(normtissue <- hpaNormalTissue())

mcpas_autoimmune = mcpas_sep %>%
  dplyr::filter(Category=='Autoimmune' | grepl('rthritis', Additional.study.details)) %>%
  mutate(Additional.study.details=ifelse(is.na(Additional.study.details), Category, Additional.study.details),
         Pathology=ifelse(grepl('rthritis', Additional.study.details), 'Rheumatoid Arthritis', Pathology)) %>%
  group_by(Epitope=Epitope.peptide) %>% summarise(Pathology=toString(sort(unique(Pathology))))

# Restrain to T1D epitope
mcpas_autoimmune = dplyr::filter(mcpas_autoimmune, Epitope=='VLFGLGFAI')

AutoimmuneViralCR = dplyr::rename(dist_selfViral_Images_melted$Top, 'Epitope'=Var1) %>% inner_join(mcpas_autoimmune) %>%
  dplyr::filter(selfViral=='Human/Viral' & Var2 %in% peptide_cdr3_db$Epitope[peptide_cdr3_db$db=='NCBIViral']) %>%
  dplyr::rename(Peptide=Var2)

AutoimmuneSelfCR = dplyr::rename(dist_selfViral_Images_melted$Top, 'Epitope'=Var1) %>% inner_join(mcpas_autoimmune) %>%
  dplyr::filter(selfViral=='Human/Human' & Var2 %in% peptide_cdr3_db$Epitope[peptide_cdr3_db$db=='GRCh38'])

# Find the closer epitope from AID epitopes
# Diabetis
AutoimmuneSelfCR = AutoimmuneSelfCR %>%
  dplyr::rename(Peptide=Var2) %>% inner_join(HuVi_netMHC) %>%
  mutate(`Gene.name`=mapIds(EnsDb.Hsapiens.v86, Identity, 'SYMBOL', 'PROTEINID'), Identity=NULL) %>% 
  group_by(Epitope, Peptide, value, ID, Var1_origin, Var2_origin, selfViral, Pathology, `Aff(nM)`, `%Rank`, `Gene.name`)

AutoimmuneSelfCR = left_join(AutoimmuneSelfCR, normtissue)

# << go to Figures_Manuscript.R
# wat_dotPlot_self <- function(AutoimmuneViralCR_, ggtitle_='Similarity between VLFGLGFAI and self epitopes', ylab_=NULL){
#   require('ggbreak')  
#   
#   AutoimmuneViralCR_ = AutoimmuneViralCR_ %>% na.omit() %>%
#     dplyr::filter(Level=='High') %>%
#     mutate(Cell.type=ifelse(Cell.type=='pancreatic endocrine cells', 'Pancreatic cells', 'Other')) %>%
#     group_by(Cell.type, Peptide, ID, Level, value) %>% summarise() %>% ungroup()
#   
#   tmp_ = bind_rows(AutoimmuneViralCR_ %>% ungroup() %>%
#                      # get the top 100, distant epitopes between them
#                      arrange(value) %>% top_n(100, wt = -value),
#                    AutoimmuneViralCR_ %>% ungroup() %>%
#                      # get the 100, closer epitopes between them
#                      arrange(value) %>% top_n(100, wt = value))
#   
#   y_l = sort(tmp_$value) %>% head(100) %>% tail(1)*1.01
#   y_u = sort(tmp_$value) %>% tail(100) %>% head(1)*0.99
#   print(c(y_l, y_u))
#   
#   tmp_ %>%
#     ggplot(aes(x=reorder(ID, value), y=value, color=Cell.type))+
#     geom_point(size=2, alpha=.5, aes(shape=Cell.type))+
#     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
#           panel.background = element_blank(), plot.title = element_text(hjust=.5, size=12),
#           axis.text.y = element_text(size=11), axis.line = element_line(), axis.title = element_text(size=12), 
#           legend.text = element_text(size=11), legend.key = element_rect(fill = 'white'), legend.position = 'bottom')+
#     scale_shape_manual(values = c(1,3))+
#     ylab(paste0('Distance ',ylab_,'\n(method=canberra)'))+ggtitle(ggtitle_)+xlab('Epitopes')+
#     scale_color_manual(values = c('Pancreatic cells'='red',Other='black'))+scale_x_discrete(expand = expansion(mult = c(.01,.01)))+
#     scale_y_break(c(y_l,y_u))+scale_y_continuous(labels = scales::comma)+
#     theme(axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank())
# }
# 
# wat_dotPlot_viral <- function(AutoimmuneViralCR_, ggtitle_='Similarity between VLFGLGFAI and viral epitopes', ylab_=NULL){
#   require('ggbreak')  
#   
#   AutoimmuneViralCR_ = AutoimmuneViralCR_ %>%
#     group_by(Peptide, ID, value) %>% summarise() %>% ungroup()
#   
#   tmp_ = bind_rows(AutoimmuneViralCR_ %>% ungroup() %>%
#                      # get the top 100, distant epitopes between them
#                      arrange(value) %>% top_n(100, wt = -value),
#                    AutoimmuneViralCR_ %>% ungroup() %>%
#                      # get the 100, closer epitopes between them
#                      arrange(value) %>% top_n(100, wt = value))
#   
#   y_l = sort(tmp_$value) %>% head(100) %>% tail(1)*1.01
#   y_u = sort(tmp_$value) %>% tail(100) %>% head(1)*0.99
#   print(c(y_l, y_u))
#   
#   tmp_ %>%
#     ggplot(aes(x=reorder(ID, value), y=value))+
#     geom_point(size=2, fill='white', shape=1)+
#     theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
#           panel.background = element_blank(), plot.title = element_text(hjust=.5, size=12),
#           axis.text.y = element_text(size=11), axis.line = element_line(), axis.title = element_text(size=12), 
#           legend.text = element_text(size=11), legend.key = element_rect(fill = 'white'), legend.position = 'bottom')+
#     ylab(paste0('Distance ',ylab_,'\n(method=canberra)'))+ggtitle(ggtitle_)+xlab('Epitopes')+
#     scale_x_discrete(expand = expansion(mult = c(.01,.01)))+
#     scale_y_break(c(y_l,y_u))+scale_y_continuous(labels = scales::comma)+
#     theme(axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank())
# }

# wat_dotPlot_self(AutoimmuneSelfCR)
# wat_dotPlot_viral(AutoimmuneViralCR)

unique(AutoimmuneSelfCR$Peptide) %>% length()
unique(AutoimmuneViralCR$Peptide) %>% length()


# Get the self Epitopes
dplyr::filter(HuVi_netMHC, Peptide=='VLFGLGFAI') %>%
  mutate(Gene=mapIds(EnsDb.Hsapiens.v86, Identity, 'SYMBOL', 'PROTEINID'))

na.omit(AutoimmuneSelfCR) %>%
  dplyr::filter(Level=='High' & Cell.type=='pancreatic endocrine cells') %>%
  ungroup() %>% slice_min(order_by = value, n=1) %>% group_by(Epitope, Peptide, Gene.name, value) %>% summarise() %>%
  formattable::formattable()
  
na.omit(AutoimmuneSelfCR) %>%
  dplyr::filter(Level=='High' & Cell.type=='pancreatic endocrine cells') %>%
  ungroup() %>% slice_max(order_by = value, n=2) %>% group_by(Epitope, Peptide, Gene.name, value) %>% summarise() %>%
  formattable::formattable()


# Get the viral Epitopes
tmp = dplyr::filter(HuVi_netMHC, Peptide %in% AutoimmuneViralCR$Peptide) %>%
  mutate(ViralID=gsub('(\\S+_\\d+)_(\\d)_.*','\\1.\\2',Identity)) %>%
  inner_join(ViralID_) %>% 
  mutate(protein=gsub('','',fastaHeader),
         ViralSpecies=gsub('.*\\[(.*)\\]','\\1',fastaHeader)) %>%
  group_by(Peptide, ViralID, ViralSpecies) %>% summarise()

AutoimmuneViralCR = left_join(AutoimmuneViralCR, tmp)

AutoimmuneViralCR %>% slice_min(order_by = value, n=2)

AutoimmuneViralCR %>% slice_max(order_by = value, n=2)



#####




###
# Similarity between self and non-self epitopes
source('bin/self_viral_similarity.R')


# NOT USED, see if there's interest: Calculate distance between points in space and charge separately----
# Source to load other factors: https://github.com/vadimnazarov/kidera-atchley, replaced by immunarch
pHLA_seu_dist_cor <- qs::qread('results/RData/pHLA_seu_dist_cor_simple_qsave.rds')
##pHLA_seu_dist_cor <- qs::qread('results/RData_Final/pHLA_seu_trimm_dist_cor_simple_qsave.rds')

cols = RColorBrewer::brewer.pal(5, 'Set1')
pHLA_seu_dist_cor %>%
  dplyr::filter(distance_seu!='minkowski' & distance_image!='minkowski') %>%
  mutate(seu_=gsub('_charge','',seu_),
         seu_=gsub('epitope','Epítopo',seu_),
         seu_=gsub('pHLA_filt','ec pHLA',seu_),
         seu_=factor(seu_, levels = c("pepBinary", "pepLevshtein", "pepAtchley", "pHLA", "ec pHLA", "Epítopo")),
         perspective=factor(perspective, levels = c('Top','upTop','downTop','leftTop','rightTop'))) %>%
  ggplot(aes(x=seu_, y=spearman_rho, fill=perspective))+geom_col(position = 'dodge')+
  facet_grid(distance_image~distance_seu)+
  #  ggtitle('Spearman correlation coefficient of the distances between pHLA images and epitopes', subtitle = 'distance between pHLA images')+
  ggtitle('Correlação entre imagem do complexo pHLA e epítopos\n', subtitle = 'Distância entre os epítopos')+
  xlab('Grupo')+ylab('Coeficiente de Spearman (rho)')+
  guides(fill=guide_legend(title = 'Perspectiva'))+scale_fill_manual(values = cols)+
  scale_y_continuous(sec.axis = sec_axis(trans = ~., name = "Distância entre as imagens", breaks = NULL))+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=12), axis.title = element_text(size=14),
        strip.text = element_text(size=12), #strip.placement = 'outside', 
        strip.background = element_rect(colour = 'black', fill='white'),
        axis.title.y.right = element_text(size=12), plot.subtitle = element_text(hjust = .5, size=12),
        axis.line = element_line(), panel.background = element_blank())
#####

# NOT USED, see if there's interest: Dimension reduction to find cluster of peptides----
# DEPRECATED: seu_l do not contain radius info, replaced by seu_l_radius
# DEPRECATED: seu_trimm_l -> seu_trimm_l_radius

# 5a) peptide binary matrix, replaced by Levenshtein distance
##seu_l$pepBinary@assays$RNA[1:5,1:5]
# 5a) peptide Levenshtein distance
seu_l$pepLevshtein@assays$RNA[1:5,1:5]
# 5b) peptides encoded using Atchley factors
seu_l$pepAtchley@assays$RNA[1:5,1:5] # Position_ATCHLEYfactor1:5
# 5c) spatial and charge data using the entire pHLA, exposed pHLA residues, or peptide
seu_l$pHLA_charge@assays$RNA[1:5,1:5]

# 5) Figure DimPlot
col_dim_ = RColorBrewer::brewer.pal(4, 'Dark2')

plot_grid(
  lapply(names(seu_l_radius)[c(2:3,5:6)], function(x) DimPlot(seu_l_radius[[x]], reduction = 'tsne')+ggtitle(x)+NoLegend()) %>%
    plot_grid(plotlist = ., nrow=1),
  
  plot_grid(
    DimPlot(seu_l_radius[[2]], group.by = 'origin', reduction = 'tsne', cols=col_dim_)+ggtitle('pepLevshtein')+NoLegend(),
    DimPlot(seu_l_radius[[3]], group.by = 'origin', reduction = 'tsne', cols=col_dim_)+ggtitle('pepAtchley')+NoLegend(),
    DimPlot(seu_l_radius[[5]], group.by = 'origin', reduction = 'tsne', cols=col_dim_)+ggtitle('ec pHLA')+NoLegend(),
    DimPlot(seu_l_radius[[6]], group.by = 'origin', reduction = 'tsne', cols=col_dim_)+ggtitle('epitope'), 
    nrow=1, rel_widths = c(1,1,1,1.5)),
  
  plot_grid(
    DimPlot(seu_l_radius[[2]], group.by = 'origin', reduction = 'umap', cols=col_dim_)+ggtitle('pepLevshtein')+NoLegend(),
    DimPlot(seu_l_radius[[3]], group.by = 'origin', reduction = 'umap', cols=col_dim_)+ggtitle('pepAtchley')+NoLegend(),
    DimPlot(seu_l_radius[[5]], group.by = 'origin', reduction = 'umap', cols=col_dim_)+ggtitle('ec pHLA')+NoLegend(),
    DimPlot(seu_l_radius[[6]], group.by = 'origin', reduction = 'umap', cols=col_dim_)+ggtitle('epitope'), 
    nrow=1, rel_widths = c(1,1,1,1.5)), nrow=3
)

# LogoPlot of the AA frequency of the clusters
change_res <- function(seu_, res_){
  Idents(seu_) = seu_[[paste0('RNA_snn_res.', res_)]]
  return(seu_)
}


plot_grid(
  # First logoplot (left)
  lapply(c(0:3), function(ident_){
    ggplot()+geom_logo(subset(change_res(seu_l_radius$pepLevshtein, 0.1), idents=ident_) %>% colnames())+
      theme_logo()+ggtitle(paste('Cluster', ident_))+theme(legend.position = 'none')
  }) %>% plot_grid(plotlist = ., ncol=1),
  
  # DimPlot
  plot_grid(
    DimPlot(change_res(seu_l_radius$pepLevshtein, 0.1), reduction = 'tsne', label=TRUE)+ggtitle('pepLevshtein')+NoLegend(),
    DimPlot(change_res(seu_l_radius$pHLA_charge_filt, 0.1), reduction = 'tsne', label=TRUE)+ggtitle('ec pHLA')+NoLegend(), nrow=2),
  
  # Second logoplot (right)
  lapply(c(0:3), function(ident_){
    ggplot()+geom_logo(subset(change_res(seu_l_radius$pHLA_charge_filt, 0.1), idents=ident_) %>% colnames())+
      theme_logo()+ggtitle(paste('Cluster', ident_))+theme(legend.position = 'none')
  }) %>% plot_grid(plotlist = ., ncol=1),
  ncol=3, rel_widths = c(.5,1,.5)
)

#####



# Select recurrent paired TRA:TRB in databases
# problema de contagem aqui, o 12 possui 12 CR mas o mesmo CDR3
recurrent_CDR3_paired <- dplyr::select(peptide_cdr3_db, -db) %>% unique() %>% na.omit() %>%
  dplyr::filter(CDR3_beta!='' & CDR3_alpha!='') %>%
  # Filter for cross-reactives
  group_by(CDR3_beta, CDR3_alpha) %>% dplyr::filter(n()>1) %>% 
  # number of epitopes recognized by CDR3 pair
  mutate(n_Epitopes=n())

# >> go to Figures_Manuscript.R
# recurrent_CDR3_paired %>%
#   group_by(n_Epitopes, CDR3_beta, CDR3_alpha) %>% summarise() %>%
#   ggplot(aes(x=n_Epitopes))+geom_bar()+
#   ggtitle('Recurrent TRA:TRB')+ylab('CDR3 a:b pair')+xlab('Epitopes recognized')+
#   scale_x_continuous(breaks = c(2:6,10,12))+scale_y_sqrt(breaks = c(1,10,200,400,600), expand = expansion(mult = c(0, .1)))+
#   theme_bw()+theme(axis.text = element_text(size=11), axis.title = element_text(size=12))

# >> go to Figures_Manuscript.R
# Alignment of the epitopes recognized by CR TCRs
# peps_ = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==12) %>% pull(Epitope)
# p1=ggplot()+geom_logo(peps_)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab('Position')+
#   theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'bottom')
# 
# peps_ = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==10) %>% pull(Epitope)
# p2=ggplot()+geom_logo(peps_)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab('Position')+
#   theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'none')
# 
# plot_grid(p1, p2, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(1.5,1))

# Select most recurrent TRA or TRB
# recurrent_CDR3_single = peptide_cdr3_db %>%
#   dplyr::select(c(Epitope, origin, description, db, CDR3_beta, CDR3_alpha)) %>% unique() %>%
#   group_by(CDR3_beta, CDR3_alpha) %>%
#   dplyr::filter(n()>1) %>%
#   reshape2::melt(id.vars=c('Epitope', 'origin', 'db','description')) %>% na.omit() %>%
#   group_by(value) %>% mutate(n=n())
# 
# recurrent_CDR3_single %>%
#   ggplot(aes(x=n))+geom_bar()+facet_wrap(variable~., scales = 'free_y')+
#   ggtitle('Most recurrent CDR3')+xlab('Number of epitopes recognized')+
#   ggtitle('TRB CDR3 recorrentes')+xlab('Epítopos reconhecidos')+
#   scale_x_continuous(breaks=c(2,5,10,28), labels = c(2,5,10,28))

dist_all_pHLA_filt = dist(t(pep_trimm_l$pHLA_charge_filt), method = 'euclidean') %>% as.matrix()
dist_all_epitope = dist(t(pep_trimm_l$epitope_charge), method = 'euclidean') %>% as.matrix()
dist_all_Images = lapply(pHLA_masks_dist$canberra, as.matrix)

dist_all_pHLA_filt_melted = melt_dist_selfViral(dist_all_pHLA_filt)
dist_all_epitope_melted = melt_dist_selfViral(dist_all_epitope)
dist_all_Images_melted = lapply(dist_all_Images, melt_dist_selfViral)


# CrossReac of epitopes from the databases
ind_12 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==12) %>% pull(Epitope)
dplyr::filter(recurrent_CDR3_paired, n_Epitopes==12) %>% mutate(pull_=paste(CDR3_beta, CDR3_alpha)) %>% pull(pull_) %>% unique()
ind_10 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==10) %>% pull(Epitope)
dplyr::filter(recurrent_CDR3_paired, n_Epitopes==10) %>% mutate(pull_=paste(CDR3_beta, CDR3_alpha)) %>% pull(pull_) %>% unique()

# median difference
median_CR <- function(dist_melted_, ind_n_Tpitopes){
  dplyr::filter(dist_melted_, Var1 %in% ind_n_Tpitopes) %>%
    mutate(cat = ifelse(Var2 %in% ind_n_Tpitopes, 'pair', 'other')) %>% group_by(cat) %>% summarise(median(value))
}

median_CR(dist_all_epitope_melted, ind_12)
median_CR(dist_all_pHLA_filt_melted, ind_12)
median_CR(dist_all_Images_melted$Top, ind_12)

median_CR(dist_all_epitope_melted, ind_10)
median_CR(dist_all_pHLA_filt_melted, ind_10)
median_CR(dist_all_Images_melted$Top, ind_10)

ind_12 %in% dist_all_epitope_melted$Var1
ind_10 %in% dist_all_epitope_melted$Var1

# >> go to Figures_Manuscript.R
# Plot
plot_CR <- function(dist_melted_, ind_n_Tpitopes, dist_method_='euclidean', ggtitle_=NULL){
 # y_boxplot = quantile(dist_melted_$value)[3]
  
  dplyr::filter(dist_melted_, Var1 %in% ind_n_Tpitopes) %>%
    mutate(cat = ifelse(Var2 %in% ind_n_Tpitopes, 'CR', 'Other'), cat=factor(cat, c('Other','CR'))) %>%
    arrange(cat) %>%
    ggplot(aes(x=reorder(ID, value), y=value, color=cat))+geom_point(size=.9, alpha=.7)+
#    geom_boxplot(aes(x=value, y=y_boxplot), outlier.shape = NA)+ # , color='black'
    xlab('Epitopes')+ylab(paste0('Distance between epitopes\n(method=', dist_method_, ')'))+
    ggtitle(ggtitle_)+
    theme(axis.text = element_blank(), axis.line = element_line(), axis.title = element_text(size=11),
          panel.background = element_blank(),
          legend.title = element_blank(), legend.key = element_rect(fill='white'), plot.title = element_text(hjust=.5))+
    scale_color_manual(values = c('CR'='black','Other'='gray'))+
    scale_x_discrete(expand = c(.01,.01))
}

# >> go to Figures_Manuscript.R
# p1=plot_CR(dist_all_Images_melted$Top, ind_12, 'canberra', 'pHLA Top images')+NoLegend()
# p2=plot_CR(dist_all_pHLA_filt_melted, ind_12, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+NoLegend()
# p3=plot_CR(dist_all_epitope_melted, ind_12, ggtitle_ = 'Epitope (tr_c+XYZ)')
# 
# p4=plot_CR(dist_all_Images_melted$Top, ind_10, 'canberra', 'pHLA Top images')+NoLegend()
# p5=plot_CR(dist_all_pHLA_filt_melted, ind_10, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+NoLegend()
# p6=plot_CR(dist_all_epitope_melted, ind_10, ggtitle_ = 'Epitope (tr_c+XYZ)')
# 
# plot_grid(plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,1,1.25)),
#           plot_grid(p4, p5, p6, nrow = 1, rel_widths = c(1,1,1.25)), nrow=2)

#####

