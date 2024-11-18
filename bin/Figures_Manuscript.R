# Script to generate the Figures for Pretti 2024
quit()
pacman::p_load(ggplot2, tibble, RColorBrewer, Seurat, cowplot, corrplot, Hmisc, immunarch, ggseqlogo, 
               qs, pheatmap, hpar, ggbreak, ggplotify, aplot)
set.seed(123)
options(dplyr.summarise.inform = FALSE)

load("results/workspace202402_light.RData")
#load("/results/workspace202302.RData")
setwd('results/Figures_Manuscript_202408/')
#

# themes
theme_ggplot_ = list(theme_bw(),
                     theme(axis.text = element_text(size=12, color = 'black'), axis.title = element_text(size=12, color = 'black'), 
                           plot.title = element_text(hjust = .5)))

#

# OK: Figure 1----
p1=recurrent_CDR3_paired %>%
  group_by(n_Epitopes, CDR3_beta, CDR3_alpha) %>% summarise() %>%
  ggplot(aes(x=n_Epitopes))+geom_bar()+
  ggtitle('Recurrent TRA:TRB')+ylab('CDR3 a:b pair')+xlab('Epitopes recognized')+
  scale_x_continuous(breaks = c(2:6,10,12))+scale_y_sqrt(breaks = c(1,10,200,400,600), expand = expansion(mult = c(0, .1)))+
  theme_bw()+
  theme(axis.text = element_text(size=13, color = 'black'), axis.title = element_text(size=14, color = 'black'), 
        plot.title = element_text(hjust = .5))

peps_12 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==12) %>% pull(Epitope)
peps_10 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==10) %>% pull(Epitope)

t1=ggdraw()+
  draw_label(paste0(peps_12, collapse = '\n'), fontfamily = 'mono', x = .3, size=15)+
  draw_label(paste0(c('','',peps_10), collapse = '\n'), fontfamily = 'mono', x = .7, size=15)+
  geom_rect(aes(xmin=.15, xmax=.45, ymin=.05, ymax=.95), fill='transparent', color='black')+
  geom_rect(aes(xmin=.55, xmax=.85, ymin=.05, ymax=.8), fill='transparent', color='black')+
  geom_segment(aes(x=.45, xend=.9, y=.9, yend=.9), arrow = arrow())+
  geom_segment(aes(x=.85, xend=1, y=.4, yend=.3), arrow = arrow())

#t1

p2=ggplot()+geom_logo(peps_12)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab('Position')+
  theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'bottom', legend.title = element_blank())

p3=ggplot()+geom_logo(peps_10)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab('Position')+
  theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'none')

# plot_grid(p1, t1,
#           plot_grid(p2, p3, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(1.5,1)),
#           nrow = 1, rel_widths = c(.7,1,1))

#
tema_fig1_=list(theme(axis.title = element_text(size=14, color = 'black'), legend.text = element_text(size=12)),
                guides(color=guide_legend(override.aes = list(size=2))))
p4=plot_grid(
  plot_CR(dist_all_Images_melted$Top, ind_12, 'canberra', 'pHLA Top images')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_pHLA_filt_melted, ind_12, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_epitope_melted, ind_12, ggtitle_ = 'Epitope (tr_c+XYZ)')+tema_fig1_,
  nrow = 1, rel_widths = c(1,.1,1,.1,1.25)
)

p5=plot_grid(
  plot_CR(dist_all_Images_melted$Top, ind_10, 'canberra', 'pHLA Top images')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_pHLA_filt_melted, ind_10, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_epitope_melted, ind_10, ggtitle_ = 'Epitope (tr_c+XYZ)')+tema_fig1_, 
  nrow = 1, rel_widths = c(1,.1,1,.1,1.25)
)

# plot_grid(plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,1,1.25)),
#           plot_grid(p4, p5, p6, nrow = 1, rel_widths = c(1,1,1.25)), nrow=2)

# final
png('_Figure1.png', width = 1850, height = 1400, res = 150)
plot_grid(plot_grid(p1, t1,
          plot_grid(p2, p3, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(1.5,1)),
          nrow = 1, rel_widths = c(1,1,1.2)),
          NULL,
          p4, 
          NULL,
          p5, nrow = 5, rel_heights = c(1,.1,1,.1,1))

dev.off()
#####

# OK: Figure 2----
cols = colorRampPalette(c('blue','white','red'))
myBreaks <- c(seq(min(pre_plot_bulek_$Top, na.rm=TRUE), 0, length.out=ceiling(100/2) + 1), 
              seq(max(pre_plot_bulek_$Top, na.rm=TRUE)/100, max(pre_plot_bulek_$Top, na.rm=TRUE), 
                  length.out=floor(100/2)))

# single heatmap
png('_Figure2.png', width = 800, height = 800, res=150)
pheatmap::pheatmap(t(pre_plot_bulek_$Top), angle_col = 0, cluster_cols = FALSE, color = cols(100),
                   breaks = myBreaks, cluster_rows = FALSE, fontsize = 14)
dev.off()
#####

# OK: Figure 3----
plot_umap_ <- function(umap_, ggtitle_=NULL){
  as.data.frame(umap_@cell.embeddings) %>%
    rownames_to_column('Epitopes') %>%
    inner_join(meta_bulek %>% rownames_to_column('Epitopes')) %>%
    mutate(`TNF release\n`=as.numeric(TNFsecretion_log10)) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, fill=`TNF release\n`))+geom_point(color='black', shape=21)+
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 1,
                         labels=c(0,3,10,30,100,320), guide = 'colorbar')+
    ggtitle(ggtitle_)+
    #    guides(color=guide_legend(title = 'Liberação de TNF'))+
    theme(axis.line = element_line(), panel.background = element_blank(), #legend.position = 'bottom',
          axis.text = element_text(size=11), axis.title = element_text(size=12))
}

# A) plot_umap_(AddMetaData(seu_pHLA_masks_bulek_l$Top, meta_bulek)@reductions$umap)
pB <- plot_umap_(umap_pHLA_masks_bulek_l[[1]])

# B-C) Correlations, dispersion plots for charge and image
tema_ = list(theme_bw(), theme(axis.text = element_text(size=12), axis.title = element_text(size=12)))

# image
pC <- plot_cor_bulek(dist_bulek_(pre_umap_pHLA_masks_bulek$Top, dist_method_ = 'canberra'),
               dist_method_ = 'canberra')+ggtitle('Image: Top')+tema_

# charge
pD <- plot_cor_bulek(dist_bulek_(seu_bulek_trimm_l_radius$epitope_charge@assays$RNA@counts, dist_method_ = 'canberra'),
                     dist_method_ = 'canberra')+ggtitle('Epitope: tr_c+XYZ+r')+tema_


png('_Figure3.png', width = 1600, height = 450, res=150)
plot_grid(pB, NULL, pC, NULL, pD, nrow = 1, rel_widths = c(1.3,.1,1,.1,1))
dev.off()
#####

# OK: Figure 4----
tema_fig4_ = list(theme(axis.text = element_text(color='black'), plot.title = element_text(hjust = .5)))
plot_umap_lee_ <- function(umap_, ggtitle_=NULL){
  as.data.frame(umap_@cell.embeddings) %>%
    rownames_to_column('Epitopes') %>%
    inner_join(meta_lee[!is.na(meta_lee$IFNy_secr_g10),] %>% rownames_to_column('Epitopes')) %>%
    mutate(`IFNy production\n(g10)\n`=IFNy_secr_g10) %>%
    ggplot(aes(x=UMAP_1, y=UMAP_2, color=`IFNy production\n(g10)\n`))+geom_point()+
    scale_color_manual(values = rev(RColorBrewer::brewer.pal(11, 'RdBu'))[c(5,7,8,10,11)])+
    ggtitle(ggtitle_)+
    theme(axis.line = element_line(), panel.background = element_blank(), #legend.position = 'bottom',
          legend.key = element_rect(fill = 'white', color='black'),
          axis.text = element_text(size=11), axis.title = element_text(size=12))+tema_fig4_
}

pB <- plot_umap_lee_(AddMetaData(seu_pHLA_masks_lee_l$Top, meta_lee)@reductions$umap, ggtitle_ = 'Image: Top')

# C-D) Correlations, dispersion plots for charge and image

# image
pC <- plot_cor_lee(dist_lee_image(pre_umap_pHLA_masks_lee$Top, dist_method_ = 'canberra'),
                   dist_method_ = 'canberra')+ggtitle('Image: Top (g10)')+tema_fig4_

# charge
pD <- plot_cor_lee(dist_lee(seu_lee_trimm_l$pHLA_charge_filt, dist_method_ = 'canberra', IFNy_ = 'IFNy_secr_t5'),
                   dist_method_ = 'canberra')+ggtitle('Epitope: tr_c+XYZ (t5)')+tema_fig4_


# D)
ymax_=ifelse(max(plot_summ$rho)<0,0,max(plot_summ$rho))
ymin_=min(plot_summ$rho)*1.1

pE = plot_summ %>%
  ggplot(aes(x=x_axis, y=rho, fill=fill))+
  geom_col(position = position_dodge2(preserve = 'single'), color='black')+
  geom_hline(yintercept = 0)+
  geom_text(aes(label='*'), dplyr::filter(plot_summ, !grepl('Image', x_axis) & p.value<.05), position = position_dodge2(width = .9, preserve = 'single'),
            vjust=1.205, size=6)+
  ylab("Correlation coefficient\n(spearman's rho)")+xlab('Groups')+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  facet_grid(~group, scales = 'free_x', space = 'free_x', switch = 'x')+
  theme(axis.line = element_line(), panel.background = element_blank(),
        strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size = 11),
        legend.text = element_text(size=12),
        axis.text = element_text(size=11), axis.title = element_text(size=12))+tema_fig4_+
  scale_y_continuous(limits = c(ymin_, ymax_))+scale_fill_manual(values = c('white','gray'))+
  guides(fill=guide_legend(title = 'TCR clone'))

# pE
png('_Figure4.png', width = 1600, height = 900, res=150)
plot_grid(plot_grid(pB, pC, NULL, pD, nrow=1, rel_widths = c(1,.75,.05,.75)), 
          NULL, 
          plot_grid(NULL, pE, rel_widths = c(.05,1)), nrow=3, rel_heights = c(1,.1,1))
dev.off()
#####

# OK: Figure 5----
plot_dist_selfViral_melted = function(dist_selfViral_melted_, subtitle_=NULL){
  dist_selfViral_melted_ %>% 
    ggplot(aes(x=selfViral, y=value))+geom_violin()+geom_boxplot(outlier.size = 0)+
    ggpubr::stat_compare_means(hjust=0)+
    #    ggtitle('Similarity between human and viral', subtitle_)+
    ggtitle(subtitle_)+
    ylab(paste0('Distance between ', subtitle_, '\n(method=canberra)'))+xlab('')+theme_bw(base_size = 16)+
    scale_y_continuous(labels = scales::label_comma())+
    theme(axis.text = element_text(color = 'black'), axis.title = element_text(color = 'black'),
          plot.title = element_text(hjust = .5, size=12))
}

p1=plot_dist_selfViral_melted(dist_selfViral_Images_melted$Top, '')


wat_dotPlot <- function(dist_HumVir_melted_, ggtitle_='Similarity between human and viral epitopes', ylab_=NULL, scale_break_=TRUE){
  require('ggbreak')  
  
  if(scale_break_){
    tmp_ = bind_rows(dist_HumVir_melted_ %>% ungroup() %>%
                       # get the top 300, distant epitopes between them
                       arrange(value) %>% top_n(300, wt = -value),
                     dist_HumVir_melted_ %>% ungroup() %>%
                       # get the 300, closer epitopes between them
                       arrange(value) %>% top_n(300, wt = value))
  }else{
    tmp_ = dist_HumVir_melted_
  }
  
  
  y_l = sort(tmp_$value) %>% head(300) %>% tail(1)*1.01
  y_u = sort(tmp_$value) %>% tail(300) %>% head(1)*0.99
  print(c(y_l, y_u))
  
  p=tmp_ %>%
    ggplot(aes(x=reorder(ID, value), y=value, color=selfViral))+geom_point(size=.7)+
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.background = element_blank(),
          axis.text.y = element_text(size=12, color = 'black'),
          axis.line = element_line(), axis.title = element_text(size=14, color = 'black'), 
          axis.title.x = element_text(hjust = .6),
          legend.title = element_blank(), legend.key = element_blank(),
          legend.text = element_text(size=12), plot.title = element_text(hjust=.6, size=14))+
    ylab(paste0('\nDistance between ',ylab_,'\n(method=canberra)'))+xlab('Ranked epitope pairs')+
    #ggtitle(ggtitle_)+
    scale_y_continuous(labels = scales::label_comma(), sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL))+
    scale_color_manual(values = c('black','red','gray'))+
    guides(color=guide_legend(override.aes = list(size=1.5)))
  
  if(scale_break_){
    p+scale_y_break(c(y_l,y_u))
  }else{
    return(p)
  }
}

# p2=wat_dotPlot(dist_selfViral_Images_melted$Top, 'pHLA Top perspective', 'images', scale_break_ = FALSE)+NoLegend()

p2=wat_dotPlot(dist_selfViral_Images_melted$Top, 'Top 300 and bottom 300 epitope pairs', 'images')+
  theme(axis.title.x = element_text(hjust = .35), axis.title = element_text(size=16), plot.title = element_text(hjust = .25, size=14),
        axis.text = element_text(color='black', size=16))

aplot::plot_list(p1, p3, nrow = 1, widths = c(1,1.35))

# p3
plot_CR_dotplot = function(CancerViralCR_, dist_selfViral_, ggtitle_=NULL, annot_=NULL){

  limits_ = quantile(CancerViralCR_$value_min)[c(1,5)]*c(.95,1.05)
  width_ = max(limits_)/50
  
  CancerViralCR_ %>% 
    reshape2::dcast(Epitope+Pathology~selfViral) %>%
    mutate(Pathology=factor(Pathology, levels=c('Neoantigen','AML','Lymphoma','Melanoma, Neoantigen','TAA'))) %>% arrange(Pathology) %>%
    # annnotate epitopes
    mutate(label_=ifelse(Epitope %in% annot_, Epitope, '')) %>%
    ggplot(aes(x=`Human/Human`, y=`Human/Viral`, label=label_))+geom_abline()+geom_point()+
    ggrepel::geom_text_repel(max.overlaps = Inf, min.segment.length = 0, size=3.5, box.padding = 1)+
    scale_y_continuous(limits = limits_, labels = scales::comma)+
    scale_x_continuous(limits = limits_, labels = scales::comma)+
    ggtitle(ggtitle_)+
    ylab('Distance to the closest viral epitope')+xlab('Distance to the closest self epitope')+
    theme_bw(base_size = 16)+theme(axis.title = element_text(color = 'black'), axis.text = element_text(color = 'black'), 
                                   plot.title = element_text(hjust = .5),
                     legend.text = element_text(size=10), legend.position = 'bottom', plot.margin = margin(5,15,5,5))+
    guides(color=guide_legend(title = element_blank()))
}

process_CancerViralCR <- function(CancerViralCR_){
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

# Annotate neoantigens and WT peptides from Zhang 2018
neoantigens_zhang <- readxl::read_xlsx('results/Sup_fromZhang2018_Neoantigens_NIHMS1508189-supplement-8.xlsx', skip = 1)

CancerViralCR_Top_noWT_ <- CancerViralCR_Top %>%
  dplyr::filter(!Epitope %in% neoantigens_zhang$`Wildtype peptide` &
                  !Var2 %in% neoantigens_zhang$`Wildtype peptide`) 

CancerViralCR_Top_noWT <- CancerViralCR_Top_noWT_ %>%
  process_CancerViralCR()

# CancerViralCR_Top_filtered_noWT <- CancerViralCR_Top_noWT %>%
#   dplyr::filter(Pathology!='Neoantigen' | Epitope %in% neoantigens_zhang$`Mutant peptide`)
# 
# CancerViralCR_Top_filtered_ <- CancerViralCR_Top_ %>%
#   dplyr::filter(Pathology!='Neoantigen' | Epitope %in% neoantigens_zhang$`Mutant peptide`)


# debug
p=plot_CR_dotplot(CancerViralCR_Top_noWT, dist_selfViral_ = dist_selfViral_Images_melted$Top, 
                ggtitle_ = 'pHLA Top perspective', annot_ = CancerViralCR_Top_$Epitope)
ggplotly(p)

annot_ = c('FLCMKALLL', 'SLAETFLET', 'FLIYLDVSV')
p3=plot_CR_dotplot(CancerViralCR_Top_noWT, dist_selfViral_ = dist_selfViral_Images_melted$Top, 
                ggtitle_ = '', annot_ = annot_)

# Figure
png('Figure_5.png', width = 2500, height = 800, res = 150)
#svg('Figure_5.svg', width = 17, height = 5)
plot_grid(aplot::plot_list(p1, p2, p3, nrow = 1, widths = c(1.05,1.3,1.05)))
dev.off()

# Rank for pHLA images
lookup_closest_epitope('FLCMKALLL') # center

lookup_closest_epitope('SLAETFLET') # bottom

lookup_closest_epitope('FLIYLDVSV') # top

#dplyr::filter(dist_selfViral_Images_melted$Top, grepl('YLLNYDLSV', ID))

lookup_closest_epitope <- function(pep_TAA_){
  df_ <- dist_selfViral_Images_melted$Top %>%
    dplyr::filter(!Var1 %in% neoantigens_zhang$`Wildtype peptide` &
                    !Var2 %in% neoantigens_zhang$`Wildtype peptide`)
  
  pep_Viral_ <- dplyr::filter(df_, grepl(pep_TAA_, ID) & selfViral=='Human/Viral') %>%
    arrange(value) %>% head(2) %>% pull(ID) %>% as.character() %>% gsub(pep_TAA_, '', .) %>% gsub('_','',.)
  
  pep_self_ <- dplyr::filter(df_, grepl(pep_TAA_, ID) & selfViral=='Human/Human') %>% 
    arrange(value) %>% head(2) %>% pull(ID) %>% as.character() %>% gsub(pep_TAA_, '', .) %>% gsub('_','',.)

  # calculate unified rank for viral and self similarity
  # no sense to combine self and viral
  t1=dplyr::filter(dist_selfViral_Images_melted$Top, Var1 %in% pep_Viral_) %>%
    mutate(rank_=rank(value, ties.method = 'min'),
           rank_=paste0(rank_, '/', length(rank_))) %>%
    dplyr::filter(Var2 == pep_TAA_)

  t2=dplyr::filter(dist_selfViral_Images_melted$Top, Var1 %in% pep_self_) %>%
    mutate(rank_=rank(value, ties.method = 'min'),
           rank_=paste0(rank_, '/', length(rank_))) %>%
    dplyr::filter(Var2 == pep_TAA_)

  print(bind_rows(t1, t2))
}

# Stats: greater or lower than slope=1
tmp_ = CancerViralCR_Top_noWT %>%
  reshape2::dcast(Epitope~selfViral, value.var = 'value_min')

table(tmp_$`Human/Human`> tmp_$`Human/Viral`)
#####




# Ressampling (Bootstrap) of dispersion plot in Figure 6 to observe if the distribution remains

CancerViralCR_Top_noWT_bootstrap <- lapply(1:1000, function(time_){
  lapply(c(100,200,300,400,500), function(size_){
    ind_ = sample(1:nrow(CancerViralCR_Top_noWT_), size = size_, replace = FALSE)
    CancerViralCR_Top_noWT_[ind_,] %>% process_CancerViralCR() %>% 
      group_by(selfViral) %>% summarise(value_min=median(value_min), boot=time_, size = size_)
  }) %>% bind_rows()
}) %>% bind_rows()

library(rstatix)
cohen_ <- CancerViralCR_Top_noWT_bootstrap %>% ungroup() %>%
  group_by(size) %>%
  cohens_d(value_min~selfViral)

cohen_$effsize = signif(cohen_$effsize, 3)

CancerViralCR_Top_noWT_bootstrap %>%
  ggplot(aes(x=selfViral, y=value_min))+
  geom_violin()+stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1))+
  facet_grid(~size)+ylab('Distance to the closest epitope')+xlab('')+
#  ggpubr::stat_compare_means()
  ggpubr::stat_pvalue_manual(cohen_, label = 'effsize', y.position = 2400)+
  scale_y_continuous(labels = scales::label_comma())+
  theme_bw(base_size = 14)+theme(strip.background = element_rect(fill = 'white'), axis.text = element_text(color='black'))

# plot_CR_dotplot(CancerViralCR_Top_noWT, dist_selfViral_ = dist_selfViral_Images_melted$Top, 
#                 ggtitle_ = 'pHLA Top perspective', annot_ = NULL)
#####


# Figure 6 heatmap of rgb images----
library(png)
FLIYLDVSV <- readPNG('data/modeled/FLIYLDVSV/FLIYLDVSV_scale10_Top.png')
YLLNYDLSV <- readPNG('data/modeled/YLLNYDLSV/YLLNYDLSV_scale10_Top.png')
YMIGTDFYV <- readPNG('data/modeled/YMIGTDFYV/YMIGTDFYV_scale10_Top.png')

pdf('Figure5D_cont.pdf', height = 5)
diff_electrostatic_pot_(FLIYLDVSV, YLLNYDLSV)
diff_electrostatic_pot_(FLIYLDVSV, YMIGTDFYV)
dev.off()

SLAETFLET <- readPNG('data/modeled/SLAETFLET/SLAETFLET_scale10_Top.png')
MMAELPFEV <- readPNG('data/modeled/MMAELPFEV/MMAELPFEV_scale10_Top.png')
MLADVVFEI <- readPNG('data/modeled/MLADVVFEI/MLADVVFEI_scale10_Top.png')


pdf('Figure5E_cont.pdf', height = 5)
diff_electrostatic_pot_(SLAETFLET, MMAELPFEV)
diff_electrostatic_pot_(SLAETFLET, MLADVVFEI)
dev.off()


FLCMKALL <- readPNG('data/modeled/FLCMKALLL/FLCMKALLL_scale10_Top.png')
FLAHVLNPV <- readPNG('data/modeled/FLAHVLNPV/FLAHVLNPV_scale10_Top.png')
FLLEKPFSV <- readPNG('data/modeled/FLLEKPFSV/FLLEKPFSV_scale10_Top.png')

pdf('Figure5F_cont.pdf', height = 5)
diff_electrostatic_pot_(FLCMKALL, FLAHVLNPV)
diff_electrostatic_pot_(FLCMKALL, FLLEKPFSV)
dev.off()


diff_electrostatic_pot_ <- function(x, y){
  # red-blue+transparency
  x = x[,,1]-x[,,2]+x[,,3]
  y = y[,,1]-y[,,2]+y[,,3]
  sub_ = x-y
#  print(quantile(sub_))
  
  cols = colorRampPalette(c('blue','white','red'))
  cols = cols(100)
  paletteLength = length(cols)
  myBreaks <- c(seq(min(sub_, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(sub_, na.rm=TRUE)/paletteLength, max(sub_, na.rm=TRUE), length.out=floor(paletteLength/2)))
  
  pheatmap(sub_, cluster_rows = FALSE, cluster_cols = FALSE, color = cols, breaks = myBreaks, legend_breaks = c(-.5,0,.5), legend_labels = c(-.5,0,.5))
}

#####


# OK: Figure S2----
tema=list(theme(axis.text = element_text(size=11, color='black'), axis.title = element_text(size=12, color='black'), axis.line = element_line(color='black'),
                panel.background = element_blank()),scale_y_continuous(expand = c(0,NA)))

# prePlot_charge_HuVi <- pep_l$epitope_charge[grepl('Charge', rownames(pep_l$epitope_charge)),rownames(meta_data_HuVi)]
# prePlot_charge_bulek <- pep_l$epitope_charge[grepl('Charge', rownames(pep_l$epitope_charge)),rownames(meta_bulek)]

p1_ <- prePlot_charge_bulek %>%
  `rownames<-`(gsub('Epitope_(\\d)_Charge','P\\1',rownames(.))) %>% t() %>%
  pheatmap(cluster_cols = F, color = colorRampPalette(brewer.pal(n = 6, name = "Greens"))(100), 
           show_rownames = FALSE, border_color = NA, angle_col = 0, silent = TRUE)

p2_ <- prePlot_charge_HuVi %>%
  `rownames<-`(gsub('Epitope_(\\d)_Charge','P\\1',rownames(.))) %>% t() %>%
  pheatmap(cluster_cols = F, color = colorRampPalette(brewer.pal(n = 6, name = "Blues"))(100), 
           show_rownames = FALSE, border_color = NA, angle_col = 0, silent = TRUE)


prePlot_spatial <- pep_l$epitope_charge[!grepl('Charge', rownames(pep_l$epitope_charge)),] %>%
  rownames_to_column('ID') %>%
  reshape2::melt() %>% 
  mutate(coord=gsub('Epitope_\\d_','',ID),
         Pos=gsub('Epitope_(\\d)_\\S+', 'P\\1', ID)) 


prePlot_spatial_bulek = dplyr::filter(prePlot_spatial, variable %in% rownames(meta_bulek)) %>%
  # absolute positions are not meaningful when comparing different peptides
  group_by(coord, Pos) %>% summarise(sd=sd(value), q1=quantile(value)[2], q3=quantile(value)[4], inter_q=q3-q1) %>%
  group_by(Pos) %>% mutate(mean_sd=mean(sd), mean_inter_q=mean(inter_q))

prePlot_spatial_HuVi = dplyr::filter(prePlot_spatial, variable %in% rownames(meta_data_HuVi)) %>%
  group_by(coord, Pos) %>% summarise(sd=sd(value), q1=quantile(value)[2], q3=quantile(value)[4], inter_q=q3-q1) %>%
  group_by(Pos) %>% mutate(mean_sd=mean(sd), mean_inter_q=mean(inter_q)) 

p3_ <- prePlot_spatial_bulek %>%
  ggplot(aes(x=Pos, fill=coord))+
  geom_col(aes(y=sd+inter_q), position = 'dodge')+
  geom_col(aes(y=inter_q), position = 'dodge', color='black')+
  geom_line(aes(y=mean_inter_q+mean_sd, group=coord))+
  ylab('Dispersion index\n(IQR+SD)')+
  tema+theme(axis.title.x = element_blank())+guides(fill=guide_legend(title = 'coordinate'))


p4_ <- prePlot_spatial_HuVi %>%
  ggplot(aes(x=Pos, fill=coord))+
  geom_col(aes(y=sd+inter_q), position = 'dodge')+
  geom_col(aes(y=inter_q), position = 'dodge', color='black')+
  geom_line(aes(y=mean_inter_q+mean_sd, group=coord))+
  ylab('Dispersion index\n(IQR+SD)')+
  tema+theme(axis.title.x = element_blank())+guides(fill=guide_legend(title = 'coordinate'))#+guides(fill='none')

# plot_grid(p3_, p4_)


png('_FigureS2.png', width = 1400, height = 800, res=150)
svg('Figure_S2.svg', width = 7.5, height = 5)
plot_grid(
  NULL, NULL, NULL,
  plot_grid(p1_$gtable, labels = 'A', label_y = 1.1), NULL, plot_grid(p3_, labels = 'B', label_y = 1.1),
  NULL, NULL, NULL,
  plot_grid(p2_$gtable, labels = 'C', label_y = 1.1), NULL, plot_grid(p4_, labels = 'D', label_y = 1.1), 
  nrow=4, rel_widths = c(1,.05,1.2), rel_heights = c(.1,1,.1,1), align = 'hv', axis = 'lr'
)
dev.off()

#####

# OK: Figure S3----
cols = colorRampPalette(c('white','grey'))

# Save images of the pHLA marked for high SD
library(ggplotify)
border_ = geom_rect(mapping = aes(xmin=0,xmax=1,ymin=0,ymax=1), color='black', fill='transparent')

figS3_l = lapply(names(pHLA_masks_sd_mat_l), function(x){
  pHLA_masks_sd_mat_l_ = pHLA_masks_sd_mat_l[[x]]
#  pHLA_masks_sd_mat_l_[50:430,100:540]
  p=pheatmap::pheatmap(pHLA_masks_sd_mat_l_[70:400,100:520], cluster_rows = FALSE, cluster_cols = FALSE, color = cols(5), 
                       silent = TRUE, main = x)
  #  png(paste0('results/Figures/pHLA_masks_sd_', x, '.png'), 640*2, 480*2, res=150)
  as.ggplot(p[[4]])+border_
  #  print(p)
#  dev.off()
})

png('_FigureS3.png', width = 1600, height = 1200, res = 150)
svg('FigureS3.svg', width = 9.5, height = 8)
plot_grid(NULL, figS3_l[[4]], NULL,
          figS3_l[[3]], figS3_l[[1]], figS3_l[[2]],
          NULL, figS3_l[[5]], NULL,
          nrow = 3)

dev.off()
#####

# OK: Figure S4----
tema_ = list(theme(legend.position = 'none', plot.title = element_text(hjust = .5), 
                   axis.text = element_text(size=12, color='black')))

legend_ = plot_umap_(umap_pHLA_masks_bulek_l[[1]]) %>% get_legend() %>% as.ggplot()
list_ = lapply(names(umap_pHLA_masks_bulek_l)[c(1:6)], function(x) plot_umap_(umap_pHLA_masks_bulek_l[[x]])+ggtitle(x)+tema_)

pA = plot_grid(plotlist = c(list_[1:3], NA,
                            rep(NA, 3), list(legend_),
                            list_[4:6], NA), ncol = 4, nrow=3, rel_heights = c(1,.05,1), rel_widths = c(1,1,1,.5))
pA

# pB
plot_bulek_summ_ = plot_bulek_summ %>% 
  mutate(group=gsub('Control','ctrl',group),
         # for displaying purposes
         group=gsub('Image\\+','Image\n+ ',as.character(group)),         
         group=factor(group, levels=c("ctrl","pHLA","ec pHLA","epitope","Image","Image\n+ ctrl",
                                      "Image\n+ pHLA","Image\n+ ec pHLA","Image\n+ epitope")))

pB <- plot_bulek_summ_ %>% 
  ggplot(aes(x=x_axis, y=rho))+
  geom_col(position = position_dodge2(preserve = 'single'), fill='white', color='black')+
  geom_hline(yintercept = 0)+
  geom_text(aes(label='*'), dplyr::filter(plot_bulek_summ_, !grepl('Image', x_axis) & p.value<.05), vjust=1.205, size=6)+
  ylab("Correlation coefficient\n(spearman's rho)")+xlab('Groups')+
  facet_grid(~group, scales = 'free_x', space = 'free_x', switch = 'x')+
  theme(axis.line = element_line(), panel.background = element_blank(),
        strip.placement = 'outside', strip.background = element_blank(), strip.text = element_text(size = 12),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.text = element_text(size=12, color='black'), axis.title = element_text(size=14, color='black'))+
  scale_y_continuous(limits = c(min(plot_bulek_summ_$rho)*1.1, max(plot_bulek_summ_$rho)))

#pB

png('_FigureS4.png', width = 2000, height = 1600, res = 150)
svg('Figure_S4.svg', width = 10.5, height = 9)
plot_grid(plot_grid(NULL, pA, NULL, rel_widths = c(.1,1,.1), nrow = 1, labels = 'A', label_size = 16), 
          NULL, plot_grid(pB, labels = 'B', label_y = 1.1, label_size = 16), nrow=3, rel_heights = c(1.2,.05,1))
dev.off()
#####

# OK: Figure S5 ou S6?----
tema_ = list(theme(axis.title = element_text(size=11, color='black'), axis.text = element_text(size=10, color='black'), 
                   plot.title = element_text(size=12, color='black'), strip.background = element_rect(fill='white'),
                   legend.text = element_text(size=9), legend.title = element_text(size=10), line = element_line(color = 'black')))

p1=HuVi_netMHC_summ %>%
  mutate(n=ifelse(n>10,10,n)) %>%
  group_by(Pep_len, source, n) %>% summarise(overRep=n()) %>%
  ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
  facet_grid(~Pep_len)+scale_y_log10(labels = scales::label_comma())+scale_x_continuous(breaks = 1:10)+
  ggtitle('Overrepresented epitopes')+
  ylab('Epitopes')+xlab('Overrepresentation')+
  theme_bw()+tema_+guides(color=guide_legend(title = 'origin'))

plot_ = HuVi_netMHC_shared %>%
  mutate(overRep=ifelse(n>10,10,n)) %>%
  group_by(overRep, Pep_len, carac) %>% summarise(n=n()) %>%
  group_by(Pep_len, carac) %>% mutate(prop=prop.table(n)) 

p2=plot_ %>%
  ggplot(aes(x=overRep, y=prop, color=carac))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = 1:10)+facet_grid(~Pep_len)+
  ggtitle('Shared epitopes')+
  ylab('Proportion of Epitopes')+xlab('Overrepresentation')+
  scale_color_manual(values = c('black','gray'))+theme_bw()+tema_+
  guides(color=guide_legend(title = 'epitope'))

HuVi_netMHC_summ_tr = HuVi_netMHC %>%
  dplyr::filter(nchar(Peptide) %in% c(9,10)) %>%
  mutate(Peptide_trimm=gsub('(\\S)\\S(.*)\\S','\\1\\2',Peptide)) %>%
  group_by(Peptide_trimm, source) %>% summarise(n=n())

p3=HuVi_netMHC_summ_tr %>%
  mutate(n=ifelse(n>10,10,n),
         Pep_len=nchar(Peptide_trimm)) %>%
  group_by(source, n, Pep_len) %>% summarise(overRep=n()) %>%
  ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
  scale_y_log10(labels = scales::label_comma())+scale_x_continuous(breaks = 1:10)+
  ggtitle('Overrepresented trimmed epitopes')+
  facet_grid(~Pep_len)+ylab('Epitopes')+xlab('Overrepresentation')+theme_bw()+tema_+
  guides(color=guide_legend(title = 'origin'))

plot_ = HuVi_netMHC_summ_tr_shared %>%
  mutate(overRep=ifelse(n>10,10,n),
         Pep_len=nchar(Peptide_trimm)) %>%
  group_by(overRep, Pep_len, carac) %>% summarise(n=n()) %>%
  group_by(Pep_len, carac) %>% mutate(prop=prop.table(n)) 

p4=plot_ %>%
  ggplot(aes(x=overRep, y=prop, color=carac))+
  geom_point()+geom_line()+
  scale_x_continuous(breaks = 1:10)+facet_grid(~Pep_len)+
  ggtitle('Shared trimmed epitopes')+
  ylab('Proportion of Epitopes')+xlab('Overrepresentation')+
  scale_color_manual(values = c('black','gray'))+theme_bw()+tema_+
  guides(color=guide_legend(title = 'epitope'))

png('_FigureS5.png', width = 1400, height = 1200, res=150)
svg('Figure_S5', width = 9)
plot_grid(
  plot_grid(plot_grid(p1, labels = 'A'), plot_grid(p2, labels = 'B'), align = 'hv', axis = 'lr', nrow = 2),
  plot_grid(plot_grid(p3, labels = 'C'), plot_grid(p4, labels = 'D'), nrow = 1), nrow = 2, rel_heights = c(2,1))
dev.off()
#####


# ------------
# After revision----
library(Biostrings)

# Sequence identity scores of epitopes recognizing cross reactive CDR3 (databases)----
Epitopes_recurrent_CDR3_paired_l = recurrent_CDR3_paired %>%
    group_by(n_Epitopes, CDR3_beta, CDR3_alpha) %>% summarise(Epitopes=toString(Epitope)) %>%
  pull(Epitopes) %>% strsplit(., ', ')

PercentIdentity <- lapply(Epitopes_recurrent_CDR3_paired_l, function(epitopes_){
  if(length(epitopes_)>2){
    sapply(epitopes_, function(x) sapply(epitopes_, function(y) pairwiseAlignment(x, y) %>% pid(.)))
  }else{
    pairwiseAlignment(epitopes_[1], epitopes_[2]) %>% pid(.)
  }
})

PercentIdentity_ = PercentIdentity
# Remove repeated values when comparing >2 epitopes
ind_ = lapply(PercentIdentity_, length) > 1
PercentIdentity_[ind_] = lapply(PercentIdentity_[ind_], function(x) x[upper.tri(x)])
PercentIdentity_[[300]]

PercentIdentity_ = unlist(PercentIdentity_)
p1_2=ggplot()+geom_histogram(aes(x=PercentIdentity_))+
  xlab('Identity score of CR epitopes')+
  theme_bw()+
  theme(axis.text = element_text(size=13, color = 'black'), axis.title = element_text(size=14, color = 'black'), 
        plot.title = element_text(hjust = .5))+
  scale_y_continuous(expand = expansion(mult = c(0, .1)))+
  scale_x_continuous(breaks = signif(0:9/9*100, 2))
#  scale_x_continuous(breaks = c(0,12,20,33,50,75,88))


# OK: After review: Figure 1
p1=recurrent_CDR3_paired %>%
  group_by(n_Epitopes, CDR3_beta, CDR3_alpha) %>% summarise() %>%
  ggplot(aes(x=n_Epitopes))+geom_bar()+
  ggtitle('Recurrent TRA:TRB')+ylab('CDR3 a:b pair')+xlab('Epitopes recognized')+
  scale_x_continuous(breaks = c(2:6,10,12))+scale_y_sqrt(breaks = c(1,10,200,400,600), expand = expansion(mult = c(0, .1)))+
  theme_bw()+
  theme(axis.text = element_text(size=13, color = 'black'), axis.title = element_text(size=14, color = 'black'), 
        plot.title = element_text(hjust = .5))

peps_12 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==12) %>% pull(Epitope)
peps_10 = dplyr::filter(recurrent_CDR3_paired, n_Epitopes==10) %>% pull(Epitope)

t1=ggdraw(xlim = c(.1,1))+
  draw_label(paste0(peps_12, collapse = '\n'), fontfamily = 'mono', x = .3, size=15)+
  draw_label(paste0(c('','',peps_10), collapse = '\n'), fontfamily = 'mono', x = .7, size=15)+
  geom_rect(aes(xmin=.15, xmax=.45, ymin=.05, ymax=.95), fill='transparent', color='black')+
  geom_rect(aes(xmin=.55, xmax=.85, ymin=.05, ymax=.8), fill='transparent', color='black')+
  geom_segment(aes(x=.45, xend=.9, y=.9, yend=.9), arrow = arrow())+
  geom_segment(aes(x=.85, xend=1, y=.4, yend=.3), arrow = arrow())


t1=ggdraw(xlim = c(.1,.9))+
  draw_label(paste0(peps_12, collapse = '\n'), fontfamily = 'mono', x = .3, size=15)+
  draw_label(paste0(c('','',peps_10), collapse = '\n'), fontfamily = 'mono', x = .65, size=15)+
  geom_rect(aes(xmin=.15, xmax=.45, ymin=.05, ymax=.95), fill='transparent', color='black')+
  geom_rect(aes(xmin=.5, xmax=.8, ymin=.05, ymax=.8), fill='transparent', color='black')+
  geom_segment(aes(x=.45, xend=.9, y=.9, yend=.9), arrow = arrow())+
  geom_segment(aes(x=.8, xend=.9, y=.4, yend=.3), arrow = arrow())

#t1

p2=ggplot()+geom_logo(peps_12)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab(NULL)+
  theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'bottom', legend.title = element_blank())+
  guides(fill=guide_legend(nrow = 2))

p3=ggplot()+geom_logo(peps_10)+theme_bw()+scale_y_continuous(expand = c(0,0))+xlab('Position')+
  theme(axis.text = element_text(size=11), axis.title = element_text(size=12), legend.position = 'none')

plot_grid(p1, t1,
          plot_grid(p2, p3, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(1.5,1)),
          NULL, p1_2,
          nrow = 1, rel_widths = c(.7,1,.7,.1,.7))

#
tema_fig1_=list(theme(axis.title = element_text(size=14, color = 'black'), legend.text = element_text(size=12)),
                guides(color=guide_legend(override.aes = list(size=2))))
p4=plot_grid(
  plot_CR(dist_all_Images_melted$Top, ind_12, 'canberra', 'pHLA Top images')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_pHLA_filt_melted, ind_12, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_epitope_melted, ind_12, ggtitle_ = 'Epitope (tr_c+XYZ)')+tema_fig1_,
  nrow = 1, rel_widths = c(1,.1,1,.1,1.25)
)

p5=plot_grid(
  plot_CR(dist_all_Images_melted$Top, ind_10, 'canberra', 'pHLA Top images')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_pHLA_filt_melted, ind_10, ggtitle_ = 'ec pHLA (tr_c+XYZ)')+tema_fig1_+NoLegend(),NULL,
  plot_CR(dist_all_epitope_melted, ind_10, ggtitle_ = 'Epitope (tr_c+XYZ)')+tema_fig1_, 
  nrow = 1, rel_widths = c(1,.1,1,.1,1.25)
)

# plot_grid(plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(1,1,1.25)),
#           plot_grid(p4, p5, p6, nrow = 1, rel_widths = c(1,1,1.25)), nrow=2)

# final
png('_afterReview_Figure1.png', width = 1850, height = 1400, res = 150)
plot_grid(plot_grid(p1, t1,
                    plot_grid(p2, p3, nrow = 2, align = 'hv', axis = 'lr', rel_heights = c(1.5,1)),
                    NULL, p1_2,
                    nrow = 1, rel_widths = c(.6,.8,.65,.05,.75)),
          NULL,
          p4, 
          NULL,
          p5, nrow = 5, rel_heights = c(1,.1,1,.1,1))

dev.off()


# Question of the reviewer about Figure S5
p1=HuVi_netMHC_summ %>%
  mutate(n=ifelse(n>=30, 30, ifelse(n>=20,20, ifelse(n>=15, 15, ifelse(n>=10, 10, n))))) %>%
  group_by(Pep_len, source, n) %>% summarise(overRep=n()) %>%
  ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
  facet_grid(~Pep_len)+scale_y_log10(labels = scales::label_comma())+
  scale_x_sqrt(breaks = c(1:10, 15, 20, 30))+
  ggtitle('Figure S5Ax: Overrepresented epitopes')+
  ylab('Epitopes')+xlab('Overrepresentation')+
  theme_bw()+tema_+guides(color=guide_legend(title = 'origin'))

p2=HuVi_netMHC_summ %>%
  group_by(Pep_len, source, n) %>% summarise(overRep=n()) %>%
  ggplot(aes(x=n, y=overRep, color=source, group=source))+geom_point()+geom_line()+
  facet_grid(~Pep_len)+scale_y_log10(labels = scales::label_comma())+
  scale_x_continuous(limits = c(0,50))+
  ggtitle('Figure S5Ay: Overrepresented epitopes')+
  ylab('Epitopes')+xlab('Overrepresentation')+
  theme_bw()+tema_+guides(color=guide_legend(title = 'origin'))

plot_grid(p1, p2, nrow = 2)


# Question of the reviewer about Figure 5
bind_rows(dist_selfViral_Images_melted$Top %>% ungroup() %>%
            # get the top 300, distant epitopes between them
            arrange(value) %>% top_n(300, wt = -value),
          dist_selfViral_Images_melted$Top %>% ungroup() %>%
            # get the 300, closer epitopes between them
            arrange(value) %>% top_n(300, wt = value)) %>%
  ggplot(aes(x=reorder(ID, value), y=value, color=selfViral))+geom_point(size=.7)+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size=12, color = 'black'),
        axis.line = element_line(), axis.title = element_text(size=14, color = 'black'), 
        axis.title.x = element_text(hjust = .6),
        legend.title = element_blank(), legend.key = element_blank(),
        legend.text = element_text(size=12), plot.title = element_text(hjust=.6, size=14))+
  ylab(paste0('\nDistance between images','\n(method=canberra)'))+xlab('Ranked epitope pairs')+
  #ggtitle(ggtitle_)+
  scale_y_continuous(labels = scales::label_comma(), sec.axis = sec_axis(~ ., breaks = NULL, labels = NULL))+
  scale_color_manual(values = c('black','red','gray'))+
  guides(color=guide_legend(override.aes = list(size=1.5)))


# <<<<<<<<<<<<
# Questions of the reviewer about consistent features of tumor antigens that ressemble self antigens
head(dist_selfViral_Images_melted$Top)
head(CancerViralCR_Top_noWT)

plot_CR_dotplot(CancerViralCR_Top_noWT, dist_selfViral_ = dist_selfViral_Images_melted$Top, 
                ggtitle_ = 'pHLA Top perspective') # , annot_ = CancerViralCR_Top_$Epitope

CancerViralCR_Top_noWT_dcast <- CancerViralCR_Top_noWT %>%
  reshape2::dcast(Epitope~selfViral, value.var = 'value_min')

p1=CancerViralCR_Top_noWT_dcast %>%
  # selecting Tumor antigens closest to self
  dplyr::filter(`Human/Human` < `Human/Viral`) %>% pull(Epitope) %>%
  gsub('^(\\w)\\w(\\w{6})\\w','\\1 \\2 ',.) %>%
  ggseqlogo()+labs(title='CloserToHuman')

p2=CancerViralCR_Top_noWT_dcast %>%
  # selecting Tumor antigens closest to self
  dplyr::filter(`Human/Human` > `Human/Viral`) %>% pull(Epitope) %>%
  gsub('^(\\w)\\w(\\w{6})\\w','\\1 \\2 ',.) %>%
  ggseqlogo()+labs(title='CloserToViral')+guides(fill='none')


plot_grid(p1, p2, align = 'h')

## Atchley scores
epitope_ = 'YMDGTMSQV'
annot_ATCHLEY <- function(epitope_){
  tibble(amino.acid=strsplit(epitope_, '') %>% unlist(),
         pos=1:9) %>% left_join(ATCHLEY, by = 'amino.acid') 
}

TumorCloserHuman_atchley_l <- CancerViralCR_Top_noWT_dcast %>%
  # selecting Tumor antigens closest to self
  dplyr::filter(`Human/Human` < `Human/Viral`) %>% pull(Epitope) %>% 
  lapply(., annot_ATCHLEY)

TumorCloserViral_atchley_l <- CancerViralCR_Top_noWT_dcast %>%
  # selecting Tumor antigens closest to self
  dplyr::filter(`Human/Human` > `Human/Viral`) %>% pull(Epitope) %>% 
  lapply(., annot_ATCHLEY)

TumorCloser_achley_l <- bind_rows(
  bind_rows(TumorCloserHuman_atchley_l)[-1] %>% mutate(class='CloserToHuman'),
  bind_rows(TumorCloserViral_atchley_l)[-1] %>% mutate(class='CloserToViral')
)


# Calculate statistics
library(rstatix)
TumorCloser_achley_l %>% group_by(pos) %>% wilcox_test(formula = f1~class) # P6, p=0.0255
TumorCloser_achley_l %>% group_by(pos) %>% wilcox_test(formula = f2~class) # P8, p=0.0372
TumorCloser_achley_l %>% group_by(pos) %>% wilcox_test(formula = f3~class)
TumorCloser_achley_l %>% group_by(pos) %>% wilcox_test(formula = f4~class)
TumorCloser_achley_l %>% group_by(pos) %>% wilcox_test(formula = f5~class)

dplyr::filter(TumorCloser_achley_l, pos==6) %>% group_by(class, pos) %>% summarise(median(f1))
dplyr::filter(TumorCloser_achley_l, pos==8) %>% group_by(class, pos) %>% summarise(median(f2))

# Plot heatmap
quantiles_ = TumorCloser_achley_l %>% group_by(pos, class) %>% 
  summarise(f1=median(f1), f2=median(f2), f3=median(f3), f4=median(f4), f5=median(f5)) %>%
  ungroup() %>% dplyr::select(3:7) %>% unlist() %>% quantile()


# Center zero in a correlation pheatmap
cols = colorRampPalette(c('blue','white','red'))
cols = cols(100)
paletteLength = length(cols)
myBreaks <- c(seq(min(quantiles_, na.rm=TRUE), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(quantiles_, na.rm=TRUE)/paletteLength, max(quantiles_, na.rm=TRUE), length.out=floor(paletteLength/2)))


plot_atchley_heat <- function(f_){
  p=TumorCloser_achley_l %>% 
    reshape2::dcast(class~pos, value.var = f_, median) %>% column_to_rownames('class') %>% 
    pheatmap(cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 0, silent = TRUE, main = f_,
             color = cols, breaks = myBreaks)
  p[[4]]
}

sapply(paste0('f', 1:5), plot_atchley_heat) %>% plot_grid(plotlist = ., ncol = 1)

# Confirm polarity in P6 for Factor 1
bind_rows(
  CancerViralCR_Top_noWT_dcast %>%
    dplyr::filter(`Human/Human` < `Human/Viral`) %>% 
    dplyr::select(Epitope) %>% mutate(class='CloserToHuman'),
  CancerViralCR_Top_noWT_dcast %>%
    dplyr::filter(`Human/Human` > `Human/Viral`) %>% 
    dplyr::select(Epitope) %>% mutate(class='CloserToViral')
) %>%
  mutate(P6=substr(Epitope, 6,6),
         chemistry=ifelse(P6 %in% c('A','F','I','L','M','P','V','W'), 'Hydrophobic',
                          ifelse(P6 %in% c('C','G','S','T','Y'), 'Polar',
                                 ifelse(P6 %in% c('D','E'), 'Acid',
                                        ifelse(P6 %in% c('H','K','R'), 'Basic',
                                               ifelse(P6 %in% c('N','Q'), 'Neutral', 'error')))))) %>%
  group_by(class, chemistry) %>% summarise(n=n()) %>%
  mutate(prop=signif(prop.table(n)*100, 2)) %>% formattable::formattable()

# ggseqlogo(data = paste0(LETTERS, collapse = ''))


# Electrostatic profile stats


bind_rows(
  CancerViralCR_Top_noWT_dcast %>%
    dplyr::filter(`Human/Human` < `Human/Viral`) %>% 
    mutate(class='CloserToHuman',
           min_dist=`Human/Human`) %>%
    dplyr::select(Epitope, class, min_dist),
  CancerViralCR_Top_noWT_dcast %>%
    dplyr::filter(`Human/Human` > `Human/Viral`) %>% 
    mutate(class='CloserToViral',
           min_dist=`Human/Viral`) %>%
    dplyr::select(Epitope, class, min_dist)
) %>%
  group_by(class) %>% summarise(median=median(min_dist))

#####


# Supplementary Table (xls file) with proposed cross-reactivities

format_table_ <- function(df_){
  as_tibble(df_) %>% mutate(ID=strsplit(ID, '_') %>% lapply(., function(x) paste0(sort(x), collapse = '_')) %>% unlist) %>%
    dplyr::filter(!duplicated(ID)) %>%
    dplyr::select(Epitope1=Epitope, Epitope2=Var2, distance=value, origin=selfViral) %>%
    unique() %>% arrange(distance)
}

supp_a <- format_table_(CancerViralCR_epitope) %>% dplyr::rename(euclidean_dist=distance)
supp_b <- format_table_(CancerViralCR_ecpHLA) %>% dplyr::rename(euclidean_dist=distance)
supp_c <- format_table_(CancerViralCR_Top) %>% dplyr::rename(canberra_dist=distance)

supp_d <- dplyr::rename(AutoimmuneSelfCR, Var2=Peptide) %>% format_table_() %>% dplyr::rename(canberra_dist=distance)

#fwrite(supp_a, 'results/Supplementary_Data_3A.csv.gz')
#fwrite(supp_b, 'results/Supplementary_Data_3B.csv.gz')


head(supp_a) %>% reshape2::dcast(Epitope1~Epitope2, value.var='euclidean_dist')

