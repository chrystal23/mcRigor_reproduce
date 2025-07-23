
library(Seurat)
library(mcRigor)

data_path <- './zman_data/'


###### Data loading

filenames = list.files(paste0(data_path, 'raw/'), pattern = '*.txt.gz')
counts = matrix(0, nrow = 52634, ncol=0)
for (ii in c(1:56)) {
  ff = filenames[grep(paste0('GSM73107', ifelse(1+ii>9, 1+ii, paste0(0,1+ii))), filenames)]
  matt = read.table(gzfile(paste0('/raw/', ff)))  
  counts = cbind(counts, matt)
}
counts = as.matrix(counts)

metadata = read.table(gzfile(paste0(data_path, 'raw/','GSE232040_metadata_q_annotated.txt.gz')),
                      fill = T, nrows = 3)
metanames = c('id', colnames(metadata),'xx')
metadata = read.table(gzfile(paste0(data_path, 'raw/','GSE232040_metadata_q_annotated.txt.gz')),
                      col.names = metanames,
                      fill = T, nrows = 21505) 
metadata$Well_ID[dim(metadata)[1]] == colnames(counts)[dim(counts)[2]]

metadata = metadata[-c(1),]
metadata = metadata[,-c(1)]

for (jj in 1:dim(metadata)[1]){
  if (metadata$xx[jj] != '') {
    cat(jj, metadata$cluster_colors[jj], metadata$sc_x[jj], '\n')
    metadata$cluster_colors[jj] = paste0(metadata$cluster_colors[jj], ' ', metadata$sc_x[jj])
    metadata[jj, c(13:15)] = metadata[jj, c(14:16)]
  }
}

metadata = metadata[,-c(16)]
all(colnames(counts) == metadata$Well_ID)
obj_singlecell_zman_all = CreateSeuratObject(counts = counts, 
                                             meta.data = metadata)

load(paste0(data_path, "mc_example/mc.T_clean.Rda"))
obj_singlecell_zman = obj_singlecell_zman_all[, names(object@mc)]
load(paste0(data_path, "mc_example/mat.T_clean.Rda"))
obj_singlecell_zman = obj_singlecell_zman[object@genes,]

obj_singlecell_zman$celltype = obj_singlecell_zman$cluster_colors
obj_singlecell_zman = NormalizeData(obj_singlecell_zman)
obj_singlecell_zman = FindVariableFeatures(obj_singlecell_zman)
obj_singlecell_zman = ScaleData(obj_singlecell_zman)
obj_singlecell_zman = RunPCA(obj_singlecell_zman)
obj_singlecell_zman = RunUMAP(obj_singlecell_zman, reduction = 'pca', dims = 1:50, seed.use = 123456)
UMAPPlot(obj_singlecell_zman, group.by = 'celltype')

saveRDS(obj_singlecell_zman, file = paste0(data_path, 'zman.rds'))   # intermediate results saved

obj_singlecell_zman = readRDS(paste0(data_path, 'zman.rds'))




###### Run mcRigor

# obtained by MetaCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "mc1_cell_membership_rna_zman.csv"), check.names = F, row.names = 1)

res = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership, weight = NULL,
                       output_file = paste0(data_path, 'mc1_tabmc_RNA_zman.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_seacells = read.csv(paste0(data_path,
                                           "seacells_cell_membership_rna_zman.csv"), check.names = F, row.names = 1)

res_seacells = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_seacells, weight = NULL,
                                output_file = paste0(data_path, 'seacells_tabmc_RNA_zman.RData'))   # intermediate results saved

# obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_supercell = read.csv(paste0(data_path,
                                            "supercell_cell_membership_rna_zman.csv"), check.names = F, row.names = 1)

res_supercell = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_supercell, weight = NULL, 
                                 output_file = paste0(data_path, 'supercell_tabmc_RNA_zman.RData'))   # intermediate results saved





###### Figure 2d & Supplementary Figure 10

load(paste0(data_path, "mc1_tabmc_RNA_zman.RData"))

TabMC$depctr_TT_div = TabMC$T_rowperm1 / TabMC$T_bothperm1
TabMC$depctr_TT_div[is.na(TabMC$depctr_TT_div)] = 1
TabMC[is.na(TabMC)] = 1
rowperm_col = grep('T_rowperm', colnames(TabMC))
bothperm_col = grep('T_bothperm', colnames(TabMC))

Thre = as.data.frame(NULL)
for (size in unique(TabMC$size)) {
  if (size == 1) next
  mcs = TabMC[TabMC$size==size,]
  thre = quantile(unlist(mcs[, rowperm_col] / mcs[, bothperm_col]), probs = 1-0.01)
  Thre = rbind(Thre, c(size, thre))
}
colnames(Thre) = c('size', 'thre')
Thre = Thre[order(Thre$size),]
a=lowess(x=Thre$size, y=Thre$thre, f=1/3)
Thre$thre = a$y
Thre$size = a$x

TabMC$testres = 'trustworthy'
for (mcid in 1:dim(TabMC)[1]) {
  if (TabMC$size[mcid] %in% Thre$size || TabMC$size[mcid] == 1) {
    if (TabMC$size[mcid] >1 && TabMC$TT_div[mcid] > Thre$thre[Thre$size == TabMC$size[mcid]]) TabMC$testres[mcid] = 'dubious'
  } else {
    temp_thre = approx(x = Thre$size, y = Thre$thre, xout = TabMC$size[mcid], 
                       yleft = min(Thre$thre), yright = max(Thre$thre))$y
    if (TabMC$TT_div[mcid] > temp_thre) TabMC$testres[mcid] = 'dubious'
  }
}
print(table(TabMC$testres))

DD <- data.frame(gamma=sort(unique(TabMC$gamma)))
DD$D <- 0
DD$ZeroRate <- 0
for (gamma in DD$gamma) {
  mcs = TabMC[TabMC$gamma==gamma & TabMC$testres == 'trustworthy',]
  dub_mc = TabMC[TabMC$gamma==gamma & TabMC$testres == 'dubious',]
  DD$D[DD$gamma == gamma] = sum(dub_mc$size) / (sum(mcs$size) + sum(dub_mc$size)) 
  DD$ZeroRate[DD$gamma == gamma] = mean(TabMC[TabMC$gamma == gamma, 'ZeroRate'])
}
DD = DD[!is.nan(DD$D),]
if (max(DD$gamma) - min(DD$gamma) > 5){
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 19 / (max(DD$gamma) - min(DD$gamma)))$y
}
DD$score = DD$ZeroRate * 0.5 + DD$D * 0.5
DD$score = 1 - DD$score
opt_gamma = DD$gamma[which.max(DD$score)]

pelbow = ggplot(DD, aes(x=gamma, y=D)) + geom_path(aes(col = 'D')) +
  annotate('text', x=opt_gamma, y=-0.03, label=opt_gamma, color='red', size = 3) +
  geom_path(aes(x = gamma, y = ZeroRate, col = 'ZeroRate')) + 
  geom_path(aes(x = gamma, y = score, col = 'Score')) +
  scale_color_manual(name = NULL, values = c('Score' = 'darkred', 'D' = 'darkblue', 'ZeroRate' = 'darkgreen')) +
  geom_vline(xintercept = opt_gamma, col = 'red', linetype = 'dashed') +
  ylab(' ') +
  coord_cartesian(ylim = c(0,1), clip = 'off')+
  theme_light()+ theme(aspect.ratio = 2/3)
suppressWarnings(elbow_legend <- cowplot::get_legend(pelbow))
pelbow = pelbow + theme(legend.position = 'none')
pelbow   # Supplementary Figure 9a

sc_membership = cell_membership[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership)

mcRigor_projection(obj_singlecell = obj_singlecell_zman,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette =  c("#D81B60",  "#6D4C41", '#440154', '#31688E', '#35B779', "#9E9D24"),
                   add_testres = T, test_stats = TabMC, 
                   dub_mc_test.label = T)   # Supplementary Figure 10b

obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_zman, sc_membership = sc_membership,
                               add_testres = T, test_stats = TabMC, Thre = Thre)


#

load(paste0(data_path, "seacells_tabmc_RNA_zman.RData"))

TabMC$depctr_TT_div = TabMC$T_rowperm1 / TabMC$T_bothperm1
TabMC$depctr_TT_div[is.na(TabMC$depctr_TT_div)] = 1
TabMC[is.na(TabMC)] = 1
rowperm_col = grep('T_rowperm', colnames(TabMC))
bothperm_col = grep('T_bothperm', colnames(TabMC))

Thre = as.data.frame(NULL)
for (size in unique(TabMC$size)) {
  if (size == 1) next
  mcs = TabMC[TabMC$size==size,]
  thre = quantile(unlist(mcs[, rowperm_col] / mcs[, bothperm_col]), probs = 1-0.01)
  Thre = rbind(Thre, c(size, thre))
}
colnames(Thre) = c('size', 'thre')
Thre = Thre[order(Thre$size),]
a=lowess(x=Thre$size, y=Thre$thre, f=2/3)
Thre$thre = a$y
Thre$size = a$x

TabMC$testres = 'trustworthy'
for (mcid in 1:dim(TabMC)[1]) {
  if (TabMC$size[mcid] %in% Thre$size || TabMC$size[mcid] == 1) {
    if (TabMC$size[mcid] >1 && TabMC$TT_div[mcid] > Thre$thre[Thre$size == TabMC$size[mcid]]) TabMC$testres[mcid] = 'dubious'
  } else {
    temp_thre = approx(x = Thre$size, y = Thre$thre, xout = TabMC$size[mcid], 
                       yleft = min(Thre$thre), yright = max(Thre$thre))$y
    if (TabMC$TT_div[mcid] > temp_thre) TabMC$testres[mcid] = 'dubious'
  }
}
print(table(TabMC$testres))

DD <- data.frame(gamma=sort(unique(TabMC$gamma)))
DD$D <- 0
DD$ZeroRate <- 0
for (gamma in DD$gamma) {
  mcs = TabMC[TabMC$gamma==gamma & TabMC$testres == 'trustworthy',]
  dub_mc = TabMC[TabMC$gamma==gamma & TabMC$testres == 'dubious',]
  DD$D[DD$gamma == gamma] = sum(dub_mc$size) / (sum(mcs$size) + sum(dub_mc$size)) 
  DD$ZeroRate[DD$gamma == gamma] = mean(TabMC[TabMC$gamma == gamma, 'ZeroRate'])
}
DD = DD[!is.nan(DD$D),]
if (max(DD$gamma) - min(DD$gamma) > 5){
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 15 / (max(DD$gamma) - min(DD$gamma)))$y
}
DD$score = DD$ZeroRate * 0.5 + DD$D * 0.5
DD$score = 1 - DD$score
opt_gamma = DD$gamma[which.max(DD$score)]

pelbow = ggplot(DD, aes(x=gamma, y=D)) + geom_path(aes(col = 'D')) +
  annotate('text', x=opt_gamma, y=-0.03, label=opt_gamma, color='red', size = 3) +
  geom_path(aes(x = gamma, y = ZeroRate, col = 'ZeroRate')) + 
  geom_path(aes(x = gamma, y = score, col = 'Score')) +
  scale_color_manual(name = NULL, values = c('Score' = 'darkred', 'D' = 'darkblue', 'ZeroRate' = 'darkgreen')) +
  geom_vline(xintercept = opt_gamma, col = 'red', linetype = 'dashed') +
  ylab(' ') +
  coord_cartesian(ylim = c(0,1), clip = 'off')+
  theme_light()+ theme(aspect.ratio = 2/3)
suppressWarnings(elbow_legend <- cowplot::get_legend(pelbow))
pelbow = pelbow + theme(legend.position = 'none')
pelbow   # Supplementary Figure 10a

sc_membership = cell_membership_seacells[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership_seacells)

mcRigor_projection(obj_singlecell = obj_singlecell_zman,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette =  c("#D81B60",  "#6D4C41", '#440154', '#31688E', '#35B779', "#9E9D24"),
                   add_testres = T, test_stats = TabMC,
                   dub_mc_test.label = T)   # Supplementary Figure 10b

obj_metacell_sea = mcRigor_buildmc(obj_singlecell = obj_singlecell_zman, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)


#
load(paste0(data_path, "supercell_tabmc_RNA_zman.RData"))

TabMC$depctr_TT_div = TabMC$T_rowperm1 / TabMC$T_bothperm1
TabMC$depctr_TT_div[is.na(TabMC$depctr_TT_div)] = 1
TabMC[is.na(TabMC)] = 1
rowperm_col = grep('T_rowperm', colnames(TabMC))
bothperm_col = grep('T_bothperm', colnames(TabMC))

Thre = as.data.frame(NULL)
for (size in unique(TabMC$size)) {
  if (size == 1) next
  mcs = TabMC[TabMC$size==size,]
  thre = quantile(unlist(mcs[, rowperm_col] / mcs[, bothperm_col]), probs = 1-0.01)
  Thre = rbind(Thre, c(size, thre))
}
colnames(Thre) = c('size', 'thre')
Thre = Thre[order(Thre$size),]
a=lowess(x=Thre$size, y=Thre$thre, f=1/6)
Thre$thre = a$y
Thre$size = a$x

TabMC$testres = 'trustworthy'
for (mcid in 1:dim(TabMC)[1]) {
  if (TabMC$size[mcid] %in% Thre$size || TabMC$size[mcid] == 1) {
    if (TabMC$size[mcid] >1 && TabMC$TT_div[mcid] > Thre$thre[Thre$size == TabMC$size[mcid]]) TabMC$testres[mcid] = 'dubious'
  } else {
    temp_thre = approx(x = Thre$size, y = Thre$thre, xout = TabMC$size[mcid], 
                       yleft = min(Thre$thre), yright = max(Thre$thre))$y
    if (TabMC$TT_div[mcid] > temp_thre) TabMC$testres[mcid] = 'dubious'
  }
}
print(table(TabMC$testres))

DD <- data.frame(gamma=sort(unique(TabMC$gamma)))
DD$D <- 0
DD$ZeroRate <- 0
for (gamma in DD$gamma) {
  mcs = TabMC[TabMC$gamma==gamma & TabMC$testres == 'trustworthy',]
  dub_mc = TabMC[TabMC$gamma==gamma & TabMC$testres == 'dubious',]
  DD$D[DD$gamma == gamma] = sum(dub_mc$size) / (sum(mcs$size) + sum(dub_mc$size)) 
  DD$ZeroRate[DD$gamma == gamma] = mean(TabMC[TabMC$gamma == gamma, 'ZeroRate'])
}
DD = DD[!is.nan(DD$D),]
if (max(DD$gamma) - min(DD$gamma) > 5){
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 14 / (max(DD$gamma) - min(DD$gamma)))$y
}
DD$score = DD$ZeroRate * 0.5 + DD$D * 0.5
DD$score = 1 - DD$score
opt_gamma = DD$gamma[which.max(DD$score)]

pelbow = ggplot(DD, aes(x=gamma, y=D)) + geom_path(aes(col = 'D')) +
  annotate('text', x=opt_gamma, y=-0.03, label=opt_gamma, color='red', size = 3) +
  geom_path(aes(x = gamma, y = ZeroRate, col = 'ZeroRate')) + 
  geom_path(aes(x = gamma, y = score, col = 'Score')) +
  scale_color_manual(name = NULL, values = c('Score' = 'darkred', 'D' = 'darkblue', 'ZeroRate' = 'darkgreen')) +
  geom_vline(xintercept = opt_gamma, col = 'red', linetype = 'dashed') +
  ylab(' ') +
  coord_cartesian(ylim = c(0,1), clip = 'off')+
  theme_light()+ theme(aspect.ratio = 2/3)
suppressWarnings(elbow_legend <- cowplot::get_legend(pelbow))
pelbow = pelbow + theme(legend.position = 'none')
pelbow   # Supplementary Figure 10a

sc_membership = cell_membership_supercell[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership_supercell)

mcRigor_projection(obj_singlecell = obj_singlecell_zman,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette =  c("#D81B60",  "#6D4C41", '#440154', '#31688E', '#35B779', "#9E9D24"),
                   add_testres = T, test_stats = TabMC, 
                   dub_mc_test.label = T)   # Supplementary Figure 10b

obj_metacell_sup = mcRigor_buildmc(obj_singlecell = obj_singlecell_zman, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)





###

# git clone https://github.com/kenxie7/ZmanR
# R
# require(devtools)
# devtools::install("ZmanR")

library(ZmanR)
library(metacell)

scdb_init(paste0(data_path, "mc_example"), force_reinit=T)
scfigs_init(paste0(data_path, "mc_example"))
org_id = "T_clean"

object_mat = scdb_mat(org_id)
object_mc = scdb_mc(org_id)
object_mc2d = scdb_mc2d(org_id)

our_id = 'rigMC'

load(file = paste0(data_path, 'our_annot.RData'))
our_color = our_annot$celltype
our_color[our_annot$celltype == 'CD4'] = '#FF2000'
our_color[our_annot$celltype == 'CD8'] = '#FFD39B'
our_color[our_annot$celltype == 'Treg'] = '#FF9000'
our_color[our_annot$celltype == 'chemotactic'] = '#00FFFF'
our_color[our_annot$celltype == 'cytotoxic'] = '#00FFFF'
our_color[our_annot$celltype == 'intermediate'] = '#7B596D'
our_color[our_annot$celltype == 'dysfunctional'] = '#7F19F3'

mcell_mc_add_color(our_id, our_color)

mcell_add_cgraph_from_mat_bknn(mat_id=org_id, 
                               gset_id = "test_feats", 
                               graph_id="test_graph",
                               K=150,
                               dsamp=F)

mcell_mc2d_force_knn(mc2d_id=our_id, mc_id=our_id, graph_id="test_graph")
mcell_mc2d_plot(mc2d_id=our_id)

###

cell_membership <- read.csv(paste0(data_path, 'seacells_cell_membership_rna_zman.csv'), check.names = F, row.names = 1)

sc_membership = cell_membership[['43']]
sc_membership = strsplit(sc_membership, split = '-')
sc_membership = sapply(sc_membership, function(x) as.integer(x[4])+1 )
names(sc_membership) = rownames(cell_membership)

our_mc = scdb_mc(our_id)

load(file = paste0(data_path, 'our_annot.RData'))

mc_cdf <- compute_mc_cdf(our_id, well_fcs_GBM_T_time, select_mcs = 1:57, mc_annotations = our_annot, 
                         time_points = c("12H","24H","36H"), time_for_auc = c(0,12,24,36))
options(repr.plot.width=21, repr.plot.height=3.5)

NK_mc_exprs = normalize_mc_exprs(our_id, mc_annotations=our_annot, mc_cdf$mc2d_auc_time, 
                                     select_celltypes = c("chemotactic", "cytotoxic", "intermediate", "dysfunctional"))
select_GO_mc <- get_GO_exp(NK_mc_exprs$select_exprs,gene_type = "gene_names", organism = "mouse", takelog=F)
filtered_GO_mc <- filter_exp(select_GO_mc, dispersion_threshold=0.05, threads = 1)

NK_smoothed_res = smooth_zman_trajectory(filtered_GO_mc, NK_mc_exprs, ref_k = 4)
NK_predicted_res = predict_expression_along_time(our_id, NK_smoothed_res, out_len = 36, loess_degree=1, loess_span = .9)

nk_corr = calculate_corr_genes(our_id, NK_smoothed_res, "spearman")
signif(nk_corr$correlation[match(c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), 
                                 names(nk_corr$correlation))],3)
signif(nk_corr$pvalues[match(c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), names(nk_corr$correlation))],3)

nk_corr$p_val_adj = p.adjust(nk_corr$pvalues, method = "BH")
mcdeg = which(nk_corr$p_val_adj<0.05)

traj_plot <- plot_smoothed_trajectory(NK_smoothed_res, mc_cdf)
gene_plot <- plot_zman_genes_heatmap(NK_predicted_res, NK_smoothed_res, 
                                     up_regulated_genes = c('Pmepa1','Tigit','Ctla2b','Ccr2','Srgn','Plac8','Ctla2a','Car2','Itgb1','Ahr','Isg20','Ski','Lamp1','Xcl1','Lgals3',
                                                            'Itga1','Cox6a2','B4galt6','Igflr1','Tcf7','Il2','Gzmc'),          # c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"),
                                     down_regulated_genes = c('Itga4','Ccl5','Ccl3','Il12rb2','Anxa6','Gzmb','Gzma','Nkg7','Il2rg','Ccl9','Tmsb4x','Prf1','Zeb2','Cd48',
                                                              'Itgb2','Serpinb9b','Cma1','Lgals1','Klrg1','S1pr5','Cx3cr1','Ifitm10'),     #  c("Ccl3","Prf1","Gzma", "Gzmb","Nkg7")
                                     k = 2)
gene_plot   # Supplementary Figure 10d

###

mc_cdf_org <- compute_mc_cdf(org_id, well_fcs_GBM_T_time, select_mcs = 1:37, mc_annotations = GBM_T_mc_annotations, time_points = c("12H","24H","36H"), time_for_auc = c(0,12,24,36))
options(repr.plot.width=21, repr.plot.height=3.5)

NK_mc_exprs_org = normalize_mc_exprs(org_id, mc_annotations=GBM_T_mc_annotations, mc_cdf_org$mc2d_auc_time, 
                                         select_celltypes = c("chemotactic", "cytotoxic", "intermediate", "dysfunctional"))
select_GO_mc_org <- get_GO_exp(NK_mc_exprs_org$select_exprs,gene_type = "gene_names", organism = "mouse", takelog=F)
filtered_GO_mc_org <- filter_exp(select_GO_mc_org, dispersion_threshold=0.05, threads = 1)

NK_smoothed_res_org = smooth_zman_trajectory(filtered_GO_mc_org, NK_mc_exprs_org, ref_k = 4)
NK_predicted_res_org = predict_expression_along_time(org_id, NK_smoothed_res_org, out_len = 36, loess_degree=1, loess_span = .9)

nk_corr_org = calculate_corr_genes(org_id, NK_smoothed_res_org, "spearman")
signif(nk_corr_org$correlation[match(c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), names(nk_corr_org$correlation))],3)
signif(nk_corr_org$pvalues[match(c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"), names(nk_corr_org$correlation))],3)

nk_corr_org$p_val_adj = p.adjust(nk_corr_org$pvalues, method = "BH")
mcdeg_org = which(nk_corr_org$p_val_adj<0.05)

traj_plot <- plot_smoothed_trajectory(NK_smoothed_res_org, mc_cdf_org)
gene_plot <- plot_zman_genes_heatmap(NK_predicted_res_org, NK_smoothed_res_org, 
                                     up_regulated_genes = c('Pmepa1','Tigit','Ctla2b','Ccr2','Srgn','Plac8','Ctla2a','Car2','Itgb1','Ahr','Isg20','Ski','Lamp1','Xcl1','Lgals3',
                                                            'Itga1','Cox6a2','B4galt6','Igflr1','Tcf7','Il2','Gzmc'),          # c("Tigit","Xcl1","Pmepa1","Igflr1", "Gzmc"),
                                     down_regulated_genes = c('Itga4','Ccl5','Ccl3','Il12rb2','Anxa6','Gzmb','Gzma','Nkg7','Il2rg','Ccl9','Tmsb4x','Prf1','Zeb2','Cd48',
                                                              'Itgb2','Serpinb9b','Cma1','Lgals1','Klrg1','S1pr5','Cx3cr1','Ifitm10'),     #  c("Ccl3","Prf1","Gzma", "Gzmb","Nkg7")
                                     k = 2)
gene_plot  # Supplementary Figure 10d

####

p1 = ggplot(mc_cdf$mc_auc_time[mc_cdf$mc_auc_time$celltype %in% c('chemotactic','cytotoxic','intermediate','dysfunctional'),],
            aes(x=time, y=cums, group=mc,color=celltype)) + 
  geom_line() + scale_color_viridis_d(direction = 1) + 
  theme(text = element_text(size = 20),legend.position = "none") + 
  ylab('% of cells') + xlab("time (hours)")+geom_abline(intercept =  0.01412,     slope= 0.02810   ,
                                                        color="red",linetype = "dashed", alpha = .5)  +
  ggtitle(paste0('area = ', c(25.85), ', cor = ', c(0.967))) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text=element_text(size=12), axis.title = element_text(size=12))
max(mc_cdf$mc_auc_time$auc) - min(mc_cdf$mc_auc_time$auc)

mc_cdf$mc_auc_time$cTET = 1 - (mc_cdf$mc_auc_time$auc - 0) / 36
max(mc_cdf$mc_auc_time$cTET) - min(mc_cdf$mc_auc_time$cTET)

p2 = ggplot(mc_cdf_org$mc_auc_time[mc_cdf_org$mc_auc_time$celltype %in% c('chemotactic','cytotoxic','intermediate','dysfunctional'),],
            aes(x=time, y=cums, group=mc,color=celltype)) +
  geom_line() + scale_color_viridis_d(direction = 1) +
  theme(text = element_text(size = 20), legend.title = element_blank()) + 
  ylab('% of cells') + xlab("time (hours)")+geom_abline(intercept =  0.01412,     slope= 0.02810   ,
                                                        color="red",linetype = "dashed", alpha = .5) +
  ggtitle(paste0('area = ', c(19.38), ', cor = ', c(0.954))) +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        axis.text=element_text(size=12), axis.title = element_text(size=12))
max(mc_cdf_org$mc_auc_time$auc) - min(mc_cdf_org$mc_auc_time$auc)

mc_cdf_org$mc_auc_time$cTET = 1 - (mc_cdf_org$mc_auc_time$auc -0) / 36
max(mc_cdf_org$mc_auc_time$cTET) - min(mc_cdf_org$mc_auc_time$cTET)

p1+p2   # Figure 2d, top

###

obj_singlecell_zman[['Gating']] = object_mat@cell_metadata[["Gating"]][match(colnames(obj_singlecell_zman), rownames(object_mat@cell_metadata))]
Idents(obj_singlecell_zman) = 'time_assignment'
markers = FindMarkers(obj_singlecell_zman, ident.1 = c('36H'), ident.2 = '12H')
markers_genes = rownames(markers[markers$p_val_adj<0.05,])
markers_genes = markers_genes[markers_genes %in% colnames(NK_predicted_res)]

pdat = data.frame(cor = c(signif(nk_corr$correlation[match(markers_genes, names(nk_corr$correlation))],3), 
                          signif(nk_corr_org$correlation[match(markers_genes, names(nk_corr_org$correlation))],3)) %>% abs,
                  pval = c(signif(nk_corr$pvalues[match(markers_genes, names(nk_corr$correlation))],3),
                           signif(nk_corr_org$pvalues[match(markers_genes, names(nk_corr_org$correlation))],3)),
                  group = c(rep('mcRigor', length(markers_genes)), rep('original', length(markers_genes))))
pbox1 = ggplot(pdat, aes(y = pval, x = group)) +geom_boxplot(aes(fill = group)) + 
  xlab('metacell assignment') + ylab('p-value') + theme(legend.position = "none")
pbox2 = ggplot(pdat, aes(y = cor, x = group)) +geom_boxplot(aes(fill = group)) + 
  xlab('metacell assignment') + ylab('correlation') + theme(legend.position = "none")
pbox1 + pbox2   # Supplementary Figure 10c

#
markers = markers[markers$p_val_adj<0.05,]

cor_our = nk_corr$correlation[match(rownames(markers), names(nk_corr$correlation))]
pval_our = nk_corr$pvalues[match(rownames(markers), names(nk_corr$pvalues))]
cor.test(pval_our, markers$p_val_adj, method = 'spearman')

cor_org = nk_corr_org$correlation[match(rownames(markers), names(nk_corr_org$correlation))]
pval_org = nk_corr_org$pvalues[match(rownames(markers), names(nk_corr_org$pvalues))]
cor.test(pval_org, markers$p_val_adj, method = 'spearman')

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             # axis.ticks = element_blank(), axis.text = element_blank(),
             aspect.ratio = 2.6/3)  #

ggplot(data = data.frame(y = log2(pval_our), x = log2(markers$p_val_adj)),
       mapping = aes(x=x, y=y)) +
  geom_point() + xlim(c(-50, 0))

ggplot(data = data.frame(y = log2(pval_org), x = log2(markers$p_val_adj)),
       mapping = aes(x=x, y=y))+
  geom_point() + xlim(c(-50, 0))


gene_plot <- plot_zman_genes_heatmap(NK_predicted_res, NK_smoothed_res, 
                                         up_regulated_genes = markers_genes[1:50],
                                         down_regulated_genes = markers_genes[51:75],
                                         k = 2,
                                         title = paste0('our spearman correlation between metacell and single cell results: ',
                                                        0.517))
gene_plot_org <- plot_zman_genes_heatmap(NK_predicted_res_org, NK_smoothed_res_org, 
                                             up_regulated_genes = markers_genes[1:50],
                                             down_regulated_genes = markers_genes[51:75],
                                             k = 2,
                                             row_ord = hord,
                                             title = paste0('original spearman correlation between metacell and single cell results: ',
                                                            0.452))
gene_plot + gene_plot_org   # Figure 2d, bottom



