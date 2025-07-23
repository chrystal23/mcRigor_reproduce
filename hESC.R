
library(Seurat)
library(mcRigor)

data_path <- './hESC_data/'


###### Data loading

counts = read.csv(file = paste0(data_path, 'GSE75748_sc_cell_type_ec.csv'), row.names = 1)

counts_bulk = read.csv(file = paste0(data_path, 'GSE75748_bulk_cell_type_ec.csv'), row.names = 1)

metadata_bulk = data.frame(celltype = colnames(counts_bulk))
a = strsplit(metadata_bulk$celltype, split = '_')
metadata_bulk$celltype = sapply(a, function(x) x[1])
rownames(metadata_bulk) = colnames(counts_bulk)

obj_singlecell_hESC_bulk_all = CreateSeuratObject(counts = counts_bulk,
                                                  meta.data = metadata_bulk)

obj_singlecell_hESC_bulk = obj_singlecell_hESC_bulk_all[, obj_singlecell_hESC_bulk_all$celltype %in% c('H1', 'DEC')]
obj_singlecell_hESC_bulk = NormalizeData(obj_singlecell_hESC_bulk)

saveRDS(obj_singlecell_hESC_bulk, file = paste0(data_path, 'hESC_bulk.rds'))   # intermediate results saved

#
metadata = data.frame(celltype = colnames(counts))
a = strsplit(metadata$celltype, split = '_')
metadata$celltype = sapply(a, function(x) x[1])
metadata$batch = sapply(a, function(x) x[2])
a = strsplit(metadata$batch, split = '\\.')
metadata$batch = sapply(a, function(x) x[1])

obj_singlecell_hESC_all = CreateSeuratObject(counts = counts, meta.data = metadata)

obj_singlecell_hESC = obj_singlecell_hESC_all[, obj_singlecell_hESC_all$celltype %in% c('H1', 'DEC')]

obj_singlecell_hESC = NormalizeData(obj_singlecell_hESC)
obj_singlecell_hESC = FindVariableFeatures(obj_singlecell_hESC)
obj_singlecell_hESC = ScaleData(obj_singlecell_hESC)
obj_singlecell_hESC = RunPCA(obj_singlecell_hESC)
obj_singlecell_hESC = RunUMAP(obj_singlecell_hESC, reduction = 'pca', dims = 1:30, seed.use = 123456)
UMAPPlot(obj_singlecell_hESC, group.by = 'celltype')
UMAPPlot(obj_singlecell_hESC, group.by = 'batch')

saveRDS(obj_singlecell_hESC, file = paste0(data_path, 'hESC.rds'))   # intermediate results saved

obj_singlecell_hESC = readRDS(file = paste0(data_path, 'hESC.rds'))



###### Run mcRigor

# obtained by MetaCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "mc1_cell_membership_rna_hESC.csv"), check.names = F, row.names = 1)

res = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership, weight = NULL,
                       output_file = paste0(data_path, 'mc1_tabmc_RNA_hESC.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_seacells = read.csv(paste0(data_path,
                                           "seacells_cell_membership_rna_hESC.csv"), check.names = F, row.names = 1)

res_seacells = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_seacells, weight = NULL,
                                output_file = paste0(data_path, 'seacells_tabmc_RNA_hESC.RData'))   # intermediate results saved

# obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_supercell = read.csv(paste0(data_path,
                                            "supercell_cell_membership_rna_hESC.csv"), check.names = F, row.names = 1)

res_supercell = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_supercell, weight = NULL, 
                                 output_file = paste0(data_path, 'supercell_tabmc_RNA_hESC.RData'))   # intermediate results saved




###### Supplementary Figure 9

load(paste0(data_path, "mc1_tabmc_RNA_hESC.RData"))

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
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 30 / (max(DD$gamma) - min(DD$gamma)))$y
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

mcRigor_projection(obj_singlecell = obj_singlecell_hESC,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette = c("#D4E157","#D81B60"),
                   add_testres = T, test_stats = TabMC, Thre = Thre,
                   dub_mc_test.label = T)   # Supplementary Figure 9b

obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_hESC, sc_membership = sc_membership,
                               add_testres = T, test_stats = TabMC, Thre = Thre)


#

load(paste0(data_path, "seacells_tabmc_RNA_hESC.RData"))

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
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 7 / (max(DD$gamma) - min(DD$gamma)))$y
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

sc_membership = cell_membership_seacells[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership_seacells)

mcRigor_projection(obj_singlecell = obj_singlecell_hESC,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette = c("#D4E157","#D81B60"),
                   add_testres = T, test_stats = TabMC, 
                   dub_mc_test.label = T)   # Supplementary Figure 9b

obj_metacell_sea = mcRigor_buildmc(obj_singlecell = obj_singlecell_hESC, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)

#
load(paste0(data_path, "supercell_tabmc_RNA_hESC.RData"))

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
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 10 / (max(DD$gamma) - min(DD$gamma)))$y
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

sc_membership = cell_membership_supercell[[opt_gamma]]
names(sc_membership) = rownames(cell_membership_supercell)

mcRigor_projection(obj_singlecell = obj_singlecell_hESC,
                   sc_membership = sc_membership,
                   color_field = 'celltype',
                   cpalette = c("#D4E157","#D81B60"),
                   add_testres = T, test_stats = TabMC, 
                   dub_mc_test.label = T)   # Supplementary Figure 9b

obj_metacell_sup = mcRigor_buildmc(obj_singlecell = obj_singlecell_hESC, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)




###### Figure 2c

obj_bulk_hESC = readRDS(file = paste0(data_path, 'hESC_bulk.rds'))

Idents(obj_bulk_hESC) = 'celltype'

p_thre = 0.05

DEgenes_MAST = FindMarkers(obj_bulk_hESC, 
                           ident.1 = 'H1', ident.2 = 'DEC', test.use = 'MAST', min.cells.group = 2)
degenes_MAST = rownames(DEgenes_MAST)[DEgenes_MAST$p_val_adj < p_thre]


##
cell_membership <- read.csv(paste0(data_path, 'mc1_cell_membership_rna_hESC.csv'), check.names = F, row.names = 1)

Gammas = c(colnames(cell_membership), '1')

ResDE_mc1 = list()

for (tgamma in Gammas) {
  
  if (tgamma == '1') obj_metacell = obj_singlecell_hESC
  else {
    sc_membership = cell_membership[[tgamma]]
    names(sc_membership) = rownames(cell_membership)
    
    sc_membership = sc_membership[colnames(obj_singlecell_hESC)]
    
    obj_metacell = rigorMC_buildmc(obj_singlecell_hESC, sc_membership, 
                                   aggregate_method = 'sum')
  }
  
  Idents(obj_metacell) = 'celltype'
  mc_DEgenes_MAST = FindMarkers(obj_metacell,
                                ident.1 = 'H1', ident.2 = 'DEC', test.use = 'MAST', min.cells.group = 1)
  mc_deg_MAST = rownames(mc_DEgenes_MAST)[mc_DEgenes_MAST$p_val_adj < p_thre]
  
  gtp_MAST = intersect(mc_deg_MAST, degenes_MAST)
  gfp_MAST = setdiff(mc_deg_MAST, degenes_MAST)
  gfn_MAST = setdiff(degenes_MAST, mc_deg_MAST)
  
  gprecision_MAST = length(gtp_MAST) / (length(gtp_MAST) + length(gfp_MAST))
  grecall_MAST = length(gtp_MAST) / (length(gtp_MAST) + length(gfn_MAST))
  gFF_MAST = 2 * gprecision_MAST * grecall_MAST / (gprecision_MAST + grecall_MAST)
  
  ResDE_mc1[[tgamma]] = list(gprecision_MAST = gprecision_MAST, grecall_MAST = grecall_MAST, gFF_MAST = gFF_MAST,
                             mc_DEgenes_MAST = mc_DEgenes_MAST, mc_deg_MAST = mc_deg_MAST)
}

save(ResDE_mc1, file = paste0(data_path, 'ResDE_mc1.RData'))   # intermediate results saved, similar for seacells and supercell


cut_bins <- function(x, bin_width = 5){
  
  x$gamma = as.integer(x$gamma)
  metrics = setdiff(colnames(x), 'gamma')
  
  lwb = floor(min(x$gamma) / bin_width) * bin_width
  upb = ceiling(max(x$gamma) / bin_width) * bin_width
  bins = cut(x$gamma, breaks = seq(lwb, upb, by = bin_width))
  
  x = x[!is.na(bins),]
  
  xlist = lapply(levels(bins), function(ll) {
    x$gamma[bins == ll]
  })
  
  psg = seq(lwb+1, upb, by=1)
  psbins = cut(psg, breaks = seq(lwb, upb, by = bin_width))
  pslist = lapply(levels(psbins), function(ll) {
    psg[psbins == ll]
  })
  
  newx = data.frame(gamma = sapply(pslist, mean),
                    bin = levels(bins))
  for (mm in metrics){
    newx[[mm]] = sapply(xlist, function(y) mean(x[[mm]][x$gamma %in% y]))
  }
  
  return(newx)
}


##
sc_ff = c(1, NA, ResDE_seacells[['1']]$gFF_MAST)
sea = sapply(ResDE_seacells, function(x) x$gFF_MAST)
sup = sapply(ResDE_supercell, function(x) x$gFF_MAST)
mmc = sapply(ResDE_mc1, function(x) x$gFF_MAST)
pdat_sea = data.frame(gamma = as.integer(names(ResDE_seacells)),
                      fscore = sea)
pdat_sea = cut_bins(pdat_sea)
pdat_sea = rbind(pdat_sea, sc_ff)
pdat_mmc = data.frame(gamma = as.integer(names(ResDE_mc1)),
                      fscore = mmc)
pdat_mmc = cut_bins(pdat_mmc)
pdat_mmc = pdat_mmc[pdat_mmc$gamma<50,]
pdat_mmc = rbind(pdat_mmc, sc_ff)
pdat_sup = data.frame(gamma = as.integer(names(ResDE_supercell)),
                      fscore = sup)
pdat_sup = cut_bins(pdat_sup)
pdat_sup = rbind(pdat_sup, sc_ff)

ppff = ggplot(pdat_sea, aes(x = gamma, y = fscore)) + 
  geom_line(aes(col = 'SEACells'), linetype = 'longdash') + geom_point(aes(col = 'SEACells')) +
  annotate('point', x = 13, y = -Inf, size = 3, pch = 17, color = 'purple') +
  geom_line(data = pdat_sup, aes(x = gamma, y = fscore, col = 'SuperCell'), linetype = 'longdash') +
  geom_point(data = pdat_sup, aes(x = gamma, y = fscore, col = 'SuperCell')) +
  annotate('point', x = 4, y = -Inf, size = 3, pch = 17, color = 'darkolivegreen4') +
  geom_line(data = pdat_mmc, aes(x = gamma, y = fscore, col = 'MetaCell'), linetype = 'longdash') +
  geom_point(data = pdat_mmc, aes(x = gamma, y = fscore, col = 'MetaCell')) +
  annotate('point', x = 6, y = -Inf, size = 3, pch = 17, color = 'blue') +
  scale_color_manual(name="",values=c('SEACells'='purple', 'SuperCell'='darkolivegreen4', 'MetaCell'='blue')) +
  geom_vline(xintercept = 13, linetype = 'dashed', col = 'red') +
  annotate('text', label = '13', x=13, y = 0, col='red', size = 4) +
  coord_cartesian(ylim = c(0, 0.45), xlim = c(1,50), clip = 'off')
ppff + theme(legend.position = 'none')   # Figure 2c, bottom left


##
library(ggvenn)

venn_dat = list(bulk = degenes_MAST,
                SEACells = ResDE[['12']]$mc_deg_MAST,
                SuperCell = ResDE_supercell[['4']]$mc_deg_MAST,
                MetaCell = ResDE_mc1[['6']]$mc_deg_MAST)

ggvenn(venn_dat, show_stats = 'c',
       fill_color = c("red", 'purple', 'darkolivegreen4', "blue"),)   # Figure 2c, bottom right



###
library(ComplexHeatmap)

pgenes = degenes[1:200]

bulk_exp1 = GetAssayData(obj_bulk_hESC, layer = 'data')[pgenes, obj_bulk_hESC$celltype == 'H1'] %>% as.matrix
bulk_exp2 = GetAssayData(obj_bulk_hESC, layer = 'data')[pgenes, obj_bulk_hESC$celltype == 'DEC'] %>% as.matrix
bulk_exp = cbind(bulk_exp1, bulk_exp2)

aggexp_bulk1 = rowMeans(bulk_exp1)
aggexp_bulk2 = rowMeans(bulk_exp2)

heat_bulk = Heatmap(bulk_exp, cluster_columns = F, show_row_dend = F, show_column_dend = F,
                    show_row_names = F, show_column_names = F,
                    col = colorRamp2(c(0,1,2,3,4), c('deepskyblue3', 'skyblue1', 'beige','orange', 'darkorange3')),
                    show_heatmap_legend = F)  
heat_bulk = draw(heat_bulk)   # Figure 2c, top right
gene_ord = row_order(heat_bulk)

#
cell_membership <- read.csv(paste0(data_path, 'mc1_cell_membership_rna_hESC.csv'), check.names = F, row.names = 1)

Gammas = as.character(colnames(cell_membership))

ResHeat_mc1 = list()

for (tgamma in Gammas) {
  
  if (tgamma == '1') {
    obj_metacell = obj_singlecell_hESC
  } else {
    sc_membership = cell_membership[[tgamma]]
    names(sc_membership) = rownames(cell_membership)
    sc_membership = sc_membership[colnames(obj_singlecell_hESC)]
    obj_metacell = mcRigor_buildmc(obj_singlecell_hESC, sc_membership, 
                                   aggregate_method = 'sum')
  }
  
  mc_exp1 = GetAssayData(obj_metacell, layer = 'data')[pgenes, obj_metacell$celltype == 'H1'] %>% as.matrix
  mc_exp2 = GetAssayData(obj_metacell, layer = 'data')[pgenes, obj_metacell$celltype == 'DEC'] %>% as.matrix
  mc_exp = cbind(mc_exp1, mc_exp2)
  
  aggexp_mc1 = rowMeans(mc_exp1)
  aggexp_mc2 = rowMeans(mc_exp2)
  
  heat_mc = Heatmap(mc_exp, row_order = gene_ord, cluster_columns = F, show_column_dend = F,
                    show_row_names = F, show_column_names = F,
                    col = colorRamp2(c(0,1,2,3,4), c('deepskyblue3', 'skyblue1', 'beige','orange', 'darkorange3')),
                    show_heatmap_legend = F)
  heat_mc   # Figure 2c, top right
  
  cor_cor = cor.test(cor(t(bulk_exp)), cor(t(mc_exp)), method = 'pearson')$estimate
  cor_exp1 = cor.test(aggexp_bulk1, aggexp_mc1, method = 'pearson')$estimate
  cor_exp2 = cor.test(aggexp_bulk2, aggexp_mc2, method = 'pearson')$estimate
  cor_exp = cor.test(c(aggexp_bulk1, aggexp_bulk2), c(aggexp_mc1, aggexp_mc2), method = 'pearson')$estimate
  
  ResHeat_mc1[[tgamma]] = list(cor_cor = cor_cor, cor_exp = cor_exp,
                               cor_exp1 = cor_exp1, cor_exp2 = cor_exp2)   # similar for seacells and supercell
  
}

save(ResHeat_seacells, ResHeat_supercell, ResHeat_mc1, file = paste0(data_path, 'ResHeat.RData'))   # intermediate results saved

#
load(file = paste0(data_path, 'ResHeat.RData'))

sc_dat = c(1, NA, ResHeat_seacells[['1']]$cor_exp)

COR = sapply(ResHeat_seacells, function(x) x$cor_exp)
pdat = data.frame(gamma = as.integer(names(ResHeat_seacells)),
                  cor = COR)
pdat = cut_bins(pdat, bin_width = 5)
pdat = rbind(pdat, sc_dat)

CORsup = sapply(ResHeat_supercell, function(x) x$cor_exp)
pdatsup = data.frame(gamma = as.integer(names(ResHeat_supercell)),
                     cor = CORsup)
pdatsup = cut_bins(pdatsup, bin_width = 5)
pdatsup = rbind(pdatsup, sc_dat)

CORmc1 = sapply(ResHeat_mc1, function(x) x$cor_exp)
pdatmc1 = data.frame(gamma = as.integer(names(ResHeat_mc1)),
                     cor = CORmc1)
pdatmc1 = cut_bins(pdatmc1, bin_width = 5)
pdatmc1 = pdatmc1[pdatmc1$gamma <=55,]
pdatmc1 = rbind(pdatmc1, sc_dat)

ppexp = ggplot(pdat, aes(x = gamma, y = cor)) + 
  geom_line(aes(col = 'SEACells'), linetype = 'longdash') + geom_point(aes(col = 'SEACells')) +
  annotate('point', x = 13, y = -Inf, size = 3, pch = 17, color = 'purple') +
  geom_line(data = pdatsup, aes(x = gamma, y = cor, col = 'SuperCell'), linetype = 'longdash') +
  geom_point(data = pdatsup, aes(x = gamma, y = cor, col = 'SuperCell')) +
  annotate('point', x = 4, y = -Inf, size = 3, pch = 17, color = 'darkolivegreen4') +
  geom_line(data = pdatmc1, aes(x = gamma, y = cor, col = 'MetaCell'), linetype = 'longdash') +
  geom_point(data = pdatmc1, aes(x = gamma, y = cor, col = 'MetaCell')) +
  annotate('point', x = 6, y = -Inf, size = 3, pch = 17, color = 'blue') +
  scale_color_manual(name="",values=c('SEACells'='purple', 'SuperCell'='darkolivegreen4', 'MetaCell'='blue')) +
  geom_vline(xintercept = 13, linetype = 'dashed', col = 'red') +
  coord_cartesian(ylim = c(0.76, 0.81), xlim = c(1,50), clip = 'off')
ppexp + theme(legend.position = 'none')   # Figure 2c, top left


