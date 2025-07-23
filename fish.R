
library(Seurat)
library(SeuratData)
library(mcRigor)

data_path <- './fish_data/'


###### Data loading & preprocessing

dropseq <- readRDS(paste0(data_path, "melanoma_dropseq.rds"))  # gene by cell

obj_singlecell_fish = CreateSeuratObject(counts = dropseq, 
                                         assay = 'RNA',
                                         min.cells = 0,
                                         min.features = 0)

obj_singlecell_fish <- NormalizeData(obj_singlecell_fish)
obj_singlecell_fish <- FindVariableFeatures(obj_singlecell_fish, nfeatures = 2000)
obj_singlecell_fish <- ScaleData(obj_singlecell_fish)
obj_singlecell_fish = RunPCA(obj_singlecell_fish)
ElbowPlot(obj_singlecell_fish)
obj_singlecell_fish = RunUMAP(obj_singlecell_fish, dims = 1:15)
UMAPPlot(obj_singlecell_fish)

obj_singlecell_fish = FindNeighbors(obj_singlecell_fish, dims = 1:15)
obj_singlecell_fish = FindClusters(obj_singlecell_fish, resolution = 0.5)

obj_singlecell_fish[['louvain_type']] = obj_singlecell_fish[['seurat_clusters']]

saveRDS(obj_singlecell_fish, file = paste0(data_path, 'fish.rds'))

#
obj_singlecell_fish = readRDS(file = paste0(data_path, 'fish.rds'))


### smRNA FISH data

fish <- read.table(paste0(data_path, "fishSubset.txt"), header = TRUE, row.names = 1)

fish.filt <- fish[fish[, "GAPDH"] < quantile(fish[, "GAPDH"], 0.9) & 
                    fish[, "GAPDH"] > quantile(fish[, "GAPDH"], 0.1), ]
fish.norm <- sweep(fish.filt, 1, fish.filt[, "GAPDH"]/mean(fish.filt[, "GAPDH"]), "/")

n.genes <- nrow(dropseq)
n.cells <- ncol(dropseq)
gene.names <- rownames(dropseq)
cell.names <- colnames(dropseq)

genes <- which(gene.names %in% colnames(fish.norm))

sf <- colSums(dropseq)/mean(colSums(dropseq))

normalize.gapdh <- function(x) {
  x.filt <- x[, x["GAPDH", ] < quantile(x["GAPDH", ], 0.9) &
                x["GAPDH", ] > quantile(x["GAPDH", ], 0.1)]
  x.norm <- sweep(x.filt, 2, x.filt["GAPDH", ]/mean(x.filt["GAPDH", ]), "/")
}

dropseq.filt <- dropseq[genes, dropseq["GAPDH", ] < quantile(dropseq["GAPDH", ], 0.9) &
                          dropseq["GAPDH", ] > quantile(dropseq["GAPDH", ], 0.1)]
dropseq.norm <- normalize.gapdh(dropseq[genes, ])



###### Run mcRigor

# obtained by MetaCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "mc1_cell_membership_rna_fish.csv"), check.names = F, row.names = 1)

res = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership, weight = NULL,
                     output_file = paste0(data_path, 'mc1_tabmc_RNA_fish.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_seacells = read.csv(paste0(data_path,
                                  "seacells_cell_membership_rna_fish.csv"), check.names = F, row.names = 1)

res_seacells = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_seacells, weight = NULL,
                     output_file = paste0(data_path, 'seacells_tabmc_RNA_fish.RData'))   # intermediate results saved

# obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_supercell = read.csv(paste0(data_path,
                                  "supercell_cell_membership_rna_fish.csv"), check.names = F, row.names = 1)

res_supercell = mcRigor_OPTIMIZE(obj_singlecell = obj_singlecell, cell_membership = cell_membership_supercell, weight = NULL, 
                     output_file = paste0(data_path, 'supercell_tabmc_RNA_fish.RData'))   # intermediate results saved




######  Figure 2b & Supplementary Figure 8

load(paste0(data_path, "mc1_tabmc_RNA_fish.RData"))

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
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 12 / (max(DD$gamma) - min(DD$gamma)))$y
}
DD$score = DD$ZeroRate * 0.3 + DD$D * 0.7
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
pelbow   # Supplementary Figure 8a

sc_membership = cell_membership[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership)

mcRigor_projection(obj_singlecell = obj_singlecell_fish,
                   sc_membership = sc_membership,
                   color_field = 'orig.ident',
                   cpalette =  c('#86ABD4'),
                   add_testres = T, test_stats = TabMC,
                   dub_mc_test.label = T)   # Supplementary Figure 8b

obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_fish, sc_membership = sc_membership,
                               add_testres = T, test_stats = TabMC, Thre = Thre)


#

load(paste0(data_path, "seacells_tabmc_RNA_fish.RData"))

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
DD$score = DD$ZeroRate * 0.3 + DD$D * 0.7
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
pelbow   # Supplementary Figure 8a

sc_membership = cell_membership_seacells[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership_seacells)

mcRigor_projection(obj_singlecell = obj_singlecell_fish,
                   sc_membership = sc_membership,
                   color_field = 'orig.ident',
                   cpalette =  c('#86ABD4'),
                   add_testres = T, test_stats = TabMC,
                   dub_mc_test.label = T)   # Supplementary Figure 8b

obj_metacell_sea = mcRigor_buildmc(obj_singlecell = obj_singlecell_fish, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)

#
load(paste0(data_path, "supercell_tabmc_RNA_fish.RData"))

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
  DD$D = lowess(x = DD$gamma, y = DD$D, f = 16 / (max(DD$gamma) - min(DD$gamma)))$y
}
DD$score = DD$ZeroRate * 0.7 + DD$D * 0.3
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
pelbow   # Supplementary Figure 8a

sc_membership = cell_membership_supercell[[as.character(opt_gamma)]]
names(sc_membership) = rownames(cell_membership_supercell)

mcRigor_projection(obj_singlecell = obj_singlecell_fish,
                   sc_membership = sc_membership,
                   color_field = 'orig.ident',
                   cpalette =  c('#86ABD4'),
                   add_testres = T, test_stats = TabMC,
                   dub_mc_test.label = T)   # Supplementary Figure 8b

obj_metacell_sup = mcRigor_buildmc(obj_singlecell = obj_singlecell_fish, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thre)



##

normalize.gapdh <- function(x) {
  x.filt <- x[, x["GAPDH", ] < quantile(x["GAPDH", ], 0.9) &
                x["GAPDH", ] > quantile(x["GAPDH", ], 0.1)]
  x.norm <- sweep(x.filt, 2, x.filt["GAPDH", ]/mean(x.filt["GAPDH", ]), "/")
}

genes <- intersect(colnames(fish), rownames(obj_singlecell_fish))

fish.norm = fish.norm[,genes]
sc.count = GetAssayData(obj_singlecell_fish, layer = 'counts')[genes,]
sc.norm = normalize.gapdh(sc.count) %>% as.matrix %>% t
mc.count = GetAssayData(obj_metacell, layer = 'counts')[genes,]
mc.norm = normalize.gapdh(mc.count) %>% as.matrix %>% t
sea.count = GetAssayData(obj_metacell_sea, layer = 'counts')[genes,]
sea.norm = normalize.gapdh(sea.count) %>% as.matrix %>% t
sup.count = GetAssayData(obj_metacell_sup, layer = 'counts')[genes,]
sup.norm = normalize.gapdh(sup.count) %>% as.matrix %>% t

mult.factor_sc <- sapply(genes, function(x)
  mean(fish.norm[, x], na.rm = TRUE)/mean(sc.norm[, x]))
mult.factor_mc <- sapply(genes, function(x)
  mean(fish.norm[, x], na.rm = TRUE)/mean(mc.norm[, x]))
mult.factor_sea <- sapply(genes, function(x)
  mean(fish.norm[, x], na.rm = TRUE)/mean(sea.norm[, x]))
mult.factor_sup <- sapply(genes, function(x)
  mean(fish.norm[, x], na.rm = TRUE)/mean(sup.norm[, x]))


##
for (method in c('mc1', 'seacells', 'supercell')){
  
  cell_membership = read.csv(paste0(data_path, method, "_cell_membership_rna_fish.csv"), check.names = F, row.names = 1)
  load(paste(data_path, method, "_tabmc_RNA_fish.RData"))
  
  rigorMC_res = mcRigor_threshold(TabMC, test_cutoff = 0.05)
  Thret = rigorMC_res$threshold
  
  genes <- intersect(colnames(fish), rownames(obj_singlecell_fish))
  
  sc.count = GetAssayData(obj_singlecell_fish, layer = 'counts')[genes,]
  sc.count = t(as.matrix(sc.count))
  
  zeropdat = data.frame(fish = apply(fish[,genes], 2, function(x) length(which(abs(x)<=0.000001)) / length(which(!is.na(x)))),
                        sc = apply(sc.count, 2, function(x) length(which(abs(x)<=0.000001)) / length(which(!is.na(x)))))
  
  for (tgamma in colnames(cell_membership)){
    
    sc_membership = cell_membership[[tgamma]]
    names(sc_membership) = rownames(cell_membership)
    
    obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_fish, sc_membership = sc_membership,
                                   add_testres = T, test_stats = TabMC, Thre = Thret)
    
    mc.count = GetAssayData(obj_metacell, layer = 'counts')[genes,]
    mc.count = t(as.matrix(mc.count))
    
    binary_mc_sc = sc.count
    for (j in 1:dim(binary_mc_sc)[1]) {
      scid = rownames(binary_mc_sc)[j]
      mcid = sc_membership[scid]
      if (!is.na(mcid) && mcid != '') binary_mc_sc[j,] = mc.count[mcid,]
    }
    sc_testres = obj_metacell@misc$cell_membership$testres_sc
    sc_testres = sc_testres[match(rownames(binary_mc_sc), 
                                  rownames(obj_metacell@misc$cell_membership))]
    
    zeropdat[[tgamma]] = apply(binary_mc_sc, 2, function(x) length(which(abs(x)<=0.000001)) / length(which(!is.na(x))))
  }
  
  save(zeropdat, file = paste0(data_path, 'zero_', method, '.RData'))   # intermediate results saved
  
}


##
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
load(paste0(data_path, "zero_seacells.RData"))
zero_sea = zeropdat
load(paste0(data_path, "zero_mc1.RData"))
zero_mc1 = zeropdat
load(paste0(data_path, "zero_supercell.RData"))
zero_sup = zeropdat

ppdat = data.frame(gamma = c(7:100))
ppdat$sea = colMeans(zero_sea[, as.character(ppdat$gamma)]) 
ppdat$mc1 = colMeans(zero_mc1[, as.character(ppdat$gamma)]) 
ppdat$sup = colMeans(zero_sup[, as.character(ppdat$gamma)]) 

ppdat = cut_bins(ppdat)

pp = ggplot(ppdat, aes(x = gamma)) + 
  geom_path(aes(y = sea, color = 'SEACells'), linetype = 'longdash') + 
  geom_point(aes(y = sea, col = 'SEACells')) + 
  annotate('point', x = 44, y = -Inf, size = 3, pch = 17, color = 'purple') +   
  geom_path(aes(y = mc1, color = 'MetaCell'), linetype = 'longdash') + 
  geom_point(aes(y = mc1, color = 'MetaCell')) +
  annotate('point', x = 24, y = -Inf, size = 3, pch = 17, color = 'blue') +   
  geom_path(aes(y = sup, color = 'SuperCell'), linetype = 'longdash') + 
  geom_point(aes(y = sup, color = 'SuperCell')) +
  annotate('point', x = 56, y = -Inf, size = 3, pch = 17, color = 'darkolivegreen4') +   
  # geom_vline(xintercept = seacells_res$optimized$gamma, color = 'purple') + 
  scale_color_manual(name="",values=c('SEACells'='purple', 'SuperCell'='darkolivegreen4', 'MetaCell'='blue')) +
  geom_abline(slope = 0, intercept = mean(zeropdat[, 'fish']), color = 'red') +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_cartesian(clip = 'off')
pp   # Figure 2b, top


##
load(paste0(data_path, "zero_seacells.RData"))
zero_sea = zeropdat
load(paste0(data_path, "zero_mc1.RData"))
zero_mc1 = zeropdat
load(paste0(data_path, "zero_supercell.RData"))
zero_sup = zeropdat

gene = 'FGFR1'  # 'RUNX2' 'KDM5B'

ppdat = data.frame(gamma = c(7:100))
ppdat$sea = zero_sea[gene, as.character(ppdat$gamma)] %>% unlist
ppdat$mc1 = zero_mc1[gene, as.character(ppdat$gamma)] %>% unlist
ppdat$sup = zero_sup[gene, as.character(ppdat$gamma)] %>% unlist

ppdat = cut_bins(ppdat)

pp = ggplot(ppdat, aes(x = gamma)) + 
  geom_path(aes(y = sea, color = 'SEACells'), linetype = 'longdash') + 
  geom_point(aes(y = sea, col = 'SEACells')) + 
  annotate('point', x = 44, y = -Inf, size = 3, pch = 17, color = 'purple') +
  geom_path(aes(y = mc1, color = 'MetaCell'), linetype = 'longdash') + 
  geom_point(aes(y = mc1, color = 'MetaCell')) +
  annotate('point', x = 24, y = -Inf, size = 3, pch = 17, color = 'blue') +
  geom_path(aes(y = sup, color = 'SuperCell'), linetype = 'longdash') + 
  geom_point(aes(y = sup, color = 'SuperCell')) +
  annotate('point', x = 56, y = -Inf, size = 3, pch = 17, color = 'darkolivegreen4') +
  # geom_vline(xintercept = seacells_res$optimized$gamma, color = 'purple') + 
  scale_color_manual(name="",values=c('SEACells'='purple', 'SuperCell'='darkolivegreen4', 'MetaCell'='blue')) +
  geom_abline(slope = 0, intercept = zeropdat[gene, 'fish'], color = 'red') +
  ggtitle(gene) +
  theme(plot.title = element_text(hjust = 0.5), axis.title.x = element_blank(), axis.title.y = element_blank()) +
  coord_cartesian(clip = 'off')
pp   # Figure 2b, bottom & Supplementary Figure 8c



## Supplementary Figure 8d

gene = 'LMNA'

par(mfrow = c(1, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.04), 
     xlim = c(-10, 500), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 500, by = 100))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext(gene, side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, gene], na.rm = TRUE)$bw
lines(density(fish.norm[, gene], na.rm = TRUE), lwd = 4, lty = 1, col = "red")
lines(density(sc.norm[, gene]*mult.factor_sc[gene], bw = dens.bw), lwd = 2, col = 'black')
lines(density(mc.norm[, gene]*mult.factor_mc[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'blue')
lines(density(sea.norm[, gene]*mult.factor_sea[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'purple')
lines(density(sup.norm[, gene]*mult.factor_sup[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'darkolivegreen4')

#
gene = 'CCNA2'

par(mfrow = c(1, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.25), 
     xlim = c(-2, 60), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 60, by = 20))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext(gene, side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, gene], na.rm = TRUE)$bw
lines(density(fish.norm[, gene], na.rm = TRUE), lwd = 4, lty = 1, col = "red")
lines(density(sc.norm[, gene]*mult.factor_sc[gene], bw = dens.bw), lwd = 2, col = 'black')
lines(density(mc.norm[, gene]*mult.factor_mc[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'blue')
lines(density(sea.norm[, gene]*mult.factor_sea[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'purple')
lines(density(sup.norm[, gene]*mult.factor_sup[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'darkolivegreen4')

#
gene = 'FGFR1'

par(mfrow = c(1, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.36), 
     xlim = c(-3, 30), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 30, by = 10))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext(gene, side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, gene], na.rm = TRUE)$bw
lines(density(fish.norm[, gene], na.rm = TRUE), lwd = 4, lty = 1, col = "red")
lines(density(sc.norm[, gene]*mult.factor_sc[gene], bw = dens.bw), lwd = 2, col = 'black')
lines(density(mc.norm[, gene]*mult.factor_mc[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'blue')
lines(density(sea.norm[, gene]*mult.factor_sea[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'purple')
lines(density(sup.norm[, gene]*mult.factor_sup[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'darkolivegreen4')

#
gene = 'RUNX2'

par(mfrow = c(1, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 2.2), 
     xlim = c(-0.3, 6), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 6, by = 2))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext(gene, side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, gene], na.rm = TRUE)$bw
lines(density(fish.norm[, gene], na.rm = TRUE), lwd = 4, lty = 1, col = "red")
lines(density(sc.norm[, gene]*mult.factor_sc[gene], bw = dens.bw), lwd = 2, col = 'black')
lines(density(mc.norm[, gene]*mult.factor_mc[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'blue')
lines(density(sea.norm[, gene]*mult.factor_sea[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'purple')
lines(density(sup.norm[, gene]*mult.factor_sup[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'darkolivegreen4')

#
gene = 'KDM5B'

par(mfrow = c(1, 1), cex.main = 1.5, mar = c(1, 2, 1.5, 0) + 0.1, oma = c(3, 2, 0, 2), 
    mgp = c(3.5, 1, 0),
    cex.lab = 1.5, font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)
plot(0, type = "n", ylab = "", xlab = "", cex = 1.5, ylim = c(0, 0.8), 
     xlim = c(-1, 20), lwd = 2, pch = 5, axes = FALSE, main = "")
axis(1, seq(0, 20, by = 5))
axis(2, labels = FALSE, lwd.ticks = 0)
par(las = 0)
mtext("Density", side = 2, line = 1, cex = 1.5)
mtext(gene, side = 3, line = -2, cex = 1.5, font = 3)
dens.bw <- density(fish.norm[, gene], na.rm = TRUE)$bw
lines(density(fish.norm[, gene], na.rm = TRUE), lwd = 4, lty = 1, col = "red")
lines(density(sc.norm[, gene]*mult.factor_sc[gene], bw = dens.bw), lwd = 2, col = 'black')
lines(density(mc.norm[, gene]*mult.factor_mc[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'blue')
lines(density(sea.norm[, gene]*mult.factor_sea[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'purple')
lines(density(sup.norm[, gene]*mult.factor_sup[gene], bw = dens.bw), lwd = 3, lty = 5, col = 'darkolivegreen4')


