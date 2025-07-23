
library(Seurat)
library(mcRigor)

data_path <- './synthetic_data/'


###### Data loading (generated from syndata_generate.R)

obj_singlecell_syn = readRDS(file = paste0(data_path, 'syn.rds'))
obj_singlecell = obj_singlecell_syn


###### Figure 1b & Supplementary Figure 2 (MetaCell part)

# obtained by MetaCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "mc1_cell_membership_rna_syn.csv"), check.names = F, row.names = 1)

res_mc1 = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership, 
                         output_file = paste0(data_path, 'mc1_tabmc_RNA_syn.RData'))   # intermediate results saved


load(paste0(data_path, "mc1_tabmc_RNA_syn.RData"))

mc1_res = mcRigor_threshold(TabMC, pur_metric = 'metacell_purity')

library(aricode)

TabMC1 = mc1_res$TabMC

true_dub = TabMC1$metacell_purity
true_dub = ifelse(true_dub<1, 'dubious', 'trustworthy')
est_dub = TabMC1$mcRigor
ARI(true_dub, est_dub)
table(true_dub[est_dub == 'dubious'])
table(true_dub[est_dub == 'trustworthy']) 
recall = table(est_dub[true_dub == 'dubious'])['dubious'] / length(which(true_dub == 'dubious'))
precision = table(true_dub[est_dub == 'dubious'])['dubious'] / length(which(est_dub == 'dubious'))
FF = 2/(1/recall + 1/precision)
FF

mc1_res$test_plot   # Figure 1b bottom right & Supplementary Figure 2b 
mc1_res$purity_plot   # Figure 1b top right


###

tgamma = '61'

sc_membership = cell_membership[[tgamma]]
names(sc_membership) = rownames(cell_membership)

obj_singlecell = obj_singlecell_syn

sc.reduction = 'umap'
metric = 'size'

obj_singlecell[['Metacell']] = sc_membership
sc_membership = obj_singlecell$Metacell

if(is.character(sc.reduction)){
  # if sc_reduction does not exist compute pca and run UMAP:
  if(is.null(obj_singlecell@reductions[[sc.reduction]])){
    message("Low dimensionnal embessing not found in obj_singlecell")
    message("Computing PCA ...")
    obj_singlecell <- Seurat::NormalizeData(obj_singlecell)
    obj_singlecell <- Seurat::FindVariableFeatures(obj_singlecell)
    obj_singlecell <- Seurat::ScaleData(obj_singlecell)
    obj_singlecell <- Seurat::RunPCA(obj_singlecell, verbose = F)
    message("Running UMAP ...")
    obj_singlecell <- Seurat::RunUMAP(obj_singlecell, reduction = "pca", dims = c(1:20), verbose = F)
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[["umap"]])
  } else{
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])
  } 
} else{
  scCoord <- sc.reduction
}  

centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot

obj_metacell = mcRigor_buildmc(obj_singlecell, sc_membership, doNorm = F,
                               add_testres = F)

centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]

centroids[[metric]] <- obj_metacell[[metric]][,1]

comp_metric1 = 'metacell_purity'
comp_metric2 = 'TT_div'
centroids[[comp_metric1]] = obj_metacell[[comp_metric1]][,1]
centroids[[comp_metric2]] = TabMC[centroids$sc_membership, ][[comp_metric2]]

cor.test(centroids[[comp_metric1]], centroids[[comp_metric2]], method = 'pearson')

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             aspect.ratio = 2.3/3)  

p3 = ggplot(data = centroids, aes_string(x = comp_metric1, y = comp_metric2)) +
  geom_point()

p3 + theme(legend.position = 'none')   # Supplementary Figure 2a


p1 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric1, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'purity', # limits = c(0.3, 1),
                        option = 'C') +
  guides(size = 'none')
p1

p2 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric2, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'mcDiv', # limits = c(1, 2),
                        option = 'C', direction = -1) +
  guides(size = 'none')
p2

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             axis.ticks = element_blank(), axis.text = element_blank(),
             aspect.ratio = 2.3/3)  #

p1+p2   # Figure 1b left




###### Supplementary Figure 1d-f & Supplementary Figure (other)

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "seacells_cell_membership_rna_syn.csv"), check.names = F, row.names = 1)

res_seacells = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership, 
                         output_file = paste0(data_path, 'seacells_tabmc_RNA_syn.RData'))   # intermediate results saved

load(paste0(data_path, "seacells_tabmc_RNA_syn.RData"))

seacells_res = mcRigor_threshold(TabMC, pur_metric = 'metacell_purity')

library(aricode)

TabMC1 = seacells_res$TabMC

true_dub = TabMC1$metacell_purity
true_dub = ifelse(true_dub<1, 'dubious', 'trustworthy')
est_dub = TabMC1$mcRigor
ARI(true_dub, est_dub)
table(true_dub[est_dub == 'dubious'])
table(true_dub[est_dub == 'trustworthy']) 
recall = table(est_dub[true_dub == 'dubious'])['dubious'] / length(which(true_dub == 'dubious'))
precision = table(true_dub[est_dub == 'dubious'])['dubious'] / length(which(est_dub == 'dubious'))
FF = 2/(1/recall + 1/precision)
FF

seacells_res$test_plot   # Supplementary Figure 1d bottom right & Supplementary Figure 2b 
seacells_res$purity_plot   # Supplementary Figure 1d top right


###

tgamma = '60'

sc_membership = cell_membership[[tgamma]]
names(sc_membership) = rownames(cell_membership)

sc.reduction = 'umap'
metric = 'size'

obj_singlecell[['Metacell']] = sc_membership
sc_membership = obj_singlecell$Metacell

if(is.character(sc.reduction)){
  # if sc_reduction does not exist compute pca and run UMAP:
  if(is.null(obj_singlecell@reductions[[sc.reduction]])){
    message("Low dimensionnal embessing not found in obj_singlecell")
    message("Computing PCA ...")
    obj_singlecell <- Seurat::NormalizeData(obj_singlecell)
    obj_singlecell <- Seurat::FindVariableFeatures(obj_singlecell)
    obj_singlecell <- Seurat::ScaleData(obj_singlecell)
    obj_singlecell <- Seurat::RunPCA(obj_singlecell, verbose = F)
    message("Running UMAP ...")
    obj_singlecell <- Seurat::RunUMAP(obj_singlecell, reduction = "pca", dims = c(1:20), verbose = F)
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[["umap"]])
  } else{
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])
  } 
} else{
  scCoord <- sc.reduction
}  

centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot

obj_metacell = mcRigor_buildmc(obj_singlecell, sc_membership, doNorm = F,
                               add_testres = F)

centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]

centroids[[metric]] <- obj_metacell[[metric]][,1]

comp_metric1 = 'metacell_purity'
comp_metric2 = 'TT_div'
centroids[[comp_metric1]] = obj_metacell[[comp_metric1]][,1]
centroids[[comp_metric2]] = TabMC[centroids$sc_membership, ][[comp_metric2]]

cor.test(centroids[[comp_metric1]], centroids[[comp_metric2]], method = 'pearson')

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             aspect.ratio = 2.3/3)  

p3 = ggplot(data = centroids, aes_string(x = comp_metric1, y = comp_metric2)) +
  geom_point()

p3 + theme(legend.position = 'none')   # Supplementary Figure 2a


p1 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric1, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'purity', # limits = c(0.3, 1),
                        option = 'C') +
  guides(size = 'none')
p1

p2 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric2, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'mcDiv', # limits = c(1, 2),
                        option = 'C', direction = -1) +
  guides(size = 'none')
p2

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             axis.ticks = element_blank(), axis.text = element_blank(),
             aspect.ratio = 2.3/3)  #

p1+p2   # Supplementary Figure 1d left


###

# obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "supercell_cell_membership_rna_syn.csv"), check.names = F, row.names = 1)

res_supercell = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership, 
                              output_file = paste0(data_path, 'supercell_tabmc_RNA_syn.RData'))   # intermediate results saved

load(paste0(data_path, "supercell_tabmc_RNA_syn.RData"))

supercell_res = mcRigor_threshold(TabMC, pur_metric = 'metacell_purity', test_cutoff = 0.05)

library(aricode)

TabMC1 = supercell_res$TabMC

true_dub = TabMC1$metacell_purity
true_dub = ifelse(true_dub<1, 'dubious', 'trustworthy')
est_dub = TabMC1$mcRigor
ARI(true_dub, est_dub)
table(true_dub[est_dub == 'dubious'])
table(true_dub[est_dub == 'trustworthy']) 
recall = table(est_dub[true_dub == 'dubious'])['dubious'] / length(which(true_dub == 'dubious'))
precision = table(true_dub[est_dub == 'dubious'])['dubious'] / length(which(est_dub == 'dubious'))
FF = 2/(1/recall + 1/precision)
FF

supercell_res$test_plot   # Supplementary Figure 1e bottom right & Supplementary Figure 2b 
supercell_res$purity_plot   # Supplementary Figure 1e top right


###

tgamma = '60'

sc_membership = cell_membership[[tgamma]]
names(sc_membership) = rownames(cell_membership)

sc.reduction = 'umap'
metric = 'size'

obj_singlecell[['Metacell']] = sc_membership
sc_membership = obj_singlecell$Metacell

if(is.character(sc.reduction)){
  # if sc_reduction does not exist compute pca and run UMAP:
  if(is.null(obj_singlecell@reductions[[sc.reduction]])){
    message("Low dimensionnal embessing not found in obj_singlecell")
    message("Computing PCA ...")
    obj_singlecell <- Seurat::NormalizeData(obj_singlecell)
    obj_singlecell <- Seurat::FindVariableFeatures(obj_singlecell)
    obj_singlecell <- Seurat::ScaleData(obj_singlecell)
    obj_singlecell <- Seurat::RunPCA(obj_singlecell, verbose = F)
    message("Running UMAP ...")
    obj_singlecell <- Seurat::RunUMAP(obj_singlecell, reduction = "pca", dims = c(1:20), verbose = F)
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[["umap"]])
  } else{
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])
  } 
} else{
  scCoord <- sc.reduction
}  

centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot

obj_metacell = mcRigor_buildmc(obj_singlecell, sc_membership, doNorm = F,
                               add_testres = F)

centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]

centroids[[metric]] <- obj_metacell[[metric]][,1]

comp_metric1 = 'metacell_purity'
comp_metric2 = 'TT_div'
centroids[[comp_metric1]] = obj_metacell[[comp_metric1]][,1]
centroids[[comp_metric2]] = TabMC[centroids$sc_membership, ][[comp_metric2]]

cor.test(centroids[[comp_metric1]], centroids[[comp_metric2]], method = 'pearson')

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             aspect.ratio = 2.3/3)  

p3 = ggplot(data = centroids, aes_string(x = comp_metric1, y = comp_metric2)) +
  geom_point()

p3 + theme(legend.position = 'none')   # Supplementary Figure 2a


p1 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric1, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'purity', # limits = c(0.3, 1),
                        option = 'C') +
  guides(size = 'none')
p1

p2 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric2, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'mcDiv', # limits = c(1, 2),
                        option = 'C', direction = -1) +
  guides(size = 'none')
p2

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             axis.ticks = element_blank(), axis.text = element_blank(),
             aspect.ratio = 2.3/3)  #

p1+p2   # Supplementary Figure 1e left


###

# obtained by MetaQ (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "metaq_cell_membership_rna_syn.csv"), check.names = F, row.names = 1)

res_mc2 = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership, 
                              output_file = paste0(data_path, 'metaq_tabmc_RNA_syn.RData'))   # intermediate results saved

load(paste0(data_path, "metaq_tabmc_RNA_syn.RData"))

metaq_res = mcRigor_threshold(TabMC, pur_metric = 'metacell_purity')

library(aricode)

TabMC1 = metaq_res$TabMC

true_dub = TabMC1$metacell_purity
true_dub = ifelse(true_dub<1, 'dubious', 'trustworthy')
est_dub = TabMC1$mcRigor
ARI(true_dub, est_dub)
table(true_dub[est_dub == 'dubious'])
table(true_dub[est_dub == 'trustworthy']) 
recall = table(est_dub[true_dub == 'dubious'])['dubious'] / length(which(true_dub == 'dubious'))
precision = table(true_dub[est_dub == 'dubious'])['dubious'] / length(which(est_dub == 'dubious'))

FF = 2/(1/recall + 1/precision)
FF

metaq_res$test_plot   # Supplementary Figure 1f bottom right & Supplementary Figure 2b 
metaq_res$purity_plot   # Supplementary Figure 1f top right


###

tgamma = '60'

sc_membership = cell_membership[[tgamma]]
names(sc_membership) = rownames(cell_membership)

sc.reduction = 'umap'
metric = 'size'

obj_singlecell[['Metacell']] = sc_membership
sc_membership = obj_singlecell$Metacell

if(is.character(sc.reduction)){
  # if sc_reduction does not exist compute pca and run UMAP:
  if(is.null(obj_singlecell@reductions[[sc.reduction]])){
    message("Low dimensionnal embessing not found in obj_singlecell")
    message("Computing PCA ...")
    obj_singlecell <- Seurat::NormalizeData(obj_singlecell)
    obj_singlecell <- Seurat::FindVariableFeatures(obj_singlecell)
    obj_singlecell <- Seurat::ScaleData(obj_singlecell)
    obj_singlecell <- Seurat::RunPCA(obj_singlecell, verbose = F)
    message("Running UMAP ...")
    obj_singlecell <- Seurat::RunUMAP(obj_singlecell, reduction = "pca", dims = c(1:20), verbose = F)
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[["umap"]])
  } else{
    scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])
  } 
} else{
  scCoord <- sc.reduction
}  

centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot

obj_metacell = mcRigor_buildmc(obj_singlecell, sc_membership, doNorm = F,
                               add_testres = F)

centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]

centroids[[metric]] <- obj_metacell[[metric]][,1]

comp_metric1 = 'metacell_purity'
comp_metric2 = 'TT_div'
centroids[[comp_metric1]] = obj_metacell[[comp_metric1]][,1]
centroids[[comp_metric2]] = TabMC[centroids$sc_membership, ][[comp_metric2]]

cor.test(centroids[[comp_metric1]], centroids[[comp_metric2]], method = 'pearson')

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             aspect.ratio = 2.3/3)  

p3 = ggplot(data = centroids, aes_string(x = comp_metric1, y = comp_metric2)) +
  geom_point()

p3 + theme(legend.position = 'none')   # Supplementary Figure 2a


p1 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric1, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'purity', # limits = c(0.3, 1),
                        option = 'C') +
  guides(size = 'none')
p1

p2 <- ggplot2::ggplot(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2]))+ 
  ggplot2::geom_point(mapping = ggplot2::aes_string(color = comp_metric2, 
                                                    size = metric),
                      pch= 16, alpha = 0.5) +
  scale_size_continuous(limits = c(1,1000), range = c(1,12)) +
  scale_color_viridis_c(name = 'mcDiv', # limits = c(1, 2),
                        option = 'C', direction = -1) +
  guides(size = 'none')
p2

theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             strip.background = element_rect(fill="white"),
             panel.grid = element_blank(),
             axis.ticks = element_blank(), axis.text = element_blank(),
             aspect.ratio = 2.3/3)  #

p1+p2   # Supplementary Figure 1f left


