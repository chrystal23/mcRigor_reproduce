
library(Seurat)
library(SeuratData)
library(mcRigor)

data_path <- './bmcite_data/'


###### Data loading

InstallData('bmcite')

data("bmcite")
bmcite
head(bmcite@meta.data)
bmcite$celltype_simplified <- plyr::revalue(bmcite$celltype.l2, 
                                            c("CD8 Effector_1" = "Non-Naive CD8 cell",
                                              "CD8 Effector_2" = "Non-Naive CD8 cell",
                                              "CD8 Memory_1" = "Non-Naive CD8 cell",
                                              "CD8 Memory_2" = "Non-Naive CD8 cell",
                                              "CD8 Naive" = "Naive CD8 cell",
                                              "CD4 Naive" = "Naive CD4 cell",
                                              "CD4 Memory" = "Non-Naive CD4 cell",
                                              "Treg" = "Non-Naive CD4 cell",
                                              "Naive B" = "B cell",
                                              "Memory B" = "B cell",
                                              "CD56 bright NK" = "NK",
                                              "MAIT" = "Unconventional T",
                                              "gdT" = "Unconventional T"
                                            ))
bmcite <- bmcite[,-grep("Prog",bmcite$celltype_simplified)]
if(packageVersion("Seurat") >= 5) {
  bmcite[["RNA"]] <- as(object = bmcite[["RNA"]], Class = "Assay")
}

obj_singlecell <- bmcite[, bmcite$donor=='batch1']


###### Quality control & Preprocessing

obj_singlecell[['percent.mt']] <- PercentageFeatureSet(obj_singlecell, pattern = "^MT-")
View(head(obj_singlecell@meta.data, 10))

FeatureScatter(obj_singlecell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')
VlnPlot(obj_singlecell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
min(obj_singlecell$nFeature_RNA)

obj_singlecell <- subset(obj_singlecell, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & 
                           percent.mt < 10)

obj_singlecell <- NormalizeData(obj_singlecell)
obj_singlecell <- FindVariableFeatures(obj_singlecell, nfeatures = 2000)
obj_singlecell <- ScaleData(obj_singlecell, features = rownames(obj_singlecell))

celltype_colors <- c("#7E57C2", "#1E88E5", "#FFC107", "#004D40", "#9E9D24", 
                     "#F06292", "#546E7A", "#D4E157", "#76FF03", "#6D4C41",
                     "#26A69A", "#AB47BC", "#EC407A", "#D81B60", "#42A5F5",
                     "#2E7D32", "#FFA726", "#5E35B1", "#EF5350", "#3949AB",
                     'darkblue', 'yellow')
names(celltype_colors) <- unique(obj_singlecell$celltype.l1)
names(celltype_colors) <- unique(obj_singlecell$celltype.l2)

obj_singlecell <- RunPCA(obj_singlecell)
DimPlot(obj_singlecell, dims = c(1,2), reduction = 'pca', 
        group.by = 'celltype.l1', cols = celltype_colors) + theme_classic() + NoLegend()


##### Run scDEED to obtain the optimal visualization

# devtools::install_github("JSB-UCLA/scDEED")
library('scDEED')

chooseK(obj_singlecell)

tsne_sc <- scDEED(obj_singlecell, num_pc = 19, use_method = 'tsne', visualization = T)

tsne_sc$num_dubious$perplexity[tsne_sc$num_dubious$number_dubious_cells==min(tsne_sc$num_dubious$number_dubious_cells)]
# optimized perplexity is 750

obj_singlecell =  RunTSNE(obj_singlecell, dims = 1:19, perplexity = 750)

# save preprocessed single cell data
saveRDS(obj_singlecell, file = paste0(data_path, 'bmcite_mcRig.rds'))





###### Run mcRigor

obj_singlecell <- readRDS(paste0(data_path, 'bmcite_mcRig.rds'))

TSNEPlot(obj_singlecell, group.by = 'celltype_simplified') + theme(plot.title = element_blank())

# obtained by MetaCell2 in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path,
                                  "cell_membership_bmcite.csv"), check.names = F, row.names = 1)

res = mcRigor_DETECT(obj_singlecell = obj_singlecell, cell_membership = cell_membership, 
                     output_file = paste0(data_path, 'Tabmc_bmcite.RData'))   # intermediate results saved

# use the metacell partitioning under gamma = 30 as an example
sc_membership = cell_membership[['30']]
names(sc_membership) = rownames(cell_membership)

obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell, sc_membership = sc_membership,
                               add_testres = T, test_stats = res$TabMC, Thre = res$thre)
# OR (since the intermediate results has been saved as 'Tabmc_bmcite.RData')
load(paste0(data_path, 'Tabmc_bmcite.RData'))
obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell, sc_membership = sc_membership,
                               add_testres = T, test_stats = TabMC)

saveRDS(obj_metacell, file = paste0(data_path, 'obj_mc_bmcite.rds'))   # intermediate results saved

obj_metacell = readRDS(obj_metacell, file = paste0(data_path, 'obj_mc_bmcite.rds'))






###### Figure 1c

sc.reduction = 'tsne'
metric = 'size'
max_mcsize = 1000
dub_mc.label = T
pt_size = 1
sc.alpha = 0.5
mc.alpha = 1
color_field = 'celltype_simplified'

cpalette = c('#13678A', '#45C4B0', '#9AEBA3', '#BF847E', '#86ABD4', '#F2C12E', '#FACFCE', '#DAFDBA',
             "#7E57C2", "#42A5F5", "#9E9D24", "#546E7A", "#D4E157", "#76FF03", "#6D4C41", "#004D40",
             "#AB47BC", "#D81B60")

scCoord <- Seurat::Embeddings(obj_singlecell@reductions[[sc.reduction]])

sc_membership <- obj_metacell@misc$cell_membership$metacell_name

centroids <- stats::aggregate(scCoord~sc_membership, scCoord, mean) #should be taken from object slot

centroids = centroids[match(colnames(obj_metacell), centroids[,1]),]

centroids[[metric]] <- obj_metacell[[metric]][,1]

obj_metacell$testres = obj_metacell$mcRigor

centroids[['dub_mc_test']] = NA
centroids[['dub_mc_test']][obj_metacell$testres == 'dubious'] = colnames(obj_metacell)[obj_metacell$testres == 'dubious']

centroids[[color_field]] <- obj_metacell[[color_field]][,1]

scCoord <- data.frame(scCoord)
scCoord[[color_field]] <- obj_singlecell[[color_field]][,1]

fig <- ggplot2::ggplot(scCoord,
                     ggplot2::aes_string(colnames(scCoord)[1],
                                         colnames(scCoord)[2],
                                         color = color_field)) +
  ggplot2::geom_point(size=pt_size, alpha = sc.alpha) +
  ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2)))
fig  


plegend = ggpubr::get_legend(ggplot2::ggplot(scCoord,
                                             ggplot2::aes_string(colnames(scCoord)[1],
                                                                 colnames(scCoord)[2],
                                                                 color = color_field)) +
                               ggplot2::geom_point(size=pt_size, alpha = 1) +
                               ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 2))) + 
                               ggplot2::scale_color_manual(name = 'celltype', values = cpalette))

fig <- fig + 
  ggplot2::geom_point(data=centroids,
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2],
                                          fill = color_field, size = metric),
                      pch=21, color = 'black', stroke = 1, alpha = mc.alpha) +
  scale_size_continuous(limits = c(1, max_mcsize), range = c(1,12))
fig


fig <- fig + 
  ggplot2::geom_point(data=centroids[!is.na(centroids$dub_mc_test),],
                      ggplot2::aes_string(colnames(centroids)[1 + 1],
                                          colnames(centroids)[1 + 2],
                                          fill = color_field, size = metric),
                      pch=21, color = 'red', stroke = 1.5, alpha = 1) +
  scale_size_continuous(limits = c(1, max_mcsize), range = c(1,12))


fig <- fig + 
  ggplot2::scale_fill_manual(values = cpalette) +
  ggplot2::scale_color_manual(values = cpalette) +
  ggplot2::theme_classic() +
  ggplot2::theme(legend.position = 'none')

fig  # Figure 1c, top 


###

library('ComplexHeatmap')

mcid1 = 'mc30-452' 
mcid2 = 'mc30-86'

sub_sc = obj_singlecell[, sc_membership == mcid1]

sub_sc = FindVariableFeatures(sub_sc)
top100 = VariableFeatures(sub_sc)[1:50]

#
sub_sc1 = obj_singlecell[, sc_membership == mcid2]

dat_sc1 = GetAssayData(sub_sc1, layer = 'data')
dat_sc1 = t(as.matrix(dat_sc1[top100, ]))

sdd = apply(dat_sc1, 2, sd)

top100 = top100[sdd>0]

#
dat_sc = GetAssayData(sub_sc, layer = 'data')
dat_sc = t(as.matrix(dat_sc[top100, ]))

hm_dat =  Heatmap(dat_sc, 
                  # cluster_rows = F, cluster_columns = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  col = colorRamp2(c(0, 5), c('white','red')),
                  show_row_dend = T, show_column_dend = F,
                  show_heatmap_legend = F,
                  show_row_names = F, show_column_names = T)
hm_dat_ord = draw(hm_dat)
gene_colord = column_order(hm_dat_ord)
gene_roword = row_order(hm_dat_ord)

cor_sc = cor(dat_sc)
cor_sc[is.na(cor_sc)] = 0

hm_sc = Heatmap(cor_sc, 
                cluster_rows = F, cluster_columns = F,
                row_order = gene_colord, column_order = gene_colord,
                row_names_gp = gpar(fontsize = 5),
                column_names_gp = gpar(fontsize = 5),
                col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                show_row_dend = F, show_column_dend = F,
                show_heatmap_legend = F,
                show_row_names = F, show_column_names = F,
                heatmap_legend_param = list(title = 'corr'))
hm_sc_ord = draw(hm_sc)
gene_ord = row_order(hm_sc_ord)

hm_dat; hm_sc   # Figure 1c, bottom right

#
dat_sc1 = GetAssayData(sub_sc1, layer = 'data')
dat_sc1 = t(as.matrix(dat_sc1[top100, ]))

cor_sc1 = cor(dat_sc1)
cor_sc1[is.na(cor_sc1)] = 0

hm_sc1 = Heatmap(cor_sc1, 
                 cluster_rows = F, cluster_columns = F,
                 row_order = gene_colord, column_order = gene_colord,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 col = colorRamp2(c(-1,0,1), c('blue','white','red')), 
                 show_heatmap_legend = F,
                 show_row_names = F, show_column_names = F,
                 heatmap_legend_param = list(title = 'corr'))

hm_dat1 =  Heatmap(dat_sc1, 
                   # cluster_rows = F, 
                   cluster_columns = F,
                   column_order = gene_colord,
                   row_names_gp = gpar(fontsize = 5),
                   column_names_gp = gpar(fontsize = 5),
                   col = colorRamp2(c(0, 5), c('white','red')),
                   show_row_dend = T, show_column_dend = F,
                   show_heatmap_legend = F,
                   show_row_names = F, show_column_names = F)

hm_dat1; hm_sc1   # Figure 1c, bottom left
