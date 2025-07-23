
library(Seurat)
library(mcRigor)
library('ComplexHeatmap')
library(circlize)

data_path <- './covid19pbmc_data/'


###### Data loading

pbmc <- readRDS(paste0(data_path, 'blish_covid.seu.rds'))
pbmc <- UpdateSeuratObject(object = pbmc)
DefaultAssay(pbmc) <- 'RNA'

pbmc_B = pbmc[,pbmc$cell.type.coarse %in% 'B']
mean_exp = rowMeans(pbmc_B@assays$RNA@counts/pbmc_B$nCount_RNA)
genes_selected = names(sort.int(mean_exp, decreasing = T))[1:5000]

obj_singlecell = pbmc_B


###### Quality control & Preprocessing

obj_singlecell <- NormalizeData(obj_singlecell)
obj_singlecell <- FindVariableFeatures(obj_singlecell, nfeatures = 2000)
obj_singlecell <- ScaleData(obj_singlecell, features = rownames(obj_singlecell))

celltype_colors <- c("#7E57C2", "#1E88E5", "#FFC107", "#004D40", "#9E9D24", 
                     "#F06292", "#546E7A", "#D4E157", "#76FF03", "#6D4C41",
                     "#26A69A", "#AB47BC", "#EC407A", "#D81B60", "#42A5F5",
                     "#2E7D32", "#FFA726", "#5E35B1", "#EF5350", "#3949AB",
                     'darkblue', 'yellow')
names(celltype_colors) <- unique(obj_singlecell$cell.type)

obj_singlecell <- RunPCA(obj_singlecell)

UMAPPlot(obj_singlecell, group.by = 'cell.type', cols = celltype_colors)

pbmc_B_healthy = subset(obj_singlecell, subset = Status == 'Healthy')

UMAPPlot(pbmc_B_healthy)

pbmc_B_covid19 = subset(obj_singlecell, subset = Status == 'COVID')

UMAPPlot(pbmc_B_covid19)

saveRDS(pbmc_B_healthy, 
        file = paste0(data_path, "covid_pbmc_B_healthy.rds"))   # intermediate results saved

saveRDS(pbmc_B_covid19, 
        file = paste0(data_path, "covid_pbmc_B_covid19.rds"))   # intermediate results saved




###### Figure 1e & Supplementary Figure 4

obj_singlecell_covid19 = readRDS(file = paste0(data_path, "covid_pbmc_B_covid19.rds"))

obj_singlecell_healthy = readRDS(file = paste0(data_path, "covid_pbmc_B_healthy.rds"))

gamma = '30'
method = 'supercell'

{
  sub_sc = subset(obj_singlecell_covid19)
  Idents(sub_sc) = 'cell.type'
  
  gene_names = c('CD74','HLA-DRA','HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DRB5',
                 'HLA-DQA1','HLA-DMB','HLA-DQB1','HLA-DMA','HLA-DOA','HLA-DQA2',
                 'FCGR2B','HLA-DQB2',
                 'IGHM','IGLC3','IGKV4-1','IGLV2-14','IGHV4-59','CCR2','IGLV2-8',
                 'LAX1','IGHV4-39','IGHV4-61','JAK2','IGLV2-23','IGHV4-31',
                 'ADAR','EIF2AK2','MX2','IFIT3')
  all(gene_names %in% rownames(sub_sc))
  
  dat_sc = GetAssayData(sub_sc, layer = 'data')
  dat_sc = t(as.matrix(dat_sc[gene_names, ]))
  
  cor_sc = cor(dat_sc)
  
  cor_sc = ifelse(cor_sc>=0, cor_sc^(7/5), -abs(cor_sc)^(7/5))  
  
  hm_sc = Heatmap(cor_sc, 
                  cluster_rows = F, cluster_columns = F,
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                  heatmap_legend_param = list(title = 'corr'))
  hm_sc_ord = draw(hm_sc)
  gene_ord = row_order(hm_sc_ord)
  
  #
  
  obj_singlecell_healthy <- readRDS(paste0(data_path, 'covid_pbmc_B_healthy.rds'))
  
  sub_sc_healthy = subset(obj_singlecell_healthy)
  Idents(sub_sc_healthy) = 'cell.type'
  
  dat_sc_healthy = GetAssayData(sub_sc_healthy, layer = 'data')
  dat_sc_healthy = t(as.matrix(dat_sc_healthy[gene_names, ]))
  
  cor_sc_healthy = cor(dat_sc_healthy)
  
  cor_sc_healthy = ifelse(cor_sc_healthy>=0, cor_sc_healthy^(7/5), -abs(cor_sc_healthy)^(7/5))
  
  hm_sc_healthy = Heatmap(cor_sc_healthy, 
                          cluster_rows = F, cluster_columns = F,
                          show_column_names = F,
                          row_names_gp = gpar(fontsize = 5),
                          column_names_gp = gpar(fontsize = 5),
                          col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                          heatmap_legend_param = list(title = 'corr'))
  hm_sc_ord_healthy = draw(hm_sc_healthy)
  gene_ord_healthy = row_order(hm_sc_ord_healthy)
  
  hm_sc_healthy + hm_sc
  
  ##
  # obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
  cell_membership_covid19 <- read.csv(paste0(data_path, method, "_cell_membership_rna_covid19.csv"), check.names = F, row.names = 1)
  
  res = mcRigor_DETECT(obj_singlecell = obj_singlecell_covid19, cell_membership = cell_membership_covid19, 
                       test_cutoff = 0.05,
                       output_file = paste0(data_path, method, "_tabmc_RNA_covid19.RData"))   # intermediate results saved
  
  load(paste0(data_path, method, "_tabmc_RNA_covid19.RData"))
  TabMC_covid19 = TabMC
  
  sc_membership = cell_membership_covid19[[gamma]]
  names(sc_membership) = rownames(cell_membership_covid19)
  
  obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_covid19, sc_membership = sc_membership,
                                 add_testres = T, test_stats = TabMC_covid19, test_cutoff = 0.05,
                                 aggregate_method = 'sum')
  
  sub_mc = subset(obj_metacell)
  
  dat_mc <- GetAssayData(sub_mc, layer = 'data')
  dat_mc = t(as.matrix(dat_mc[gene_names,]))
  cor_mc = cor(dat_mc)
  
  cor_mc = ifelse(cor_mc>=0, cor_mc^(7/5), -abs(cor_mc)^(7/5)) 
  
  hm_mc = Heatmap(cor_mc, 
                  cluster_rows = F, cluster_columns = F,
                  # row_order = gene_ord, column_order = gene_ord,
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                  show_heatmap_legend = F)
  hm_sc + hm_mc
  
  #
  obj_metacell$dubmc = ifelse(obj_metacell$mcRigor == 'dubious', T, F)
  
  sub_mc_trust = subset(obj_metacell, subset = dubmc == F)
  
  dat_mc_trust <- GetAssayData(sub_mc_trust, layer = 'data')
  dat_mc_trust = t(as.matrix(dat_mc_trust[gene_names,]))
  cor_mc_trust = cor(dat_mc_trust)
  
  cor_mc_trust = ifelse(cor_mc_trust>=0, cor_mc_trust^(7/5), -abs(cor_mc_trust)^(7/5))
  
  hm_mc_trust = Heatmap(cor_mc_trust, 
                        cluster_rows = F, cluster_columns = F,
                        # row_order = gene_ord, column_order = gene_ord,
                        show_column_names = F,
                        row_names_gp = gpar(fontsize = 5),
                        column_names_gp = gpar(fontsize = 5),
                        col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                        show_heatmap_legend = F)
  hm_mc_trust
  hm_sc + hm_mc_trust
  
  
  ##
  # obtained by SuperCell (see "Implementing metacell partitioning methods" session in our online tutorial)
  cell_membership_healthy <- read.csv(paste0(data_path, method, "_cell_membership_rna_healthy.csv"), check.names = F, row.names = 1)
  
  res = mcRigor_DETECT(obj_singlecell = obj_singlecell_healthy, cell_membership = cell_membership_healthy, 
                       test_cutoff = 0.05,
                       output_file = paste0(data_path, method, "_tabmc_RNA_healthy.RData"))   # intermediate results saved
  
  load(paste0(data_path, method, "_tabmc_RNA_healthy.RData"))
  TabMC_healthy = TabMC
  
  sc_membership = cell_membership_healthy[[gamma]]
  names(sc_membership) = rownames(cell_membership_healthy)
  
  obj_metacell_healthy = mcRigor_buildmc(obj_singlecell = obj_singlecell_healthy, sc_membership = sc_membership,
                                         add_testres = T, test_stats = TabMC_healthy, test_cutoff = 0.05,
                                         aggregate_method = 'sum')
  
  sub_mc_healthy = subset(obj_metacell_healthy)
  
  dat_mc_healthy <- GetAssayData(sub_mc_healthy, layer = 'data')
  dat_mc_healthy = t(as.matrix(dat_mc_healthy[gene_names,]))
  cor_mc_healthy = cor(dat_mc_healthy)
  
  cor_mc_healthy = ifelse(cor_mc_healthy>=0, cor_mc_healthy^(7/5), -abs(cor_mc_healthy)^(7/5))
  
  hm_mc_healthy = Heatmap(cor_mc_healthy, 
                          cluster_rows = F, cluster_columns = F,
                          # row_order = gene_ord, column_order = gene_ord,
                          show_column_names = F,
                          row_names_gp = gpar(fontsize = 5),
                          column_names_gp = gpar(fontsize = 5),
                          col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                          show_heatmap_legend = F)
  hm_sc_healthy + hm_mc_healthy
  
  #
  obj_metacell_healthy$dubmc = ifelse(obj_metacell_healthy$mcRigor == 'dubious', T, F)
  
  sub_mc_trust_healthy = subset(obj_metacell_healthy, subset = dubmc == F)
  
  dat_mc_trust_healthy <- GetAssayData(sub_mc_trust_healthy, layer = 'data')
  dat_mc_trust_healthy = t(as.matrix(dat_mc_trust_healthy[gene_names,]))
  cor_mc_trust_healthy = cor(dat_mc_trust_healthy)
  
  cor_mc_trust_healthy = ifelse(cor_mc_trust_healthy>=0, cor_mc_trust_healthy^(7/5), -abs(cor_mc_trust_healthy)^(7/5))
  
  hm_mc_trust_healthy = Heatmap(cor_mc_trust_healthy, 
                                cluster_rows = F, cluster_columns = F,
                                # row_order = gene_ord, column_order = gene_ord,
                                show_column_names = F,
                                row_names_gp = gpar(fontsize = 5),
                                column_names_gp = gpar(fontsize = 5),
                                col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                                show_heatmap_legend = F)
  hm_mc_trust_healthy
  hm_sc_healthy + hm_mc_trust_healthy
  
}

hm_sc_healthy + hm_mc_trust_healthy + hm_mc_healthy   # Figure 1e, top left

hm_sc + hm_mc_trust + hm_mc   # Figure 1e, bottom left


#
gg1 = 'IGHV4-61'
gg2 = 'IGKV4-1' 
gene1 = GetAssayData(obj_singlecell_healthy, layer = 'data')[gg1,]
gene2 = GetAssayData(obj_singlecell_healthy, layer = 'data')[gg2,]

pointdat = data.frame(gene1 = gene1,
                      gene2 = gene2,
                      celltype = obj_singlecell_healthy$cell.type)

ggplot(pointdat, aes(x=gene1, y=gene2)) + 
  geom_point(size = 1, col = '#811EFF') +
  xlab('') + ylab('') 

mcgene1 = GetAssayData(obj_metacell_healthy, layer = 'data')[gg1,]
mcgene2 = GetAssayData(obj_metacell_healthy, layer = 'data')[gg2,]

mcpointdat = data.frame(gene1 = mcgene1,
                        gene2 = mcgene2,
                        celltype = obj_metacell_healthy$cell.type,
                        dubmc = obj_metacell_healthy$dubmc) 


mcdid = which(mcpointdat$dubmc)

mcpp = ggplot(data = mcpointdat[mcdid,], aes(x=gene1, y=gene2)) + 
  geom_point(mapping = aes(x=gene1, y=gene2), 
             col = '#F11E6B', size = 3) +
  geom_point(data = mcpointdat[-mcdid,], size = 3, col = '#1229D9') +
  geom_point(data = mcpointdat[mcdid,], col = '#F11E6B', size = 3) +
  xlab(gg1) + ylab(gg2) 
mcpp   # Figure 1e, right





###### Supplementary Figure 4

obj_singlecell_covid19 = readRDS(file = paste0(data_path, "covid_pbmc_B_covid19.rds"))

obj_singlecell_healthy = readRDS(file = paste0(data_path, "covid_pbmc_B_healthy.rds"))


###

gamma = '30'
method = 'seacells'   # method = "mc2", "mc1" for Supplementary Figure 4b-c

{
  sub_sc = subset(obj_singlecell_covid19)
  Idents(sub_sc) = 'cell.type'
  
  gene_names = c('CD74','HLA-DRA','HLA-DRB1','HLA-DPB1','HLA-DPA1','HLA-DRB5',
                 'HLA-DQA1','HLA-DMB','HLA-DQB1','HLA-DMA','HLA-DOA','HLA-DQA2',
                 'FCGR2B','HLA-DQB2',
                 'IGHM','IGLC3','IGKV4-1','IGLV2-14','IGHV4-59','CCR2','IGLV2-8',
                 'LAX1','IGHV4-39','IGHV4-61','JAK2','IGLV2-23','IGHV4-31',
                 'ADAR','EIF2AK2','MX2','IFIT3')
  all(gene_names %in% rownames(sub_sc))
  
  dat_sc = GetAssayData(sub_sc, layer = 'data')
  dat_sc = t(as.matrix(dat_sc[gene_names, ]))
  
  cor_sc = cor(dat_sc)
  
  cor_sc = ifelse(cor_sc>=0, cor_sc^(7/5), -abs(cor_sc)^(7/5))
  
  hm_sc = Heatmap(cor_sc, 
                  cluster_rows = F, cluster_columns = F,
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                  heatmap_legend_param = list(title = 'corr'))
  hm_sc_ord = draw(hm_sc)
  gene_ord = row_order(hm_sc_ord)
  
  #
  
  obj_singlecell_healthy <- readRDS(paste0(data_path, 'covid_pbmc_B_healthy.rds'))
  
  sub_sc_healthy = subset(obj_singlecell_healthy)
  Idents(sub_sc_healthy) = 'cell.type'
  
  dat_sc_healthy = GetAssayData(sub_sc_healthy, layer = 'data')
  dat_sc_healthy = t(as.matrix(dat_sc_healthy[gene_names, ]))
  
  cor_sc_healthy = cor(dat_sc_healthy)
  
  cor_sc_healthy = ifelse(cor_sc_healthy>=0, cor_sc_healthy^(7/5), -abs(cor_sc_healthy)^(7/5))
  
  hm_sc_healthy = Heatmap(cor_sc_healthy, 
                          cluster_rows = F, cluster_columns = F,
                          show_column_names = F,
                          row_names_gp = gpar(fontsize = 5),
                          column_names_gp = gpar(fontsize = 5),
                          col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                          heatmap_legend_param = list(title = 'corr'))
  hm_sc_ord_healthy = draw(hm_sc_healthy)
  gene_ord_healthy = row_order(hm_sc_ord_healthy)
  
  hm_sc_healthy + hm_sc
  
  ##
  # obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
  cell_membership_covid19 <- read.csv(paste0(data_path, method, "_cell_membership_rna_covid19.csv"), check.names = F, row.names = 1)
  
  res = mcRigor_DETECT(obj_singlecell = obj_singlecell_covid19, cell_membership = cell_membership_covid19, 
                       test_cutoff = 0.05,
                       output_file = paste0(data_path, method, "_tabmc_RNA_covid19.RData"))   # intermediate results saved
  
  load(paste0(data_path, method, "_tabmc_RNA_covid19.RData"))
  TabMC_covid19 = TabMC
  
  sc_membership = cell_membership_covid19[[gamma]]
  names(sc_membership) = rownames(cell_membership_covid19)
  
  obj_metacell = mcRigor_buildmc(obj_singlecell = obj_singlecell_covid19, sc_membership = sc_membership,
                                 add_testres = T, test_stats = TabMC_covid19, test_cutoff = 0.05,
                                 aggregate_method = 'sum')
  
  sub_mc = subset(obj_metacell)
  
  dat_mc <- GetAssayData(sub_mc, layer = 'data')
  dat_mc = t(as.matrix(dat_mc[gene_names,]))
  cor_mc = cor(dat_mc)
  
  cor_mc = ifelse(cor_mc>=0, cor_mc^(7/5), -abs(cor_mc)^(7/5)) 
  
  hm_mc = Heatmap(cor_mc, 
                  cluster_rows = F, cluster_columns = F,
                  # row_order = gene_ord, column_order = gene_ord,
                  show_column_names = F,
                  row_names_gp = gpar(fontsize = 5),
                  column_names_gp = gpar(fontsize = 5),
                  col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                  show_heatmap_legend = F)
  hm_sc + hm_mc
  
  #
  obj_metacell$dubmc = ifelse(obj_metacell$mcRigor == 'dubious', T, F)
  
  sub_mc_trust = subset(obj_metacell, subset = dubmc == F)
  
  dat_mc_trust <- GetAssayData(sub_mc_trust, layer = 'data')
  dat_mc_trust = t(as.matrix(dat_mc_trust[gene_names,]))
  cor_mc_trust = cor(dat_mc_trust)
  
  cor_mc_trust = ifelse(cor_mc_trust>=0, cor_mc_trust^(7/5), -abs(cor_mc_trust)^(7/5))
  
  hm_mc_trust = Heatmap(cor_mc_trust, 
                        cluster_rows = F, cluster_columns = F,
                        # row_order = gene_ord, column_order = gene_ord,
                        show_column_names = F,
                        row_names_gp = gpar(fontsize = 5),
                        column_names_gp = gpar(fontsize = 5),
                        col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                        show_heatmap_legend = F)
  hm_mc_trust
  hm_sc + hm_mc_trust
  
  
  ##
  # obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
  cell_membership_healthy <- read.csv(paste0(data_path, method, "_cell_membership_rna_healthy.csv"), check.names = F, row.names = 1)
  
  res = mcRigor_DETECT(obj_singlecell = obj_singlecell_healthy, cell_membership = cell_membership_healthy, 
                       test_cutoff = 0.05,
                       output_file = paste0(data_path, method, "_tabmc_RNA_healthy.RData"))   # intermediate results saved
  
  load(paste0(data_path, method, "_tabmc_RNA_healthy.RData"))
  TabMC_healthy = TabMC
  
  sc_membership = cell_membership_healthy[[gamma]]
  names(sc_membership) = rownames(cell_membership_healthy)
  
  obj_metacell_healthy = mcRigor_buildmc(obj_singlecell = obj_singlecell_healthy, sc_membership = sc_membership,
                                         add_testres = T, test_stats = TabMC_healthy, test_cutoff = 0.05,
                                         aggregate_method = 'sum')
  
  sub_mc_healthy = subset(obj_metacell_healthy)
  
  dat_mc_healthy <- GetAssayData(sub_mc_healthy, layer = 'data')
  dat_mc_healthy = t(as.matrix(dat_mc_healthy[gene_names,]))
  cor_mc_healthy = cor(dat_mc_healthy)
  
  cor_mc_healthy = ifelse(cor_mc_healthy>=0, cor_mc_healthy^(7/5), -abs(cor_mc_healthy)^(7/5))
  
  hm_mc_healthy = Heatmap(cor_mc_healthy, 
                          cluster_rows = F, cluster_columns = F,
                          # row_order = gene_ord, column_order = gene_ord,
                          show_column_names = F,
                          row_names_gp = gpar(fontsize = 5),
                          column_names_gp = gpar(fontsize = 5),
                          col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                          show_heatmap_legend = F)
  hm_sc_healthy + hm_mc_healthy
  
  #
  obj_metacell_healthy$dubmc = ifelse(obj_metacell_healthy$mcRigor == 'dubious', T, F)
  
  sub_mc_trust_healthy = subset(obj_metacell_healthy, subset = dubmc == F)
  
  dat_mc_trust_healthy <- GetAssayData(sub_mc_trust_healthy, layer = 'data')
  dat_mc_trust_healthy = t(as.matrix(dat_mc_trust_healthy[gene_names,]))
  cor_mc_trust_healthy = cor(dat_mc_trust_healthy)
  
  cor_mc_trust_healthy = ifelse(cor_mc_trust_healthy>=0, cor_mc_trust_healthy^(7/5), -abs(cor_mc_trust_healthy)^(7/5))
  
  hm_mc_trust_healthy = Heatmap(cor_mc_trust_healthy, 
                                cluster_rows = F, cluster_columns = F,
                                # row_order = gene_ord, column_order = gene_ord,
                                show_column_names = F,
                                row_names_gp = gpar(fontsize = 5),
                                column_names_gp = gpar(fontsize = 5),
                                col = colorRamp2(c(-1,0,1), c('blue','white','red')),
                                show_heatmap_legend = F)
  hm_mc_trust_healthy
  hm_sc_healthy + hm_mc_trust_healthy
  
}

hm_sc_healthy + hm_mc_trust_healthy + hm_mc_healthy   # Supplementary Figure 4a, top middle

hm_sc + hm_mc_trust + hm_mc   # Supplementary Figure 4a, bottom middle


##
gg1 = 'IGHV4-61'  
gg2 = 'IGKV4-1' 
gene1 = GetAssayData(obj_singlecell_healthy, layer = 'data')[gg1,]
gene2 = GetAssayData(obj_singlecell_healthy, layer = 'data')[gg2,]

pointdat = data.frame(gene1 = gene1,
                      gene2 = gene2,
                      celltype = obj_singlecell_healthy$cell.type)

ggplot(pointdat, aes(x=gene1, y=gene2)) + 
  geom_point(size = 1, col = '#811EFF') +
  xlab('') + ylab('') 

mcgene1 = GetAssayData(obj_metacell_healthy, layer = 'data')[gg1,]
mcgene2 = GetAssayData(obj_metacell_healthy, layer = 'data')[gg2,]

mcpointdat = data.frame(gene1 = mcgene1,
                        gene2 = mcgene2,
                        celltype = obj_metacell_healthy$cell.type,
                        dubmc = obj_metacell_healthy$dubmc) 


mcdid = which(mcpointdat$dubmc)

mcpp = ggplot(data = mcpointdat[mcdid,], aes(x=gene1, y=gene2)) + 
  geom_point(mapping = aes(x=gene1, y=gene2), 
             col = '#F11E6B', size = 3) +
  geom_point(data = mcpointdat[-mcdid,], size = 3, col = '#1229D9') +
  geom_point(data = mcpointdat[mcdid,], col = '#F11E6B', size = 3) +
  xlab(gg1) + ylab(gg2) 
mcpp   # Supplementary Figure 4a, right


##
sc_membership = cell_membership_healthy[[gamma]]
names(sc_membership) = rownames(cell_membership_healthy)

pproj_healthy = mcRigor_projection(obj_singlecell = obj_singlecell_healthy,
                                   sc_membership = sc_membership,
                                   color_field = 'cell.type',
                                   add_testres = T, test_stats = TabMC, test_cutoff = 0.05,
                                   dub_mc_test.label = T)
pproj_healthy   # Supplementary Figure 4a, top left


#
sc_membership = cell_membership_covid19[[gamma]]
names(sc_membership) = rownames(cell_membership_covid19)

pproj_covid19 = mcRigor_projection(obj_singlecell = obj_singlecell_covid19,
                                   sc_membership = sc_membership,
                                   color_field = 'cell.type',
                                   add_testres = T, test_stats = TabMC, test_cutoff = 0.05,
                                   dub_mc_test.label = T)
pproj_covid19   # Supplementary Figure 4a, bottom left




