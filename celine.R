
library(Seurat)
library(mcRigor)

data_path <- './celine_data/'


###### Data loading

cell_cycle_genes = read.table(file = paste0(data_path, 'regev_lab_cell_cycle_genes.txt'))

s_genes = cell_cycle_genes[[1]][1:43]
g2m_genes = cell_cycle_genes[[1]][44:97] 

load(paste0(data_path, "sincell_with_class_5cl.RData"))
obj = as.Seurat(sce_sc_10x_5cl_qc)

obj_singlecell_A549 = readRDS(file = paste0(data_path, '/A549_sce.rds'))
obj_singlecell_A549 = as.Seurat(obj_singlecell_A549)

features = union(rownames(obj_singlecell_A549), c(s_genes, g2m_genes))
features = intersect(features, rownames(obj))

obj_singlecell_A549 = subset(obj, cell_line == 'A549' & demuxlet_cls == 'SNG')
obj_singlecell_A549 = RenameAssays(obj_singlecell_A549, assay.name = names(obj_singlecell_A549@assays)[1], 
                                   new.assay.name = 'RNA')
obj_singlecell_A549 = obj_singlecell_A549[features,]

obj_singlecell_HCC827 = subset(obj, cell_line == 'HCC827' & demuxlet_cls == 'SNG')
obj_singlecell_HCC827 = RenameAssays(obj_singlecell_HCC827, assay.name = names(obj_singlecell_HCC827@assays)[1], 
                                     new.assay.name = 'RNA')
obj_singlecell_HCC827 = obj_singlecell_HCC827[features,]

obj_singlecell_H2228 = subset(obj, cell_line == 'H2228' & demuxlet_cls == 'SNG')
obj_singlecell_H2228 = RenameAssays(obj_singlecell_H2228, assay.name = names(obj_singlecell_H2228@assays)[1], 
                                    new.assay.name = 'RNA')
obj_singlecell_H2228 = obj_singlecell_H2228[features,]

##
library(AnnotationHub)
library(ensembldb)
library(magrittr)

ah <- AnnotationHub()
ahDb <- query(ah, 
              pattern = c("Homo sapiens", "EnsDb"), 
              ignore.case = TRUE)
id <- ahDb %>% mcols() %>% rownames() %>% tail(n = 1)
edb <- ah[[id]]

annotations <- genes(edb, 
                     return.type = "data.frame")
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, seq_name, gene_biotype, description)

#
obj_singlecell_jukart = readRDS(file = paste0(data_path, 'jukart_sce.rds'))
obj_singlecell_jukart = as.Seurat(obj_singlecell_jukart)

counts = GetAssayData(obj_singlecell_jukart, layer = 'counts')
counts = as.matrix(counts)
counts = counts[rownames(counts) %in% annotations$gene_id, ]
rownames(counts) = annotations$gene_name[match(rownames(counts), annotations$gene_id)]
counts = counts[rownames(counts) != '',]
features = union(rownames(counts), c(s_genes, g2m_genes))

counts_jukart = Read10X(data.dir = paste0(data_path, 'jukart'))
obj_singlecell_jukart = CreateSeuratObject(counts = counts_jukart)
features = intersect(features, rownames(obj_singlecell_jukart))
obj_singlecell_jukart = obj_singlecell_jukart[features,]

#
obj_singlecell_hek293t = readRDS(file = paste0(data_path, 'hek293t_sce.rds'))
obj_singlecell_hek293t = as.Seurat(obj_singlecell_hek293t)

counts = GetAssayData(obj_singlecell_hek293t, layer = 'counts')
counts = as.matrix(counts)
counts = counts[rownames(counts) %in% annotations$gene_id, ]
rownames(counts) = annotations$gene_name[match(rownames(counts), annotations$gene_id)]
counts = counts[rownames(counts) != '',]
features = union(rownames(counts), c(s_genes, g2m_genes))

counts_hek293t = Read10X(data.dir = paste0(data_path, 'hek293t'))
obj_singlecell_hek293t = CreateSeuratObject(counts = counts_hek293t)
features = intersect(features, rownames(obj_singlecell_hek293t))
obj_singlecell_hek293t = obj_singlecell_hek293t[features,]



###### Identify cell cycle phases & other preprocessing

obj_singlecell_A549 = CellCycleScoring(obj_singlecell_A549, s.features = s_genes, g2m.features = g2m_genes)
obj_singlecell_HCC827 = CellCycleScoring(obj_singlecell_HCC827, s.features = s_genes, g2m.features = g2m_genes)
obj_singlecell_H2228 = CellCycleScoring(obj_singlecell_H2228, s.features = s_genes, g2m.features = g2m_genes)
obj_singlecell_jukart = CellCycleScoring(obj_singlecell_jukart, s.features = s_genes, g2m.features = g2m_genes)
obj_singlecell_hek293t = CellCycleScoring(obj_singlecell_hek293t, s.features = s_genes, g2m.features = g2m_genes)

#
obj_singlecell_A549 <- NormalizeData(obj_singlecell_A549)
obj_singlecell_A549 <- FindVariableFeatures(obj_singlecell_A549, nfeatures = 2000)
obj_singlecell_A549 <- ScaleData(obj_singlecell_A549)
obj_singlecell_A549 = RunPCA(obj_singlecell_A549)
ElbowPlot(obj_singlecell_A549)
obj_singlecell_A549 = RunUMAP(obj_singlecell_A549, dims = 1:50)
UMAPPlot(obj_singlecell_A549, group.by='Phase')

obj_singlecell_HCC827 <- NormalizeData(obj_singlecell_HCC827)
obj_singlecell_HCC827 <- FindVariableFeatures(obj_singlecell_HCC827, nfeatures = 2000)
obj_singlecell_HCC827 <- ScaleData(obj_singlecell_HCC827)
obj_singlecell_HCC827 = RunPCA(obj_singlecell_HCC827)
ElbowPlot(obj_singlecell_HCC827)
obj_singlecell_HCC827 = RunUMAP(obj_singlecell_HCC827, dims = 1:50)
UMAPPlot(obj_singlecell_HCC827, group.by='Phase')

obj_singlecell_H2228 <- NormalizeData(obj_singlecell_H2228)
obj_singlecell_H2228 <- FindVariableFeatures(obj_singlecell_H2228, nfeatures = 2000)
obj_singlecell_H2228 <- ScaleData(obj_singlecell_H2228)
obj_singlecell_H2228 = RunPCA(obj_singlecell_H2228)
ElbowPlot(obj_singlecell_H2228)
obj_singlecell_H2228 = RunUMAP(obj_singlecell_H2228, dims = 1:50)
UMAPPlot(obj_singlecell_H2228, group.by='Phase')

obj_singlecell_jukart <- NormalizeData(obj_singlecell_jukart)
obj_singlecell_jukart <- FindVariableFeatures(obj_singlecell_jukart, nfeatures = 2000)
obj_singlecell_jukart <- ScaleData(obj_singlecell_jukart)
obj_singlecell_jukart = RunPCA(obj_singlecell_jukart)
ElbowPlot(obj_singlecell_jukart)
obj_singlecell_jukart = RunUMAP(obj_singlecell_jukart, dims = 1:50)
UMAPPlot(obj_singlecell_jukart, group.by='Phase')

obj_singlecell_hek293t <- NormalizeData(obj_singlecell_hek293t)
obj_singlecell_hek293t <- FindVariableFeatures(obj_singlecell_hek293t, nfeatures = 2000)
obj_singlecell_hek293t <- ScaleData(obj_singlecell_hek293t)
obj_singlecell_hek293t = RunPCA(obj_singlecell_hek293t)
ElbowPlot(obj_singlecell_hek293t)
obj_singlecell_hek293t = RunUMAP(obj_singlecell_hek293t, dims = 1:50)
UMAPPlot(obj_singlecell_hek293t, group.by='Phase')



###### Run mcRigor

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_A549 = read.csv(paste0(data_path,
                                  "seacells_cell_membership_A549.csv"), check.names = F, row.names = 1)
res_A549 = mcRigor_DETECT(obj_singlecell = obj_singlecell_A549, cell_membership = cell_membership_A549, 
                     output_file = paste0(data_path, 'seacells_tabmc_RNA_A549.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_H2228 = read.csv(paste0(data_path,
                                       "seacells_cell_membership_H2228.csv"), check.names = F, row.names = 1)
res_H2228 = mcRigor_DETECT(obj_singlecell = obj_singlecell_H2228, cell_membership = cell_membership_H2228, 
                          output_file = paste0(data_path, 'seacells_tabmc_RNA_H2228.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_HCC827 = read.csv(paste0(data_path,
                                       "seacells_cell_membership_HCC827.csv"), check.names = F, row.names = 1)
res_HCC827 = mcRigor_DETECT(obj_singlecell = obj_singlecell_HCC827, cell_membership = cell_membership_HCC827, 
                          output_file = paste0(data_path, 'seacells_tabmc_RNA_HCC827.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_hek293t = read.csv(paste0(data_path,
                                       "seacells_cell_membership_hek293t.csv"), check.names = F, row.names = 1)
res_hek293t = mcRigor_DETECT(obj_singlecell = obj_singlecell_hek293t, cell_membership = cell_membership_hek293t, 
                          output_file = paste0(data_path, 'seacells_tabmc_RNA_hek293t.RData'))   # intermediate results saved

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership_jukart = read.csv(paste0(data_path,
                                       "seacells_cell_membership_jukart.csv"), check.names = F, row.names = 1)
res_jukart = mcRigor_DETECT(obj_singlecell = obj_singlecell_jukart, cell_membership = cell_membership_jukart, 
                          output_file = paste0(data_path, 'seacells_tabmc_RNA_jukart.RData'))   # intermediate results saved




###### Figure 1d & Supplementary Figure 3

load(paste0(data_path, "seacells_tabmc_RNA_A549.RData"))
seacells_res_A549 = mcRigor_threshold(TabMC = TabMC[TabMC$gamma<=50,], test_cutoff = 0.05)

cell_membership <- read.csv(paste0(data_path, "seacells_cell_membership_rna_A549.csv"), check.names = F, row.names = 1)

tgamma = 30

sc_membership = cell_membership[[as.character(tgamma)]]
names(sc_membership) = rownames(cell_membership)

obj_mc_seacells_A549 = mcRigor_buildmc(obj_singlecell = obj_singlecell_A549, sc_membership = sc_membership,
                                       aggregate_method = 'mean',
                                       add_testres = T, 
                                       test_stats = seacells_res_A549$TabMC, Thre = seacells_res_A549$threshold)

ggplot(data = obj_mc_seacells_A549[[]], aes(x = mcRigor, y=Phase_purity)) + geom_violin()   # Supplementary Figure 3a

#
ggene = 'HMGB2' # 'AURKB' 'GINS2' 'MCMC6'

ggene_tab = data.frame(expr = GetAssayData(obj_mc_seacells_A549, layer = 'counts')[ggene,],
                       mcRigor = obj_mc_seacells_A549$mcRigor,
                       phase = obj_mc_seacells_A549$Phase)
ggene_dat = data.frame(expr = c(ggene_tab$expr[ggene_tab$mcRigor == 'trustworthy'],
                                ggene_tab$expr,
                                GetAssayData(obj_singlecell_A549, layer = 'counts')[ggene,]),
                       group = c(rep('trustworthy', length(which(ggene_tab$mcRigor == 'trustworthy'))),
                                 rep('all', length(ggene_tab$expr)),
                                 rep('sc', ncol(obj_singlecell_A549))),
                       phase = c(paste0(ggene_tab$phase[ggene_tab$mcRigor == 'trustworthy'], ' (trustworthy)'),
                                 paste0(ggene_tab$phase, ' (all)'),
                                 paste0(obj_singlecell_A549$Phase, ' (single cell)')))
ggene_dat$phase = factor(ggene_dat$phase, 
                         levels = c('G1 (single cell)', 'S (single cell)', 'G2M (single cell)',
                                    'G1 (all)', 'S (all)', 'G2M (all)', 
                                    'G1 (trustworthy)', 'S (trustworthy)', 'G2M (trustworthy)'))
ggene_dat$expr = log10(ggene_dat$expr + 1)


p = ggplot(data = ggene_dat, 
            aes(x = phase, y = expr, fill = phase)) +
  geom_violin() + 
  scale_fill_manual(values = c('#4F1309', '#015001', '#012962',
                               '#9E2611', '#348934', '#3474cf', 
                               '#FF7C6D', '#74f689', '#61e2ff')) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p   # Figure 1d & Supplementary Figure 3b





