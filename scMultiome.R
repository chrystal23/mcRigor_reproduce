
library(Seurat)
library(mcRigor)

library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)

library(GenomicRanges)
library(SummarizedExperiment)

data_path <- './scMultiome_data/'



###### Data loading

anno = GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
anno
seqlevels(anno)
seqlevels(anno) <- paste0('chr', seqlevels(anno))

#
obj_singlecell_rna <- anndata::read_h5ad(paste0(data_path, 'cd34_multiome_rna.h5ad'))
obj_singlecell_rna <- CreateSeuratObject(counts = Matrix::t(obj_singlecell_rna$X),
                                         meta.data = obj_singlecell_rna$obs)
obj_singlecell_rna = NormalizeData(obj_singlecell_rna)
obj_singlecell_rna = FindVariableFeatures(obj_singlecell_rna)
obj_singlecell_rna = ScaleData(obj_singlecell_rna)
obj_singlecell_rna = RunPCA(obj_singlecell_rna)
obj_singlecell_rna = RunUMAP(obj_singlecell_rna, reduction = 'pca', dims = 1:30, seed.use = 123456)
UMAPPlot(obj_singlecell_rna, group.by = 'celltype')

#
obj_singlecell_atac <- anndata::read_h5ad(paste0(data_path, 'cd34_multiome_atac.h5ad'))
assay_atac <- CreateChromatinAssay(counts = Matrix::t(obj_singlecell_atac$X),
                                   sep = c(':', '-'))
obj_singlecell_atac <- CreateSeuratObject(counts = assay_atac,
                                          assay = 'peaks',
                                          meta.data = obj_singlecell_atac$obs)
obj_singlecell_atac = RunTFIDF(obj_singlecell_atac)
obj_singlecell_atac <- FindTopFeatures(obj_singlecell_atac, min.cutoff = 'q0')
obj_singlecell_atac <- RunSVD(obj_singlecell_atac)
obj_singlecell_atac <- RunUMAP(obj_singlecell_atac, reduction = 'lsi', dims = 2:30, seed.use = 123456)
UMAPPlot(obj_singlecell_atac, group.by = 'celltype')




###### Functions

prepare_multiome_obj <- function(atac_obj, rna_obj, sc_membership = NULL, 
                                 mc_by_assay = c('ATAC', 'RNA'),
                                 aggregate_reduction = T,
                                 test_cutoff = 0.01,
                                 n_bins_for_gc = 50) {
  
  mc_by_assay = match.arg(mc_by_assay)
  
  if (is.null(sc_membership)) stop('Please provide the metacell membership of single cells.')
  if (is.null(names(sc_membership))) warning('The metacell memberships are not matched to the single cell ids.')
  
  # Subset of cells common to ATAC and RNA
  common_cells <- intersect(colnames(atac_obj), colnames(rna_obj))
  if (length(common_cells) != dim(atac_obj)[2]) {
    cat("Warning: The number of cells in RNA and ATAC objects are different. Only the common cells will be used.\n")
  }
  atac_mod_obj <- atac_obj[, common_cells]
  rna_mod_obj <- rna_obj[, common_cells]
  sc_membership = sc_membership[colnames(atac_mod_obj)]
  
  # Generate metacell matrices: discard the single cells (outliers?)
  metacells <- unique(sc_membership)
  metacells <- metacells[table(sc_membership)[metacells] > 1]
  
  cat("Building Metacells...\n")
  cat(" ATAC\n")
  atac_mc_obj <- mcRigor_buildmc(obj_singlecell = atac_mod_obj, sc_membership = sc_membership,
                                 assay_type = 'ATAC', doNorm = F, aggregate_method = 'sum')
  atac_mc_obj <- atac_mc_obj[, metacells]
  
  if (aggregate_reduction) {
    svd <- atac_mod_obj@reductions$lsi@cell.embeddings
    summ_svd <- aggregate(svd, list(sc_membership), mean)
    rownames(summ_svd) = summ_svd[,1]
    atac_mc_lsi = Seurat::CreateDimReducObject(embeddings = as.matrix(summ_svd[colnames(atac_mc_obj), -c(1)]),
                                               assay = 'peaks')
    atac_mc_obj@reductions$lsi = atac_mc_lsi
    
    umap <- atac_mod_obj@reductions$umap@cell.embeddings
    summ_umap <- aggregate(umap, list(sc_membership), mean)
    rownames(summ_umap) = summ_umap[,1]
    atac_mc_umap = Seurat::CreateDimReducObject(embeddings = as.matrix(summ_umap[colnames(atac_mc_obj), -c(1)]),
                                                assay = 'peaks')
    atac_mc_obj@reductions$umap = atac_mc_umap
  }
  
  atac_mc_obj = Signac::RunTFIDF(atac_mc_obj)
  
  cat(" RNA\n")
  rna_mc_obj <- mcRigor_buildmc(obj_singlecell = rna_mod_obj, sc_membership = sc_membership,
                                assay_type = 'RNA', doNorm = F, aggregate_method = 'sum')
  rna_mc_obj <- rna_mc_obj[, metacells]
  rna_mc_obj <- Seurat::NormalizeData(rna_mc_obj)
  
  if (aggregate_reduction) {
    pca <- rna_mod_obj@reductions$pca@cell.embeddings
    summ_pca <- aggregate(pca, list(sc_membership), mean)
    rownames(summ_pca) = summ_pca[,1]
    rna_mc_pca = Seurat::CreateDimReducObject(embeddings = as.matrix(summ_pca[colnames(rna_mc_obj), -c(1)]),
                                              assay = 'RNA')
    rna_mc_obj@reductions$pca = rna_mc_pca
    
    umap <- rna_mod_obj@reductions$umap@cell.embeddings
    summ_umap <- aggregate(umap, list(sc_membership), mean)
    rownames(summ_umap) = summ_umap[,1]
    rna_mc_umap = Seurat::CreateDimReducObject(embeddings = as.matrix(summ_umap[colnames(rna_mc_obj), -c(1)]),
                                               assay = 'RNA')
    rna_mc_obj@reductions$umap = rna_mc_umap
  }
  
  #
  test_res = mcRigor_threshold(TabMC, test_cutoff = test_cutoff)
  Thre = test_res$threshold
  
  atac_mc_obj$testres = 'trustworthy'
  for (mcid in colnames(atac_mc_obj)) {
    if (TabMC[mcid, 'size'] != atac_mc_obj@meta.data[mcid, 'size']) print(mcid)
    if (TabMC[mcid, 'size'] > 1 && 
        TabMC[mcid, 'TT_div'] > Thre$thre[Thre$size == TabMC[mcid, 'size']]) atac_mc_obj$testres[mcid] = 'dubious'
    
  }
  print(table(atac_mc_obj$testres))
  
  rna_mc_obj[['testres']] = atac_mc_obj[['testres']]
  
  return(list(atac_mc_obj = atac_mc_obj, rna_mc_obj = rna_mc_obj))
}



######

get_gene_peak_correlations <- function(atac_mc_obj, rna_mc_obj, gene_set = NULL,
                                       cor_by_layer = c('data', 'counts'),
                                       GeneAnno = NULL,
                                       span = 100000, 
                                       doPlot = F,
                                       n_jobs = 1) {
  
  cor_by_layer = match.arg(cor_by_layer)
  
  if (is.null(GeneAnno)) stop('Please provide gene annotations.')
  
  atac_exprs <- GetAssayData(atac_mc_obj, layer = cor_by_layer)
  rna_exprs <-GetAssayData(rna_mc_obj, layer = cor_by_layer)
  peaks_gr <- atac_mc_obj[['peaks']]@ranges
  names(peaks_gr) <- rownames(atac_mc_obj)
  
  if (is.null(gene_set)) {
    use_genes <- rownames(rna_mc_obj)
  } else {
    use_genes <- gene_set
  }
  
  gene_peak_correlations <- lapply(use_genes, function(gene) {
    peaks_correlations_per_gene(gene, atac_exprs, rna_exprs, atac_mc_obj, 
                                peaks_gr, GeneAnno, span, doPlot)
  })
  
  names(gene_peak_correlations) <- use_genes
  
  return(gene_peak_correlations)
}


peaks_correlations_per_gene <- function(gene, atac_exprs, rna_exprs, atac_mc_obj, 
                                        peaks_gr, GeneAnno, span = 100000, 
                                        doPlot = F,
                                        verbose = T, n_rand_sample = 100) {
  
  if (verbose) cat(gene, '== ')
  
  gene_anno <- GeneAnno[GeneAnno$gene_name == gene]
  if (length(gene_anno) == 0) {
    return(0)
  }
  longest_transcript <- gene_anno[which.max(width(gene_anno))]
  start <- start(longest_transcript) - span
  end <- end(longest_transcript) + span
  
  gene_gr <- GRanges(seqnames = seqnames(longest_transcript), ranges = IRanges(start = start, end = end))
  gene_peaks <- subsetByOverlaps(peaks_gr, gene_gr)
  
  if (length(gene_peaks) == 0) {
    return(0)
  }
  
  trust_id = colnames(atac_mc_obj)[atac_mc_obj$testres == 'trustworthy']
  
  X_trust <- atac_exprs[names(gene_peaks), trust_id, drop = F]
  rna_rank_trust = rank(rna_exprs[gene, trust_id])
  cors_trust <- apply(X_trust, 1, function(x) cor(rank(x), rna_rank_trust))
  # cors_trust <- apply(X_trust, 1, function(x) cor(x, rna_exprs[gene, trust_id]))
  names(cors_trust) <- rownames(X_trust)
  
  reg_peak_trust = names(cors_trust)[which.max(cors_trust)]
  
  if (doPlot) {
    pdat_trust = data.frame(atac = atac_exprs[reg_peak_trust,], rna = rna_exprs[gene,], 
                            celltype = rna_mc_obj$celltype, testres = atac_mc_obj$testres)
    splot_trust = ggplot(pdat_trust, aes(x=rna, y=atac)) + 
      geom_point(mapping = aes(col = celltype), size=2) + 
      geom_point(data = pdat_trust[pdat_trust$testres == 'dubious',], mapping = aes(x = rna, y = atac), 
                 col = 'black', shape = 1, size = 3)
    # splot_trust
  }
  
  X_all <- atac_exprs[names(gene_peaks), , drop = F]
  rna_rank_all = rank(rna_exprs[gene, ])
  cors_all <- apply(X_all, 1, function(x) cor(rank(x), rna_rank_all))
  # cors_all <- apply(X_all, 1, function(x) cor(x, rna_exprs[gene,]))
  names(cors_all) <- rownames(X_all)
  
  reg_peak_all = names(cors_all)[which.max(cors_all)]
  # cat(reg_peak_all, cors_all[reg_peak_all])
  
  if (doPlot) {
    pdat_all = data.frame(atac = atac_exprs[reg_peak_all,], rna = rna_exprs[gene,], 
                          celltype = rna_mc_obj$celltype, testres = atac_mc_obj$testres)
    splot_all = ggplot(pdat_all, aes(x=rna, y=atac)) + geom_point(mapping = aes(col = celltype), size=2) + 
      geom_point(data = pdat_all[pdat_all$testres == 'dubious',], mapping = aes(x = rna, y = atac), 
                 col = 'black', shape = 1, size = 3)
    scatter_plot = list(splot_trust = splot_trust, splot_all = splot_all)
    
  } else scatter_plot = NULL
  
  metric = abs(cors_trust - cors_all)
  metric = metric[cors_trust > 0.5 | cors_all > 0.5]
  
  quick_cors = c(cors_trust[reg_peak_trust], cors_all[reg_peak_all])
  names(quick_cors) = c('trust', 'all')
  
  return(list(corsDF = data.frame(trust = cors_trust, all = cors_all), 
              reg_peak = c(reg_peak_trust, reg_peak_all), 
              quick_cors = quick_cors,
              scatter_plot = scatter_plot,
              metric = max(metric)))
}


get_gene_peak_assocations <- function(corr_ranges, gene_pool = NULL, rna_mc_obj = NULL, 
                                      pval_cutoff = 0.1, cor_cutoff = 0.1) {
  
  if (is.null(gene_pool) && is.null(rna_mc_obj)) stop('Cannot determine which genes to use.')
  if (is.null(gene_pool)) {
    rna_mc_obj = FindVariableFeatures(rna_mc_obj, nfeatures = 1000)
    gene_pool = VariableFeatures(rna_mc_obj)
  }
  
  peak_counts <- sapply(gene_pool, function(g) {
    if (!(g %in% corr_ranges$gene)) 0
    gene_gr = corr_ranges[corr_ranges$gene == g]
    length(which(gene_gr$pvalue < pval_cutoff & gene_gr$score > cor_cutoff))
  })
  
  return(peak_counts)
}



get_gene_scores <- function(atac_mc_obj, corr_ranges, gene_pool = NULL, rna_mc_obj = NULL, 
                            pval_cutoff = 0.1, cor_cutoff = 0.1) {
  
  if (is.null(gene_pool) && is.null(rna_mc_obj)) stop('Cannot determine which genes to use.')
  if (is.null(gene_pool)) {
    rna_mc_obj = FindVariableFeatures(rna_mc_obj, nfeatures = 1000)
    gene_pool = VariableFeatures(rna_mc_obj)
  }
  
  gene_scores <- matrix(0.0, ncol = ncol(atac_mc_obj), nrow = length(gene_pool))
  rownames(gene_scores) <- gene_pool
  colnames(gene_scores) = colnames(atac_mc_obj)
  
  for (gene in rownames(gene_scores)) {
    if (!(gene %in% corr_ranges$gene)) next
    gene_gr <-  corr_ranges[corr_ranges$gene == gene]
    gene_peaks <- gene_gr$peak[gene_gr$pvalue < pval_cutoff & gene_gr$score > cor_cutoff]
    if (length(gene_peaks) == 0) next
    gene_scores[gene,] <- colSums(GetAssayData(atac_mc_obj[gene_peaks,], layer = 'data') * 
                                    gene_gr$score[gene_gr$peak %in% gene_peaks])
    cat(gene, '== ')
  }
  
  return(gene_scores)
}


splot_peakgene_corr <- function(atac_mc_obj, rna_mc_obj, gene, peak = NULL,
                                cor_by_layer = c('data', 'counts'),
                                corr_ranges) {
  
  cor_by_layer = match.arg(cor_by_layer)
  atac_exprs = GetAssayData(atac_mc_obj, layer = cor_by_layer)
  rna_exprs = GetAssayData(rna_mc_obj, layer = cor_by_layer)
  
  gene_ranges = corr_ranges[corr_ranges$gene == gene]
  
  if (is.null(peak)) {
    
    reg_peak_trust = gene_ranges$peak[which.max(gene_ranges$score)]
    reg_peak_trust1 = gene_ranges$peak[which.min(gene_ranges$pvalue)]
    
    pdat_trust = data.frame(atac = atac_exprs[reg_peak_trust,], atac1 = atac_exprs[reg_peak_trust1,],
                            rna = rna_exprs[gene,], 
                            celltype = rna_mc_obj$celltype, testres = atac_mc_obj$testres)
    
    gene_ranges_all = all_ranges[all_ranges$gene == gene]
    round(gene_ranges_all$score[gene_ranges_all$peak == reg_peak_trust], 4)
    
    splot_trust = ggplot(pdat_trust, aes(x=rna, y=atac)) + 
      geom_point(mapping = aes(col = celltype), size=2) + 
      geom_point(data = pdat_trust[pdat_trust$testres == 'dubious',], mapping = aes(x = rna, y = atac), 
                 # col = 'black', shape = 1, size = 5) +
                 col = 'red', shape = 1, size = 2.9, stroke = 1.5) +
      ggtitle(paste0(gene, '  ', round(gene_ranges$score[gene_ranges$peak == reg_peak_trust], 4), '\n ', reg_peak_trust)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 14))+
      ylab('ATAC') + xlab('RNA')
    
  } else {
    
    pdat_trust = data.frame(atac = atac_exprs[peak,],
                            rna = rna_exprs[gene,], 
                            celltype = rna_mc_obj$celltype, testres = atac_mc_obj$testres)
    
    splot_trust = ggplot(pdat_trust, aes(x=rna, y=atac)) + 
      geom_point(mapping = aes(col = celltype), size=2) + 
      geom_point(data = pdat_trust[pdat_trust$testres == 'dubious',], mapping = aes(x = rna, y = atac), 
                 col = 'black', shape = 1, size = 5) +
      ggtitle(paste0(gene, '  ', round(gene_ranges$score[gene_ranges$peak == peak],4), '\n ', peak)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 14))+
      ylab('ATAC') + xlab('RNA')
    
  }
  
  return(splot_trust)
  
}





###### Run metacell partitioning & mcRigor

# obtained by SEACells in python (see "Implementing metacell partitioning methods" session in our online tutorial)
cell_membership = read.csv(paste0(data_path, "seacells_cell_membership_atac.csv"), check.names = F, row.names = 1)

res = mcRigor_DETECT(obj_singlecell = obj_singlecell_atac, cell_membership = cell_membership, 
                     assay_type = 'ATAC',
                     output_file = paste0(data_path, 'seacells_tabmc_ATAC.RData'))   # intermediate results saved



###### Derive highly correlated peaks for each gene at single-cell and metacell level

load(paste0(data_path, "seacells_tabmc_ATAC.RData"))

gene_set = c( 'GATA2', 'GATA1', 'PRDX1', 'TAL1', 'MPO', 'IRF8','LYZ', 'KLF1',
              "AVP", "MIR181A1HG", "CRHBP" ,"PCDH9" ,"GPC5", "PRKG1",    
              "ANK1", "XACT", "RYR3", "TFR2", "SLC40A1", 'SPINK2')

sc_GenePeak_Signac = list()

both_mc_obj = obj_singlecell_atac
both_mc_obj[['RNA']] = obj_singlecell_rna[['RNA']]

Annotation(both_mc_obj) = anno
both_mc_obj <- RegionStats(both_mc_obj, genome = BSgenome.Hsapiens.UCSC.hg38)

both_mc_obj <- LinkPeaks(
  object = both_mc_obj,
  peak.assay = "peaks",
  expression.assay = "RNA",
  genes.use = gene_set,
  peak.slot = 'data',
  expression.slot = 'data',
  method = 'spearman',
  distance = 100000
)


GenePeakRes_Signac = list(all = Links(both_mc_obj))
sc_GenePeak_Signac[['1']] = GenePeakRes_Signac

save(sc_GenePeak_Signac, file = paste0(data_path, 'sc_GenePeak_Signac.RData'))   # intermediate results saved


##
seacells_GenePeak_Signac = list()

for(tgamma in as.character(c(90, 75, 50, 5))) {  
  
  named_membership = cell_membership[[tgamma]]
  names(named_membership) = rownames(cell_membership)
  
  #
  multiome_obj = prepare_multiome_obj(atac_obj = obj_singlecell_atac, rna_obj = obj_singlecell_rna,
                                      sc_membership = named_membership, mc_by_assay = 'ATAC',
                                      test_cutoff = 0.05)
  atac_mc_obj = multiome_obj$atac_mc_obj
  rna_mc_obj = multiome_obj$rna_mc_obj
  
  both_mc_obj = atac_mc_obj
  both_mc_obj[['RNA']] = rna_mc_obj[['RNA']]
  
  Annotation(both_mc_obj) = anno
  both_mc_obj <- RegionStats(both_mc_obj, genome = BSgenome.Hsapiens.UCSC.hg38)
  
  both_mc_obj_trust = subset(both_mc_obj, subset = testres == 'trustworthy')
  
  both_mc_obj <- LinkPeaks(
    object = both_mc_obj,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = gene_set,
    peak.slot = 'data',
    expression.slot = 'data',
    method = 'spearman',
    distance = 100000
  )
  
  both_mc_obj_trust <- LinkPeaks(
    object = both_mc_obj_trust,
    peak.assay = "peaks",
    expression.assay = "RNA",
    genes.use = gene_set,
    peak.slot = 'data',
    expression.slot = 'data',
    method = 'spearman',
    distance = 100000
  )
  
  GenePeakRes_Signac = list(trust = Links(both_mc_obj_trust),   # HCPs obtained using trustworthy metacells
                            all = Links(both_mc_obj))   # HCPs obtained using all metacells
  seacells_GenePeak_Signac[[tgamma]] = GenePeakRes_Signac
  
}

save(seacells_GenePeak_Signac, file = paste0(data_path, 'seacells_GenePeak_Signac_test95.RData'))   # intermediate results saved





###### Figure 1f & Supplementary Figure 5-6

load(paste0(data_path, 'seacells_GenePeak_Signac_test95.RData'))

tgamma = '90'
trust_ranges = seacells_GenePeak_Signac[[tgamma]]$trust
all_ranges = seacells_GenePeak_Signac[[tgamma]]$all
trust_ranges = subsetByOverlaps(trust_ranges, all_ranges, type = 'equal')
all_ranges = subsetByOverlaps(all_ranges, trust_ranges, type = 'equal')
trust_ranges
all_ranges
bbdat = data.frame(score = c(trust_ranges$score,
                             all_ranges$score),
                   pval = c(trust_ranges$pvalue,
                            all_ranges$pvalue),
                   gene = c(trust_ranges$gene,
                            all_ranges$gene),
                   group = c(rep('trust', length(trust_ranges)),
                             rep('all', length(all_ranges)))
)

bbdat_select = bbdat[bbdat$gene %in% gene_set,]
bbdat_select$gene = factor(bbdat_select$gene, levels = gene_set)
bbdat_select$score = abs(bbdat_select$score)
ggplot(bbdat_select, aes(x=gene, y=score)) +
  geom_boxplot(aes(fill = group)) +
  scale_fill_discrete(name = '') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))   # Figure 1f, top left

#
gene = 'GATA2'
splot_peakgene_corr(atac_mc_obj, rna_mc_obj, 
                    gene = gene,
                    corr_ranges = trust_ranges)  # Figure 1f, bottom left 1

gene = 'TAL1'
splot_peakgene_corr(atac_mc_obj, rna_mc_obj, 
                    gene = gene,
                    corr_ranges = trust_ranges)  # Figure 1f, bottom left 2


##
obj_singlecell_atac_rep1 = obj_singlecell_atac[, grep('rep1', colnames(obj_singlecell_atac))]
cell_names = colnames(obj_singlecell_atac_rep1)
cell_names = sapply(strsplit(cell_names, split = '#'), function(x) x[2])
obj_singlecell_atac_rep1 = RenameCells(obj_singlecell_atac_rep1, new.names = cell_names)

fragments_rep1 = CreateFragmentObject(path = paste0(data_path, 'BM_CD34_Rep1_atac_fragments.tsv.gz'),
                                      cells = colnames(obj_singlecell_atac_rep1))
head(fragments_rep1)
Fragments(obj_singlecell_atac_rep1) = fragments_rep1

#
both_sc_rep1 = obj_singlecell_atac_rep1

cell_names = colnames(both_sc_rep1)
cell_names = paste0('cd34_multiome_rep1#', cell_names)
both_sc_rep1 = RenameCells(both_sc_rep1, new.names = cell_names)

obj_singlecell_rna_rep1 = obj_singlecell_rna[, grep('rep1', colnames(obj_singlecell_rna))]
both_sc_rep1[['RNA']] = obj_singlecell_rna_rep1[['RNA']]

Annotation(both_sc_rep1) = anno
both_sc_rep1 <- RegionStats(both_sc_rep1, genome = BSgenome.Hsapiens.UCSC.hg38)  # Compute GC content... for Links computation

#
obj_singlecell_atac_rep2 = obj_singlecell_atac[, grep('rep2', colnames(obj_singlecell_atac))]
cell_names = colnames(obj_singlecell_atac_rep2)
cell_names = sapply(strsplit(cell_names, split = '#'), function(x) x[2])
obj_singlecell_atac_rep2 = RenameCells(obj_singlecell_atac_rep2, new.names = cell_names)

fragments_rep2 = CreateFragmentObject(path = paste0(data_path, 'BM_CD34_Rep2_atac_fragments.tsv.gz'),
                                      cells = colnames(obj_singlecell_atac_rep2))
head(fragments_rep2)
Fragments(obj_singlecell_atac_rep2) = fragments_rep2

#
cell_names = colnames(obj_singlecell_atac_rep1)
cell_names = paste0('cd34_multiome_rep1#', cell_names)
obj_singlecell_atac_rep1 = RenameCells(obj_singlecell_atac_rep1, new.names = cell_names)

cell_names = colnames(obj_singlecell_atac_rep2)
cell_names = paste0('cd34_multiome_rep2#', cell_names)
obj_singlecell_atac_rep2 = RenameCells(obj_singlecell_atac_rep2, new.names = cell_names)

obj_singlecell_atac_combined = merge(x = obj_singlecell_atac_rep1, 
                                     y = list(obj_singlecell_atac_rep2))
all(colnames(obj_singlecell_atac_combined) == colnames(obj_singlecell_atac))

#
both_sc_combined = obj_singlecell_atac_combined

both_sc_combined[['RNA']] = obj_singlecell_rna[['RNA']]

Annotation(both_sc_combined) = anno
both_sc_combined <- RegionStats(both_sc_combined, genome = BSgenome.Hsapiens.UCSC.hg38)  # Compute GC content... (added in meta.features) for Links computation

#
gene = 'GATA2'
CoveragePlot(object = both_sc_combined,
             assay = 'peaks',
             region = gene,
             features = gene,
             expression.assay = 'RNA',
             links = F,
             extend.downstream = 100000,
             extend.upstream = 100000)

#
load(paste0(data_path, "sc_GenePeak_Signac.RData"))

Links(both_sc_combined) = sc_GenePeak_Signac[['1']]$all

link_sc = CoveragePlot(object = both_sc_combined,
             assay = 'peaks',
             region = gene,
             features = gene,
             expression.assay = 'RNA',
             group.by = 'celltype',
             links = T,
             extend.downstream = 100000,
             extend.upstream = 100000)

#
load(paste0(data_path, "seacells_GenePeak_Signac_test95.RData"))

Links(both_sc_combined) = seacells_GenePeak_Signac[['75']]$all 

link_all = CoveragePlot(object = both_sc_combined,
             assay = 'peaks',
             region = gene,
             features = gene,
             expression.assay = 'RNA',
             group.by = 'celltype',
             links = T,
             extend.downstream = 100000,
             extend.upstream = 100000)

Links(both_sc_combined) = seacells_GenePeak_Signac[['75']]$trust

link_trust = CoveragePlot(object = both_sc_combined,
             assay = 'peaks',
             region = gene,
             features = gene,
             expression.assay = 'RNA',
             group.by = 'celltype',
             links = T,
             extend.downstream = 100000,
             extend.upstream = 100000)

link_sc + link_trust + link_all   # Figure 1f, right & Supplementary Figure 5b


##
Links(both_sc_combined) = seacells_GenePeak_Signac[['90']]$trust 

link_sg_mc = CoveragePlot(object = both_sc_combined,
                       assay = 'peaks',
                       region = gene,
                       features = gene,
                       expression.assay = 'RNA',
                       group.by = 'celltype',
                       links = T,
                       extend.downstream = 10000,
                       extend.upstream = 10000)

Links(both_sc_combined) = seacells_GenePeak_Signac[['5']]$all 

link_sg = CoveragePlot(object = both_sc_combined,
             assay = 'peaks',
             region = gene,
             features = gene,
             expression.assay = 'RNA',
             group.by = 'celltype',
             links = T,
             extend.downstream = 10000,
             extend.upstream = 10000)

Links(both_sc_combined) = sc_GenePeak_Signac[['1']]$all 

link_sg_sc = CoveragePlot(object = both_sc_combined,
                       assay = 'peaks',
                       region = gene,
                       features = gene,
                       expression.assay = 'RNA',
                       group.by = 'celltype',
                       links = T,
                       extend.downstream = 10000,
                       extend.upstream = 10000)

link_sg_mc + link_sg + link_sg_sc   # Supplementary Figure 6a


##
gene = 'TAL1'

Links(both_sc_combined) = seacells_GenePeak_Signac[['90']]$trust 

link_sg_mc1 = CoveragePlot(object = both_sc_combined,
                          assay = 'peaks',
                          region = gene,
                          features = gene,
                          expression.assay = 'RNA',
                          group.by = 'celltype',
                          links = T,
                          extend.downstream = 10000,
                          extend.upstream = 10000)

Links(both_sc_combined) = seacells_GenePeak_Signac[['5']]$all 

link_sg1 = CoveragePlot(object = both_sc_combined,
                       assay = 'peaks',
                       region = gene,
                       features = gene,
                       expression.assay = 'RNA',
                       group.by = 'celltype',
                       links = T,
                       extend.downstream = 10000,
                       extend.upstream = 10000)

Links(both_sc_combined) = sc_GenePeak_Signac[['1']]$all 

link_sg_sc1 = CoveragePlot(object = both_sc_combined,
                          assay = 'peaks',
                          region = gene,
                          features = gene,
                          expression.assay = 'RNA',
                          group.by = 'celltype',
                          links = T,
                          extend.downstream = 10000,
                          extend.upstream = 10000)

link_sg_mc1 + link_sg1 + link_sg_sc1   # Supplementary Figure 6b


##

sc_membership = cell_membership[['80']]
names(sc_membership) = rownames(cell_membership)

pp = mcRigor_projection(obj_singlecell = obj_singlecell_atac, sc_membership = sc_membership,
                        color_field = 'celltype',  axis_lab = F,
                        dub_mc_test.label = T, test_stats = TabMC, test_cutoff = 0.05)
pp   # Supplementary Figure 5a


