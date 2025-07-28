library(Seurat)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(dplyr)
library(patchwork)
library(clusterProfiler)
library(TxDb.Sscrofa.UCSC.susScr11.refGene)
library(org.Ss.eg.db)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(scran)
library(ArchR)
library(tidyverse)
library(parallel)
library(pheatmap)
library(do)
library(RColorBrewer)
library(ggsci)
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))
#构建参考基因组
genomeAnnotation1 <- createGenomeAnnotation(genome = BSgenome.Sscrofa.UCSC.susScr11)
geneAnnotation1 <- createGeneAnnotation(TxDb = TxDb.Sscrofa.UCSC.susScr11.refGene, OrgDb = org.Ss.eg.db)
addArchRChrPrefix(chrPrefix = FALSE)

#建立路径
setwd("/home/wangliyuan/scATAC/integ_all_inall")
integ_all <- loadArchRProject(path = "/home/wangliyuan/scATAC/integ_all_inall")

options(future.globals.maxSize = 3000 * 1024^2)
addArchRThreads(threads = 40)
pbmc_rna <- readRDS("/home/wangliyuan/scRNA_datas/scRNA_integ_all_202204/true_combine_merge/combined_All_afterName.rds")
pbmc <- readRDS("pbmc_ATAC.rds")
DefaultAssay(pbmc_rna) <- "RNA" 
transfer.anchors <- FindTransferAnchors(reference = pbmc_rna,query = pbmc, reduction = 'lsiproject')

predicted.labels <- TransferData(anchorset = transfer.anchors,refdata = pbmc_rna$cell_type,weight.reduction = pbmc[['lsi']], dims = 1:30)
pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)
p1 <- DimPlot(object = pbmc, label = TRUE, group.by = 'Clusters', repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("pbmc_ATCA_Clusters.pdf", plot = p1, width = 10, height = 8) 
plot5 <- DimPlot(object = pbmc, label = TRUE, group.by='Sample', repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("pbmc_ATCA__sample.tiff", plot = plot5, width = 12, height = 8)
p2 <- DimPlot(object = pbmc, label = TRUE, group.by = 'predicted.id', repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("pbmc_ATCA_predicted.id.pdf", plot = p2, width = 15, height = 8) 
pbmc$Tech <- "ATAC"
pbmc_rna$Tech <- "RNA"
saveRDS(pbmc,file="pbmc_ATAC_lsiproject.rds")


#scATAC_based
dir.create("scATAC_based")
refdata <- GetAssayData(pbmc_rna, assay = "RNA", slot = "data")[genesUse, ]
imputation <- TransferData(anchorset = transfer.anchors, refdata = refdata,	weight.reduction = pbmc[["lsi"]], dim = 1:30)
pbmc[["RNA"]] <- imputation
DefaultAssay(pbmc) <- 'RNA'
coembed <- merge(x = pbmc_rna, y = pbmc)
coembed <- ScaleData(coembed, features = genesUse, do.scale = FALSE)
coembed <- RunPCA(coembed, features = genesUse, verbose = FALSE)
coembed <- RunUMAP(coembed, dims = 1:30)
coembed <- RunTSNE(coembed, dims = 1:30)
coembed$Merged_cluster <- ifelse(!is.na(coembed$cell_type), coembed$cell_type, coembed$predicted.id)
#coembed better than integrated
#coembed$Merged_cluster <- Replace(data = coembed$Merged_cluster, from = "Monocyte",to = "Kupffer cells")
coembed$Merged_cluster <- Replace(data = coembed$Merged_cluster, from = "NKs and γδ cell",to = "NK cells")
#coembed$Merged_cluster <- Replace(data = coembed$Merged_cluster, from = c("[^[:alnum:]///' ]", " TUBB6"),to = "")
coembed$Merged_cluster <- Replace(data = coembed$Merged_cluster, from = "Naïve T cell",to = "Naive T cell")

p1 <- DimPlot(coembed, group.by = "Tech", reduction = "umap", repel = TRUE, cols = colpalettes, raster=FALSE) 
ggsave("coembed_Tech_umap.pdf", plot = p1, width = 10, height = 8) 
p2 <- DimPlot(coembed, group.by = "Merged_cluster", reduction = "umap", label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("umap_merged_cluster.pdf", plot = p2, width = 15, height = 8) 
p3 <- DimPlot(coembed, group.by = "Tech", reduction = "tsne") 
ggsave("coembed_Tech_tsne.pdf", plot = p3, width = 10, height = 8) 
p4 <- DimPlot(coembed, group.by = "Merged_cluster", reduction = "tsne", label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("tsne_merged_cluster.pdf", plot = p4, width = 15, height = 8) 
saveRDS(coembed,file="coembed.rds")

#scRNA_based
DefaultAssay(pbmc) <- 'GeneScore'
reference.list <- c(pbmc_rna, pbmc)
names(reference.list) <- c("RNA", "ATAC")
rna_atac.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
rna_atac_integrated <- IntegrateData(anchorset = rna_atac.anchors, dims = 1:30)
rna_atac_integrated <- ScaleData(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- RunPCA(object = rna_atac_integrated, verbose = F)
rna_atac_integrated <- FindNeighbors(object = rna_atac_integrated, dims = 1:30)
rna_atac_integrated <- FindClusters(object = rna_atac_integrated, resolution = 2)
rna_atac_integrated <- RunUMAP(object = rna_atac_integrated, reduction = "pca", dims = 1:30)
rna_atac_integrated <- RunTSNE(object = rna_atac_integrated, reduction = "pca", dims = 1:30)
rna_atac_integrated$Merged_cluster <- ifelse(!is.na(rna_atac_integrated$cell_type), rna_atac_integrated$cell_type, rna_atac_integrated$predicted.id)

p1 <- DimPlot(rna_atac_integrated, group.by = "Tech", reduction = "umap")
ggsave("rna_atac_integrated_Tech_umap.pdf", plot = p1, width = 10, height = 8) 
p2 <- DimPlot(rna_atac_integrated, group.by = "Merged_cluster", reduction = "umap", label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("rna_atac_integrated_umap_merged_cluster.pdf", plot = p2, width = 15, height = 8) 
p3 <- DimPlot(rna_atac_integrated, group.by = "Tech", reduction = "tsne")
ggsave("rna_atac_integrated_Tech_tsne.pdf", plot = p3, width = 10, height = 8)
p4 <- DimPlot(rna_atac_integrated, group.by = "Merged_cluster", reduction = "tsne", label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("rna_atac_integrated_tsne_merged_cluster.pdf", plot = p4, width = 15, height = 8) 
saveRDS(rna_atac_integrated,file="rna_atac_integrated.rds")


ATAC_coembed<-subset(coembed, Tech %in% "ATAC")
#ATAC_coembed$Sample <- factor(ATAC_coembed$Sample, levels = c("liver_E65", "liver_D1", "liver_Y1", "liver_Y3"))
#ATAC_coembed$Merged_cluster <- factor(ATAC_coembed$Merged_cluster, levels = c("Endothelia cell", "Kupffer cells", "Macrophage", "Hepatic stellate cells", "Cholangiocyte", "Myofibroblast", "Hepatocytes", "Erythroblasts", "NK cells", "B cells", "Neutrophils", "Naive T cell", "Hematopoietic stem cell", "Plasma cell", "Mast cell"))
plote <- DimPlot(ATAC_coembed, reduction = "umap", group.by='Merged_cluster', label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("ATAC_coembed_umap.tiff", plot = plote, width = 18, height = 8) 
plot5 <- DimPlot(ATAC_coembed, reduction = "umap", group.by='Sample', repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("ATAC_coembed_umap_sample.tiff", plot = plot5, width = 12, height = 8)
plote <- DimPlot(ATAC_coembed, reduction = "tsne", group.by='Merged_cluster', label = TRUE, repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("ATAC_coembed_tsne.tiff", plot = plote, width = 18, height = 8) 
plot5 <- DimPlot(ATAC_coembed, reduction = "tsne", group.by='Sample', repel = TRUE, cols = colpalettes, raster=FALSE)
ggsave("ATAC_coembed_tsne_sample.tiff", plot = plot5, width = 12, height = 8)

Idents(ATAC_coembed) <- ATAC_coembed$Merged_cluster
markers.to.plot <- c("CD79B", "CD19", "ZNF70", "LUM", "MGP", "MFAP5", "CD5", "CD3E", "CD3D", "PECAM1", "VWF", "CDH5", "CALCR", "P2RY12", "BAG3", "PROCR", "SCARA5", "CCL26", "FHL5", "KCNJ8", "TBX2", "CD8A", "GNLY", "CD244", "FCN1", "DUSP2", "PLBD1", "STMN2", "HES6", "NEUROD2", "CD209", "LYVE1", "FGFR4", "DLK1", "MSTN", "MYOD1", "LRP2", "NOX4", "UPP2", "IRS1", "SPOCK3", "ICA1", "UMOD", "TMEM52B", "SCNN1A", "MYL1", "PENK", "PCSK5", "EREG", "LY6D", "SHH", "TDRP", "MFSD2A", "EPHA4", "HBZ", "HBM", "HEMGN", "GPR84", "NRG1", "S100A9", "MARCO", "PRPS2", "P2RY13", "ITGA9", "CCDC32", "DLL1", "CLDN18", "PGC", "VSIG2", "CRISP3", "PGLYRP1", "CD177", "PNLIPRP1", "CLPS", "CPA1", "AQP5", "PHEROC", "SCGB1D1", "SOX2", "PCDH8", "LHFPL3", "NPPC", "SGCG", "CHRDL2", "BMP3", "NINL", "KCTD3", "MMRN1", "RELN", "REEP1", "SCGB3A1", "CLIC6", "PIGR", "XCR1", "NAPSA", "WARS", "MYOC", "AOX1", "PDGFRA", "RBP7", "SPRY4", "TIMP4", "MTUS2", "CHD5", "MYOM1", "MTMR3", "PRR5", "NDRG4", "CLEC4G", "NPR3", "WNT2", "ACAN", "CNMD", "CHAD", "FABP6", "GUCA2A", "LGALS2", "GP1BB", "MRPS23", "VARS", "CELA2A", "GP2", "DNASE1", "MYH7", "MYL3", "MYOZ2", "MYL4", "TNNT2", "ANKRD1", "TNFRSF17", "LYPD2", "TNFRSF13B", "JCHAIN", "TXNDC11", "BHLHA15", "MPZ", "SCN7A", "NGFR", "RHEX", "USP3", "BLK", "ASAP3", "SOCS5", "SFRP1", "MYH4", "MYBPC2", "ACTN3", "BACH2", "FAM129C", "LMTK2", "ANO3", "PAPPA", "AMIGO2", "ADIPOQ", "PLIN1", "PPP1R3C", "CLEC1A", "GCH1", "TNFRSF11A", "GPR75", "SLCO1C1", "RRAGD", "CYP7A1", "DUOX2", "ACKR4", "BMP6", "TNMD", "CDHR1", "CHGA", "CHGB", "SCG2", "TNFAIP3", "TRAF1", "SFTPA1", "CD8B", "S100Z", "CLEC12A", "HDC", "GRAP2", "HPGDS", "MKRN1", "VNN2", "OSCAR", "PRTN3", "ELANE", "MPO", "FGB", "PLG", "FGG", "PCED1B", "ZC3HAV1", "RUNX3")
plotE11 <- DotPlot(ATAC_coembed, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
ggsave("ATAC_coembed_MarkergenesDots.pdf", plot = plotE11, width = 40, height = 10, limitsize = FALSE)

cM <- confusionMatrix(paste0(ATAC_coembed$Clusters), paste0(ATAC_coembed$Merged_cluster))
write.csv(cM, file = "ATAC_coembed_Clusters_Merged_cluster.csv")
cellType <- ATAC_coembed$Merged_cluster
cellType1 <- as.matrix(cellType)
cellnames <- colnames(ATAC_coembed)
integ_all$cellType1 <- rownames(integ_all)
integ_all$cellType1 <- mapLabels(integ_all$cellType1, newLabels = cellType1, oldLabels = cellnames)
integ_all <- saveArchRProject(ArchRProj = integ_all)


markersGS <- getMarkerFeatures(ArchRProj = integ_all, useMatrix = "GeneScoreMatrix",  groupBy = "cellType1", bias = c("TSSEnrichment", "log10(nFrags)"), testMethod = "wilcoxon", logFile = createLogFile("getMarkerFeatures") )
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
write.csv(markerList, file = "markergene_cellType1.csv")
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "integ_all_GeneScores-Marker-Heatmap_cellType1", width = 8, height = 6, ArchRProj = integ_all, addDOC = FALSE)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.5 & Log2FC >= 1")
write.csv(markerList, file = "widemarkergene_cellType1.csv")
heatmapGS2 <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 1", 
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS2, name = "integ_all_GeneScores-Marker-Heatmap-wide_cellType1", width = 8, height = 6, ArchRProj = integ_all, addDOC = FALSE)

