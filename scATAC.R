library(ArchR)
library(TxDb.Sscrofa.UCSC.susScr11.refGene)
library(org.Ss.eg.db)
library(BSgenome.Sscrofa.UCSC.susScr11)
library(scran)

# 1. 初始化设置 ------------------------------------------------------------
# 加载参考基因组
integ_liver <- loadArchRProject(path = "/vol1/1Digest_liver/1Digest_liver_ATAC")
genomeAnnotation1 <- getGenomeAnnotation(ArchRProj = integ_liver)
geneAnnotation1 <- createGeneAnnotation(TxDb = TxDb.Sscrofa.UCSC.susScr11.refGene, 
                                      OrgDb = org.Ss.eg.db)
addArchRChrPrefix(chrPrefix = FALSE)

# 设置工作目录
setwd("/vol1/integ_all_inall")

# 2. 数据输入 --------------------------------------------------------------
sample_info <- list(
  E65 = c("Heart", "Liver", "Spleen", "Lung", "Kidney", "Pancreas", "PBMC", 
          "Cerebrum", "Cerebellum", "Dorsal_muscle", "Leg_muscle", "Adipose"),
  D1 = c("Heart", "Liver", "Spleen", "Lung", "Kidney", "Pancreas", "Bonemarrow", 
         "PBMC", "Cerebrum", "Cerebellum", "Dorsal_muscle", "Leg_muscle", 
         "Adipose", "Ovary", "Trachea", "Bladder"),
  Y1 = c("Heart", "Liver", "Spleen", "Lung", "Kidney", "Pancreas", "Bonemarrow",
         "PBMC", "Cerebrum", "Cerebellum", "Dorsal_muscle", "Leg_muscle",
         "Adipose", "Ovary", "Trachea", "Bladder"),
  Y3 = c("Heart", "Liver", "Spleen", "Lung", "Kidney", "Pancreas", "Bonemarrow",
         "PBMC", "Cerebrum", "Cerebellum", "Dorsal_muscle", "Leg_muscle",
         "Adipose", "Ovary", "Trachea", "Bladder")
)

sampleNames <- unlist(lapply(names(sample_info), function(time) {
  paste(sample_info[[time]], time, sep = "_")
}))

# 3. 创建Arrow文件 --------------------------------------------------------
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles1,
  sampleNames = sampleNames,
  geneAnnotation = geneAnnotation1,
  genomeAnnotation = genomeAnnotation1,
  filterTSS = 4,
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  force = TRUE,
  QCDir = file.path(getwd(), "QualityControl")
)

# 4. 创建ArchR项目 --------------------------------------------------------
integ_all <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = getwd(),
  geneAnnotation = geneAnnotation1,
  genomeAnnotation = genomeAnnotation1,
  copyArrows = FALSE
)

# 5. 质量控制可视化 -------------------------------------------------------
# 定义绘图函数
plot_qc <- function(proj, plot_name) {
  # TSS vs Fragments
  df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment"))
  p1 <- ggPoint(
    x = df[,1], y = df[,2], 
    colorDensity = TRUE, continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments", ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + 
      geom_vline(xintercept = 3, lty = "dashed")
  
  # Sample statistics
  p2 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", 
                  name = "TSSEnrichment", plotAs = "ridges")
  p3 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", 
                  name = "TSSEnrichment", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
  p4 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", 
                  name = "log10(nFrags)", plotAs = "ridges")
  p5 <- plotGroups(proj, groupBy = "Sample", colorBy = "cellColData", 
                  name = "log10(nFrags)", plotAs = "violin", alpha = 0.4, addBoxPlot = TRUE)
  
  # Fragment sizes and TSS profile
  p6 <- plotFragmentSizes(ArchRProj = proj)
  p7 <- plotTSSEnrichment(ArchRProj = proj)
  
  # 保存所有图形
  plotPDF(p1, name = paste0(plot_name, "_TSS-vs-Frags.pdf"), 
         ArchRProj = proj, addDOC = FALSE)
  plotPDF(p2, p3, p4, p5, name = paste0(plot_name, "_QC-Sample-Statistics.pdf"), 
         ArchRProj = proj, addDOC = FALSE, width = 4, height = 4)
  plotPDF(p6, p7, name = paste0(plot_name, "_QC-Sample-FragSizes-TSSProfile.pdf"), 
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

plot_qc(integ_all, "integ_all")

# 6. 降维分析 ------------------------------------------------------------
integ_all <- addIterativeLSI(
  ArchRProj = integ_all,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(resolution = c(0.1, 0.2, 0.4), 
                      sampleCells = 10000, 
                      n.start = 10), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

# Harmony校正
integ_all <- addHarmony(
  ArchRProj = integ_all,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample"
)

# 7. 聚类分析 ------------------------------------------------------------
integ_all <- addClusters(
  input = integ_all,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.9,
  force = TRUE
)

# 保存聚类结果
write.csv(table(integ_all$Clusters), 'Cluster.csv', row.names = FALSE)

# 聚类热图
cM <- confusionMatrix(paste0(integ_all$Clusters), paste0(integ_all$Sample))
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteRed"), 
  border_color = "black"
)
plotPDF(p, name = "heatmap of clusters.pdf", ArchRProj = integ_all, 
       addDOC = FALSE, width = 5, height = 5)

# 8. 可视化 --------------------------------------------------------------
# 定义可视化函数
plot_embeddings <- function(proj, reduc, plot_name) {
  # UMAP/TSNE by Sample and Clusters
  p1 <- plotEmbedding(proj, colorBy = "cellColData", name = "Sample", embedding = reduc)
  p2 <- plotEmbedding(proj, colorBy = "cellColData", name = "Clusters", embedding = reduc)
  
  # 对齐并保存图形
  ggAlignPlots(p1, p2, type = "h")
  plotPDF(p1, p2, name = paste0(plot_name, "_Plot-", reduc, "-Sample-Clusters.pdf"), 
         ArchRProj = proj, addDOC = FALSE, width = 5, height = 5)
}

# 添加UMAP和tSNE
integ_all <- addUMAP(integ_all, reducedDims = "IterativeLSI", name = "UMAP", 
                    nNeighbors = 30, minDist = 0.5, metric = "cosine")
integ_all <- addTSNE(integ_all, reducedDims = "IterativeLSI", name = "TSNE", 
                    perplexity = 30)
integ_all <- addUMAP(integ_all, reducedDims = "Harmony", name = "UMAPHarmony", 
                    nNeighbors = 30, minDist = 0.5, metric = "cosine")
integ_all <- addTSNE(integ_all, reducedDims = "Harmony", name = "TSNEHarmony", 
                    perplexity = 30)

# 绘制所有嵌入图
plot_embeddings(integ_all, "UMAP", "integ_all")
plot_embeddings(integ_all, "TSNE", "integ_all")
plot_embeddings(integ_all, "UMAPHarmony", "integ_all")
plot_embeddings(integ_all, "TSNEHarmony", "integ_all")

# 9. 标记基因分析 --------------------------------------------------------
markersGS <- getMarkerFeatures(
  ArchRProj = integ_all, 
  useMatrix = "GeneScoreMatrix",  
  groupBy = "Clusters", 
  bias = c("TSSEnrichment", "log10(nFrags)"), 
  testMethod = "wilcoxon"
)

# 保存标记基因结果
save_markers <- function(markers, cut_off, file_name) {
  markerList <- getMarkers(markers, cutOff = cut_off)
  write.csv(markerList, file = file_name)
  
  heatmap <- markerHeatmap(
    seMarker = markers, 
    cutOff = cut_off, 
    transpose = TRUE
  )
  ComplexHeatmap::draw(heatmap, heatmap_legend_side = "bot", annotation_legend_side = "bot")
  plotPDF(heatmap, name = gsub("\\.csv", "-Heatmap", basename(file_name)), 
         width = 8, height = 6, ArchRProj = integ_all, addDOC = FALSE)
}

save_markers(markersGS, "FDR <= 0.01 & Log2FC >= 1.25", "markergeneCluster.csv")
save_markers(markersGS, "FDR <= 0.1 & Log2FC >= 1", "widemarkergeneCluster.csv")

# 10. 保存项目 -----------------------------------------------------------
integ_all <- saveArchRProject(
  ArchRProj = integ_all,
  outputDirectory = getwd(),
  load = TRUE
)