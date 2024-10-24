library(Seurat)
library(CellChat)
library(patchwork)
library(doSNOW)

options(stringsAsFactors = FALSE)
options(future.globals.maxSize = 2048 * 1024 * 1024)
future::plan("multisession", workers = 3)

seurat.obj = readRDS(seurat.file)
seurat.obj@meta.data[["CellName.CellChat"]] = seurat.obj$UCRSI.SVM
target.patients = unique(seurat.obj$Sample) # 患者或者供体
print(target.patients)
# seurat.obj$CellName.CellChat[seurat.obj$CNV == "Malignant"] = paste0(seurat.obj$CellName.CellChat[seurat.obj$CNV == "Malignant"], "(Malignant)")
# CellChatDB <- CellChatDB.human

cl <- makeCluster(5)
clusterEvalQ(cl, library(CellChat))
registerDoSNOW(cl)

pb <- txtProgressBar(max = length(target.patients), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

result_df <- foreach(patient = target.patients, .options.snow = opts, .combine='c') %dopar% {
  print(patient)
  cells = seurat.obj$Sample == patient & !(seurat.obj$CellName.CellChat %in% c("Unknow", "Unknown", "Doublet"))  # 患者或者供体
  meta <- data.frame(labels = seurat.obj$CellName.CellChat[cells], row.names = colnames(seurat.obj)[cells])
  cellchat <- createCellChat(object = seurat.obj@assays$RNA@data[, cells], meta = meta, group.by = "labels")
  cellchat@DB <- CellChatDB.human
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  results = list()
  results[[patient]] = cellchat
  # results[[patient]] = subsetCommunication(cellchat)
  return(results)
}
close(pb)
stopCluster(cl)
saveRDS(result_df, file = net.save.path)
