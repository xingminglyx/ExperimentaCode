library(Seurat)
library(reticulate)
library(aricode)
library(ggpubr)
library(aricode)
library(Matrix)
use_condaenv("tensorflow")
sklearn.svm <- import("sklearn.svm")
sklearn.calibration <- import("sklearn.calibration")
np <- import("numpy")
joblib <- import("joblib")
sklearn.metrics <- import("sklearn.metrics")
pyamg = import("pyamg")
sklearn.cluster = import("sklearn.cluster")

Zhang.seurat = readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung cancer/Zhang et al.2022/Zhang.seurat.rds")
Leader.seurat = readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung cancer/Leader et al.2021/Leader.seurat.rds")
GSE148071.seurat <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/GSE148071.seurat.rds")
E_MTAB_6149.seurat <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/E_MTAB_6149.seurat.rds")
GSE127465.seurat <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/GSE127465.seurat.rds")
GSE117570.seurat <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/GSE117570.seurat.rds")
Maynard.seurat <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/Maynard.seurat.rds")

share.genes = intersect(rownames(Zhang.seurat), rownames(Leader.seurat))
share.genes = intersect(share.genes, rownames(GSE148071.seurat))
share.genes = intersect(share.genes, rownames(E_MTAB_6149.seurat))
share.genes = intersect(share.genes, rownames(GSE127465.seurat))
share.genes = intersect(share.genes, rownames(GSE117570.seurat))
share.genes = intersect(share.genes, rownames(Maynard.seurat))
saveRDS(share.genes, "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/share_genes.rds")

share.genes = readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/share_genes.rds")

data = Zhang.seurat@assays$RNA@data
data = as.data.frame(t(as.data.frame(data[share.genes, ])))
data[["CellName"]] = as.character(Zhang.seurat$CellName)
saveRDS(data, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/Zhang_data.rds")

target = list(
  "AT1" = 1,
  "ATII" = 2,
  "Basal" = 3,
  "Cilia" = 4,
  "Club" = 5,
  "EC" = 6,
  "Fib" = 7,
  "NE" = 8,
  "B" = 9,
  "CD4" = 10,
  "CD8" = 11,
  "DC" = 12,
  "Gran" = 13,
  "Mast" = 14,
  "Mo" = 15,
  "NK" = 16,
  "Tregs" = 17
)
for(celltype in names(target)) {
  data[["CellName"]][data[["CellName"]] == celltype] = target[[celltype]]
}
data[["CellName"]] = as.integer(data[["CellName"]])


# 训练模型 ---------------------------------------------------------------------
data <- readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/Zhang_data.rds")
cellname = data[["CellName"]]
data = data[, setdiff(colnames(data), c("CellName"))]
svm.linear = sklearn.calibration$CalibratedClassifierCV(sklearn.svm$LinearSVC(), n_jobs=as.integer(1))
svm.linear$fit(data, cellname)
joblib$dump(svm.linear, "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/svm_linear.pkl")

# 预测类型 ---------------------------------------------------------------------
predict.celltype = function(svm.linear, seurat.obj.path, data.path, share.genes, threshold = 0.5) {
  target = list(
    "AT1" = 1,
    "ATII" = 2,
    "Basal" = 3,
    "Cilia" = 4,
    "Club" = 5,
    "EC" = 6,
    "Fib" = 7,
    "NE" = 8,
    "B" = 9,
    "CD4" = 10,
    "CD8" = 11,
    "DC" = 12,
    "Gran" = 13,
    "Mast" = 14,
    "Mo" = 15,
    "NK" = 16,
    "Tregs" = 17
  )
  print("读取seurat对象...")
  seurat.obj = readRDS(seurat.obj.path)
  print("读取表达矩阵...")
  if(!file.exists(data.path)) {
    data = as.data.frame(t(as.data.frame(seurat.obj@assays$RNA@data[share.genes, ])))
    saveRDS(data, file = data.path)
  } else {
    data = readRDS(data.path)
  }
  print("预测类型...")
  predicted = svm.linear$predict(data)
  for(celltype in names(target)) {
    predicted[predicted == target[[celltype]]] = celltype
  }
  prob = np$max(svm.linear$predict_proba(data), axis = as.integer(1))
  unlabeled = np$where(prob < threshold)[[1]]
  predicted[unlabeled] = 'Unknown'
  seurat.obj@meta.data[["SVM.assign"]] = predicted
  gc()
  return(seurat.obj)
}

predict.celltype.within.cluster = function(svm.linear, patient.obj, data, threshold = 0.5) {
  target = list(
    "AT1" = 1,
    "ATII" = 2,
    "Basal" = 3,
    "Cilia" = 4,
    "Club" = 5,
    "EC" = 6,
    "Fib" = 7,
    "NE" = 8,
    "B" = 9,
    "CD4" = 10,
    "CD8" = 11,
    "DC" = 12,
    "Gran" = 13,
    "Mast" = 14,
    "Mo" = 15,
    "NK" = 16,
    "Tregs" = 17
  )
  # data = as.data.frame(t(as.data.frame(patient.obj@assays$RNA@data[share.genes, ])))
  patient.obj@meta.data[["SVM.assign"]] = rep("", ncol(patient.obj))
  celltypes.assign = list()
  for(ident in 0:(max(patient.obj@meta.data[["UCRSI"]]))) {
    cells = patient.obj@meta.data[["UCRSI"]] == ident
    # predicted = svm.linear$predict(data[colnames(patient.obj)[cells], ])
    prob = svm.linear$predict_proba(data[colnames(patient.obj)[cells], ])
    predicted = svm.linear$classes_[np$argmax(prob, axis = as.integer(1)) + 1]
    
    unlabeled = np$where(np$max(prob, axis = as.integer(1)) < threshold)[[1]]
    predicted[unlabeled] = 'Unknown'
    
    pred_summary = summary(as.factor(predicted))
    pred_summary = pred_summary[order(pred_summary, decreasing = T)]
    patient.obj$SVM.assign[cells] = names(pred_summary)[1]
    celltypes.assign[[as.character(ident)]] = pred_summary[1:2]
  }
  # print(celltypes.assign)
  gc()
  return(patient.obj)
}

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung cancer/Leader et al.2021/Leader.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/Leader_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/GSE148071.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE148071_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/E_MTAB_6149.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/E_MTAB_6149_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/GSE127465.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE127465_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/GSE117570.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE117570_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/Maynard.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/Maynard_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE153935/GSE153935.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE153935_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/GSE131907.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE131907_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/GSE136831.Control.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE136831_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/PRJEB31843/PRJEB31843.seurat.rds" # 正常样本
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/PRJEB31843_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/UKIM-V/UKIM.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/UKIM_data.rds"

seurat.obj.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/GSE136831.Control.seurat.rds"
data.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/GSE136831_data.rds"




setwd("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/UCRSI")
source(paste0(getwd(), "/FeatureSelection.R"))
source(paste0(getwd(), "/UCRSI.R"))
n_neighbors = 10
min_Gene = 4
output_distance = FALSE
alpha_start = 0
alpha_end = 1
alpha_step = 0.1
nJobs = 10
svm.linear = joblib$load("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/svm_linear.pkl") 
share.genes = readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/SVM/share_genes.rds")

seurat.obj = readRDS(seurat.obj.path)
if(!file.exists(data.path)) {
  data = as.data.frame(t(as.data.frame(seurat.obj@assays$RNA@data[share.genes, ])))
  saveRDS(data, file = data.path)
} else {
  data = readRDS(data.path)
}

seurat.obj@meta.data[["UCRSI.SVM"]] = rep("", ncol(seurat.obj))
for(patient in unique(seurat.obj$Sample)) { #
  print(paste0(patient, ": ", sum(seurat.obj$Sample == patient))) #
  n_cluster = 17
  if(sum(seurat.obj$Sample == patient) < 200) {
    next
  }
  patient.obj = subset(seurat.obj, subset = Sample == patient) #
  patient.obj = NormalizeData(patient.obj)
  features = features_selection(patient.obj, nbin = max(10, n_cluster * 2))
  UCRSI.results = UCRSI(as.matrix(patient.obj@assays[["RNA"]]@data[features$feature.genes, ]),
                        minGene = min_Gene,
                        n_neighbors = n_neighbors,
                        outputDistance = output_distance, 
                        alpha_start = alpha_start, 
                        alpha_end = alpha_end, 
                        alpha_step = alpha_step, 
                        nJobs = nJobs,
                        infValue = 0x7fffffff)
  patient.obj@meta.data[["UCRSI"]] = sklearn.cluster$SpectralClustering(
                                        n_clusters=as.integer(n_cluster), 
                                        n_neighbors = as.integer(n_neighbors),
                                        affinity="precomputed",
                                        assign_labels='discretize', 
                                        n_jobs = as.integer(nJobs),
                                        eigen_solver = "amg",
                                        random_state=as.integer(1))$fit_predict(UCRSI.results$adjacency)
  
  patient.obj = predict.celltype.within.cluster(svm.linear, patient.obj, data, threshold = 0.32)
  seurat.obj$UCRSI.SVM[seurat.obj$Sample == patient] = patient.obj$SVM.assign #
  gc()
}


summary(as.factor(seurat.obj$UCRSI.SVM))
ggarrange(DimPlot(seurat.obj, group.by = "CellName", label = T, raster = T, raster.dpi = c(1024, 1024)), 
          DimPlot(seurat.obj, group.by = "UCRSI.SVM", label = T, raster = T, raster.dpi = c(1024, 1024)), 
          ncol = 2, nrow = 1)
DotPlot(seurat.obj, group.by = "UCRSI.SVM", features = lung.cancer.markers)
saveRDS(seurat.obj, seurat.obj.path)

