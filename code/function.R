library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(harmony)
library(ggplot2)
library(MLmetrics)
library(mltools)
library(xlsx)
library(scales)
library(jsonlite)

neck.cancer.markers = list(
  "Epithelial"=c("EPCAM", "KRT19", "KRT18", "KRT5", "KRT15"),
  "Endothelial"=c("CDH5", "CLDN5", "RAMP2"),
  "Fibroblast"=c("C1R", "COL1A2", "COL1A1", "COL3A1", "DCN"),
  "CD4"=c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "TRAC"),
  "CD8"=c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1"),
  "NK"=c("GNLY", "NKG7", "XCL1", "XCL2", "KLRC1", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"),
  "Tregs"=c("FOXP3", "IL2RA", "TNFRSF4"),
  "B"=c("CD79A", "CD79B", "CD19", "MS4A1"),
  "DC"=c("C1orf54", "LGALS2", "CD40", "CD80", "CD83", "CCR7"),
  "Mast"=c("KIT", "MS4A2", "PTGS1", "RGS13"),
  "Mφ"=c("FCGR1A", "CD163", "CD68", "FCGR2A", "CSF1R")
)

lung.cancer.markers = list(
  "AT1"=c("AGER", "CLIC5", "PDPN"),
  "ATII"=c("LPCAT1", "NAPSA", "PGC", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SLC34A2"),
  "Basal"=c("KRT17", "KRT5", "KRT6A"),
  "Cilia"=c("AKAP14", "ALDH3B1", "ANKRD66", "C11orf88", "C11orf97", "DNAI1"),
  "Club"=c("PIGR", "SCGB1A1", "SCGB3A1"),
  "EC"=c("CDH5", "CLDN5", "RAMP2"),
  "Fib"=c("C1R", "COL1A2", "COL1A1", "COL3A1", "DCN"),
  "NE"=c("AZGP1", "CPE", "TUBB2B"),
  "B"=c("CD19", "CD79A", "MS4A1"),
  "CD4"=c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "TRAC"),
  "CD8"=c("CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1"),
  "NK"=c("GNLY", "NKG7", "XCL1", "XCL2", "KLRC1", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"),
  "DC"=c("C1orf54", "LGALS2"),
  "Gran"=c("CD300E", "CXCL8", "EREG", "S100A12"),
  "Mast"=c("KIT", "MS4A2", "PTGS1", "RGS13"),
  "Mo"=c("FCGR1A", "CD163", "CD68", "FCGR2A", "CSF1R"),
  "Tregs"=c("FOXP3", "IL2RA", "TNFRSF4")
)

step.1 = function(seurat.obj, regress.out = FALSE) {
  seurat.obj <- NormalizeData(seurat.obj)
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = 2000)
  if(regress.out) {
    seurat.obj <- ScaleData(seurat.obj, vars.to.regress = "percent.mt")
  } else {
    seurat.obj <- ScaleData(seurat.obj)
  }
  seurat.obj <- RunPCA(seurat.obj, npcs = 30)
  seurat.obj@assays[["RNA"]]@scale.data = matrix(0)
  print(ElbowPlot(seurat.obj, ndims = 30))
  return(seurat.obj)
}

step.2 = function(seurat.obj, dims, meta.name = "Patients") {
  reduction = "pca"
  seurat.obj <- RunUMAP(seurat.obj, dims = 1:dims, verbose = T, reduction=reduction)
  print(DimPlot(seurat.obj, group.by = meta.name, label = T))
  run.harmony = readline(prompt="患者间是否存在明显的批次效应?(y/n)")
  if(!(run.harmony %in% c("y", "n"))) {
    return(seurat.obj)
  }
  if(run.harmony == "y") {
    seurat.obj <- RunHarmony(seurat.obj, meta.name, dims = 1:dims, project.dim = F)
    reduction = "harmony"
    seurat.obj <- RunUMAP(seurat.obj, dims = 1:dims, verbose = T, reduction=reduction)
  }
  seurat.obj <- FindNeighbors(seurat.obj, dims = 1:dims, reduction=reduction)
  seurat.obj <- FindClusters(seurat.obj, resolution = 1)
  return(seurat.obj)
}

step.3 = function(seurat.obj, markers) {
  cell.markers = list()
  markers.names = names(markers)
  for(i in 1:length(markers)) {
    cell.markers[[paste0(names(markers)[i], "(", i, ")")]] = markers[[names(markers)[i]]]
  }
  
  idx = "-1: Doublets\n0: Unknow\n"
  for(i in 1:length(cell.markers)) {
    idx = paste0(idx, i, ": ", markers.names[i], "\n")
  }
  assign.results = list()
  for(i in 0:(length(unique(seurat.obj$seurat_clusters)) - 1)) {
    print(DotPlot(seurat.obj, features = cell.markers, idents = i) + theme(legend.position = "none"))
    cat(idx)
    assign = readline(prompt=paste0("请输入细胞类型: "))
    if(is.na(as.integer(assign))) {
      return(seurat.obj)
    }
    if(assign == "0") {
     
    } else if(assign == "-1") {
      assign.results[[as.character(i)]] = "Doublet"
    } else {
      assign.results[[as.character(i)]] = markers.names[as.integer(assign)]
    }
    
  }
  seurat.obj@meta.data[["CellName"]] = as.character(seurat.obj@meta.data[["seurat_clusters"]])
  for(i in unique(seurat.obj@meta.data[["CellName"]])) {
    seurat.obj$CellName[seurat.obj$CellName == i] = assign.results[[i]]
  }
  return(seurat.obj)
}

step.4 = function(seurat.obj, datasets.path, GSE.id) {
  immnue.celltypes = c("CD4", "CD8", "NK", "Tregs", "B", "DC", "Mast", "Mφ")
  cells.selected = c()
  for(celltype in unique(seurat.obj$CellName)) {
    if(celltype %in% immnue.celltypes & sum(seurat.obj$CellName == celltype) > 1600) {
      cells.selected = c(cells.selected, sample(colnames(seurat.obj)[seurat.obj$CellName == celltype], 1600))
    } else {
      cells.selected = c(cells.selected, colnames(seurat.obj)[seurat.obj$CellName == celltype])
    }
  }
  seurat.obj[["CellID"]] = colnames(seurat.obj)
  seurat.obj.subset = subset(seurat.obj, subset = CellID %in% cells.selected)
  SaveH5Seurat(seurat.obj.subset, filename = paste0(datasets.path, GSE.id, ".seurat.h5Seurat"))
  Convert(paste0(datasets.path, GSE.id, ".seurat.h5Seurat"), dest = "h5ad")
  file.remove(paste0(datasets.path, GSE.id, ".seurat.h5Seurat"))
}

step.5 = function(seurat.obj, GSE.id) {
  results = read.csv(paste0("D:/Experiment/硕士毕业论文实验/infercnvpy/", GSE.id, ".csv"), header = T, row.names = 1)
  seurat.obj[["CNV"]] = rep("Normal", ncol(seurat.obj))
  seurat.obj$CNV[colnames(seurat.obj) %in% results$Malignant_Cells] = "Malignant"
  return(seurat.obj)
}

TNM = function(N.stage, patients) {
  for(patient in names(patients)) {
    if(patients[[patient]][2] %in% c("N0")) {
      N.stage$N = c(N.stage$N, patient)
    } else {
      N.stage$Y = c(N.stage$Y, patient)
    }
  }
  return(N.stage)
}

TNM.to.vector = function(patients.stage) {
  results = list(Y=c(), N=c())
  for(patient in names(patients.stage)) {
    if(patients.stage[[patient]][2] == "N0") {
      results$N = c(results$N, patient)
    } else {
      results$Y = c(results$Y, patient)
    }
  }
  return(results)
}

add.edge = function(net, source, target, strength, celltype1, celltype2, pathway, interaction.name) {
  if(sum(source %in% net$Source & target %in% net$Target) != 0) {
    idx = which(net$Source == source & net$Target == target)
    if("Weight" %in% colnames(net)) {
      net$Weight[idx] = as.numeric(net$Weight[idx]) + 1
    }
    net$Strength[idx] = max(as.numeric(net$Strength[idx]), as.numeric(strength))
  } else {
    net[nrow(net) + 1, ] = c(source, target, 1, strength, celltype1, celltype2, pathway, interaction.name)
  }
  return(net)
}

get.template.net = function(net.list, patients) {
  template.net = data.frame(matrix(nrow = 0, ncol = 8))
  colnames(template.net) = c("Source", "Target", "Weight", "Celltype1", "Celltype2", "Pathway", "InteractionName", "Source.Target")
  for(patient in patients) {
    print(patient)
    tmp = data.frame(
      Source=paste0(net.list[[patient]]$ligand, "(", net.list[[patient]]$source, ")"),
      Target=paste0(net.list[[patient]]$receptor, "(", net.list[[patient]]$target, ")"), 
      Weight=1, 
      Celltype1=as.character(net.list[[patient]]$source),
      Celltype2=as.character(net.list[[patient]]$target),
      Pathway=net.list[[patient]]$pathway_name,
      InteractionName=as.character(net.list[[patient]]$interaction_name)
    )
    tmp[["Source.Target"]] = paste0(tmp$Source, "->", tmp$Target)
    macting = match(tmp$Source.Target, template.net$Source.Target)
    template.net$Weight[macting[!is.na(macting)]] = template.net$Weight[macting[!is.na(macting)]] + 1
    template.net = rbind(template.net, tmp[is.na(macting), ])
  }
  print(nrow(template.net))
  template.net$Weight = as.numeric(template.net$Weight)
  return(template.net)
}

get.testing.net = function(net.list, patient) {
  return(data.frame(
    Source = paste0(net.list[[patient]]$ligand, "(", net.list[[patient]]$source, ")"),
    Target = paste0(net.list[[patient]]$receptor, "(", net.list[[patient]]$target, ")"),
    Celltype1 = as.character(net.list[[patient]]$source),
    Celltype2 = as.character(net.list[[patient]]$target)
  ))
}

network.to.json = function(network, path, name) {
  all.nodes = unique(c(network$Source, network$Target))
  results = list("edges" = list(), "features" = list())
  results$features = as.list(all.nodes)
  names(results$features) = 0:(length(results$features) - 1)
  for(i in 1:nrow(network)) {
    results$edges[[i]] = c(which(results$features == network$Source[i]) - 1, which(results$features == network$Target[i]) - 1)
  }
  write_json(results, pretty=F, auto_unbox=TRUE, path = paste0(path, name, ".json"))
}

radio.of.interaction2 = function(net.list, patients, source, target) {
  i = 0
  for(patient in patients) {
    loc = paste0(net.list[[patient]]$ligand, "(", net.list[[patient]]$source, ")") == source & paste0(net.list[[patient]]$receptor, "(", net.list[[patient]]$target, ")") == target
    if(sum(loc) > 0) {
      i = i + 1
    }
  }
  return(i / length(patients))
}

radio.of.interaction = function(net.list, patients, interaction.name) {
  i = 0
  for(patient in patients) {
    if(interaction.name %in% net.list[[patient]]$interaction_name) {
      i = i + 1
    }
  }
  return(i / length(patients))
}

radio.of.pathway = function(net.list, patients, pathway.name) {
  i = 0
  for(patient in patients) {
    if(pathway.name %in% net.list[[patient]]$pathway_name) {
      i = i + 1
    }
  }
  return(i / length(patients))
}
