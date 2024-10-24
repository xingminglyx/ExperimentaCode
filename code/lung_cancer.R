library(Seurat)
library(jsonlite)
library(xlsx)
library(ggplot2)
library(igraph)

# markers = list(
#   "AT1"=c("AGER", "CLIC5", "PDPN"),
#   "ATII"=c("LPCAT1", "NAPSA", "PGC", "SFTPA1", "SFTPA2", "SFTPB", "SFTPC", "SLC34A2"),
#   "Basal"=c("KRT17", "KRT5", "KRT6A"),
#   "Cilia"=c("AKAP14", "ALDH3B1", "ANKRD66", "C11orf88", "C11orf97", "DNAI1"),
#   "Club"=c("PIGR", "SCGB1A1", "SCGB3A1"),
#   "EC"=c("CDH5", "CLDN5", "RAMP2"),
#   "Fib"=c("C1R", "COL1A2", "DCN"),
#   "NE"=c("AZGP1", "CPE", "TUBB2B"),
#   "B"=c("CD19", "CD79A", "MS4A1"),
#   "CD4"=c("CD2", "CD3D", "CD3E", "CD3G", "CD4", "TRAC"),
#   "CD8"=c("PTPRC,", "CD8A", "CD8B", "GZMA", "GZMB", "GZMH", "GZMK", "PRF1"),
#   "DC"=c("C1orf54", "LGALS2"),
#   "Gran"=c("CD300E", "CXCL8", "EREG", "S100A12"),
#   "Mast"=c("KIT", "MS4A2", "PTGS1", "RGS13"),
#   "Mφ"=c("CD68", "FCGR1A", "ITGAX"),
#   "NK"=c("GNLY", "NKG7", "XCL1", "XCL2", "KLRC1", "KLRD1", "KLRF1", "FCGR3A", "NCAM1"),
#   "Tregs"=c("FOXP3", "IL2RA", "TNFRSF4")
# )
# 
patients.celltypes.summary = function(seurat.obj.path, patient.meta.name, celltype.meta.name) {
  results = list()
  seurat.obj = readRDS(seurat.obj.path)
  for(patient in unique(seurat.obj@meta.data[[patient.meta.name]])) {
    results[[patient]] = unique(seurat.obj@meta.data[[celltype.meta.name]][seurat.obj@meta.data[[patient.meta.name]] == patient])
  }
  return(results)
}
patient.celltypes = list()
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Zhang et al.2022/Zhang.seurat.rds",
                                                                    "Patients",
                                                                    "CellName"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Leader et al.2021/Leader.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/GSE148071.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/E_MTAB_6149.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/GSE127465.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/GSE117570.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE153935/GSE153935.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/PRJEB31843/PRJEB31843.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/GSE136831.seurat.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
patient.celltypes = c(patient.celltypes, patients.celltypes.summary("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/GSE131907.rds",
                                                                    "Patients",
                                                                    "UCRSI.SVM"))
saveRDS(patient.celltypes, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/patient.celltypes.rds")
# # 模板网络
# net.list = c(
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Leader et al.2021/cellchat.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Zhang et al.2022/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/df.net.list.rds"),

#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Chan et al.2021/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/df.net.list.rds"),
#   readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/df.net.list.rds")
# )

# LUAD.patients.size = list()
# LUSC.patients.size = list()
# 
# get.patients.celltypes.counts = function(net.list, patients) {
#   results = list()
#   for(patient in patients) {
#     results[[patient]] = length(unique(c(net.list[[patient]]$source, net.list[[patient]]$target)))
#   }
#   return(results)
# }
# LUAD.patients.size = get.patients.celltypes.counts(net.list, LUAD.patients)
# LUSC.patients.size = get.patients.celltypes.counts(net.list, LUSC.patients)
# 
# LUAD.patients.order = names(LUAD.patients.size)[order(unlist(LUAD.patients.size), decreasing = T)]
# LUSC.patients.order = names(LUSC.patients.size)[order(unlist(LUSC.patients.size), decreasing = T)]
# 
# print(LUAD.patients.order[1:6])
# print(LUSC.patients.order[1:6])
# 
# template.patients = list()
# for(i in 2:16) {
#   LUAD.template = get.template.net(net.list, LUAD.patients.order[1:i])
#   LUSC.template = get.template.net(net.list, LUSC.patients.order[1:i])
#   template.patients[[as.character(i)]] = c(LUAD.patients.order[1:i], LUSC.patients.order[1:i])
#   network.to.json(LUAD.template, C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/dataset/LungCancer/", paste0("LUAD_template_", i))
#   network.to.json(LUSC.template, C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/dataset/LungCancer/", paste0("LUSC_template_", i))
# }
# write_json(template.patients, pretty=F, auto_unbox=F, path = C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/dataset/template_info.json")
# 
# normal.net = readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/cellchat.rds")
# normal.net.size = get.patients.celltypes.counts(normal.net, names(normal.net))
# normal.net.order = names(normal.net.size)[order(unlist(normal.net.size), decreasing = T)]
# 
# network.to.json(normal.template, C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/dataset/", "Normal_template")
# 
# 
# 
# # 比较模板网络 -----------------------------------------------------------------
# 
# 
# 
# summary.interaction = function(template.difference, interaction.interested) {
#   interaction.interested.info = data.frame(matrix(nrow = 0, ncol = 7))
#   colnames(interaction.interested.info) = c("name", "ligand", "receptor", "celltype1", "celltype2", "pathway", "weight")
#   for(name in interaction.interested) {
#     lines = which(template.difference$Source.Target == name)
#     for(line in lines) {
#       ligand = sub("\\(.*\\)", "", template.difference$Source[line])
#       receptor = sub("\\(.*\\)", "", template.difference$Target[line])
#       celltype1 = template.difference$Celltype1[line]
#       celltype2 = template.difference$Celltype2[line]
#       pathway = template.difference$Pathway[line]
#       weight = template.difference$Weight[line]
#       interaction.interested.info[nrow(interaction.interested.info)+1, ] = c(name, ligand, receptor, celltype1, celltype2, pathway, weight)
#     }
#   }
#   interaction.interested.info = interaction.interested.info[order(interaction.interested.info$weight, decreasing = T), ]
#   return(interaction.interested.info)
# }
# 
# idenTypeUsingSummary = function(net.list, patients, Type.1.summary, Type.2.summary) {
#   results = c()
#   for(patient in patients) {
#     Type.1.score = 0
#     Type.2.score = 0
#     patient.net = get.testing.net(net.list, patient)
#     patient.net[["name"]] = paste0(patient.net$Source, "->", patient.net$Target)
#     for(interaction in patient.net[["name"]]) {
#       loc1 = which(Type.1.summary$name == interaction)
#       loc2 = which(Type.2.summary$name == interaction)
#       if(length(loc1) != 0) {
#         Type.1.score = Type.1.score + 1
#       }
#       if(length(loc2) != 0) {
#         Type.2.score = Type.2.score + 1
#       }
#     }
#     Type.1.score = Type.1.score / nrow(Type.1.summary)
#     Type.2.score = Type.2.score / nrow(Type.2.summary)
#     if(Type.1.score >= Type.2.score) {
#       results = c(results, 1)
#     } else {
#       results = c(results, 2)
#     }
#   }
#   return(results)
# }
# 
# compare.net = function(net.list, type.1.patients, type.2.patients, normal.net, normal.patients, template.size = 8) {
#   Normal.template.patients = normal.patients[1:template.size]
#   Type.1.template.patients = type.1.patients[1:template.size]
#   Type.2.template.patients = type.2.patients[1:template.size]
#   
#   Normal.template = get.template.net(normal.net, Normal.template.patients)
#   Type.1.template = get.template.net(net.list, Type.1.template.patients)
#   Type.2.template = get.template.net(net.list, Type.2.template.patients)
#   
#   Normal.template[["Source.Target"]] = paste0(Normal.template$Source, "->", Normal.template$Target)
#   Type.1.template[["Source.Target"]] = paste0(Type.1.template$Source, "->", Type.1.template$Target)
#   Type.2.template[["Source.Target"]] = paste0(Type.2.template$Source, "->", Type.2.template$Target)
#   
#   Normal.template = Normal.template[Normal.template$Weight >= (template.size / 2), ]
#   Type.1.template = Type.1.template[Type.1.template$Weight >= (template.size / 2), ]
#   Type.2.template = Type.2.template[Type.2.template$Weight >= (template.size / 2), ]
#   
#   Type.1.template.difference = Type.1.template
#   for(i in 1:nrow(Type.2.template)) {
#     loc = which(Type.1.template.difference$Source == Type.2.template$Source[i] & Type.1.template.difference$Target == Type.2.template$Target[i])
#     if(length(loc) != 0) {
#       Type.1.template.difference$Weight[loc] = Type.1.template.difference$Weight[loc] - Type.2.template$Weight[i]
#     }
#   }
#   
#   Type.2.template.difference = Type.2.template
#   for(i in 1:nrow(Type.1.template)) {
#     loc = which(Type.2.template.difference$Source == Type.1.template$Source[i] & Type.2.template.difference$Target == Type.1.template$Target[i])
#     if(length(loc) != 0) {
#       Type.2.template.difference$Weight[loc] = Type.2.template.difference$Weight[loc] - Type.1.template$Weight[i]
#     }
#   }
#   
#   Type.1.template.difference = Type.1.template.difference[Type.1.template.difference$Weight > 0, ]
#   Type.2.template.difference = Type.2.template.difference[Type.2.template.difference$Weight > 0, ]
#   
#   line.remove = c()
#   source = gsub("\\(Malignant\\)", "", Type.1.template.difference$Source)
#   target = gsub("\\(Malignant\\)", "", Type.1.template.difference$Target)
#   for(i in 1:nrow(Normal.template)) {
#     loc = which(source == Normal.template$Source[i] & target == Normal.template$Target[i])
#     if(length(loc) != 0){
#       line.remove = c(line.remove, loc)
#     }
#   }
#   Type.1.template.difference = Type.1.template.difference[-1*line.remove, ]
#   
#   line.remove = c()
#   source = gsub("\\(Malignant\\)", "", Type.2.template.difference$Source)
#   target = gsub("\\(Malignant\\)", "", Type.2.template.difference$Target)
#   for(i in 1:nrow(Normal.template)) {
#     loc = which(source == Normal.template$Source[i] & target == Normal.template$Target[i])
#     if(length(loc) != 0){
#       line.remove = c(line.remove, loc)
#     }
#   }
#   Type.2.template.difference = Type.2.template.difference[-1*line.remove, ]
#   
#   Type.1.interaction.interested = setdiff(Type.1.template.difference$Source.Target, Type.2.template.difference$Source.Target)
#   Type.2.interaction.interested = setdiff(Type.2.template.difference$Source.Target, Type.1.template.difference$Source.Target)
# 
#   Type.1.summary = summary.interaction(Type.1.template.difference, Type.1.interaction.interested)
#   Type.2.summary = summary.interaction(Type.2.template.difference, Type.2.interaction.interested)
#   
#   Type.1.summary = Type.1.summary[Type.1.summary$weight >= (template.size / 2), ]
#   Type.2.summary = Type.2.summary[Type.2.summary$weight >= (template.size / 2), ]
#   
#   Type.1.summary$weight = as.numeric(Type.1.summary$weight)
#   Type.2.summary$weight = as.numeric(Type.2.summary$weight)
# 
#   Type.1.results = idenTypeUsingSummary(net.list, type.1.patients, Type.1.summary, Type.2.summary)
#   Type.2.results = idenTypeUsingSummary(net.list, type.2.patients, Type.1.summary, Type.2.summary)
#   return(list(
#     F1.score = F1_Score(c(Type.1.results, Type.2.results), c(rep(1, length(Type.1.results)), rep(2, length(Type.2.results))), positive = 1),
#     mcc = mcc(c(Type.1.results, Type.2.results), c(rep(1, length(Type.1.results)), rep(2, length(Type.2.results)))),
#     Type.1.summary = Type.1.summary,
#     Type.2.summary = Type.2.summary
#   ))
# }
# 
# calc.radio.of.interaction = function(Type.summary, Type.1.patients, Type.2.patients) {
#   Type.summary[["Type.1.radio"]] = rep(0, nrow(Type.summary))
#   Type.summary[["Type.2.radio"]] = rep(0, nrow(Type.summary))
#   source = paste0(Type.summary$ligand, "(", Type.summary$celltype1, ")")
#   target = paste0(Type.summary$receptor, "(", Type.summary$celltype2, ")")
#   for(i in 1:nrow(Type.summary)) {
#     Type.summary[["Type.1.radio"]][i] = radio.of.interaction2(net.list, LUAD.patients.order, source[i], target[i])
#     Type.summary[["Type.2.radio"]][i] = radio.of.interaction2(net.list, LUSC.patients.order, source[i], target[i])
#   }
#   return(Type.summary)
# }
# 
# plot.score = function(results, range.of.patient) {
#   plot.data = data.frame("N"=names(results), "F1"=rep(0, length(results)), "MCC"=rep(0, length(results)))
#   for(i in 1:length(results)) {
#     plot.data[i, 2] = results[[i]]$F1.score
#     plot.data[i, 3] = results[[i]]$mcc
#   }
#   plot.data$N = as.numeric(plot.data$N)
#   
#   plot.data = data.frame(N = c(plot.data$N, plot.data$N), 
#                          Score = c(plot.data$F1, plot.data$MCC), 
#                          Type = c(rep("F1-score", nrow(plot.data)), rep("MCC", nrow(plot.data))))
#   
#   return(ggplot(plot.data, aes(x=N, y=Score, group=Type, color=Type)) + 
#            geom_line(linetype="solid", size=1.3) + geom_point(size=2) + 
#            scale_x_continuous(limits=c(min(range.of.patient), max(range.of.patient)), breaks=seq(min(range.of.patient), max(range.of.patient), 1)) + 
#            theme(legend.position = "bottom") +
#            xlab("患者数") + ylab("评分"))
# }
# 
# plot.network = function(source, target, edge.weights) {
#   data = data.frame(source=source, target=target)
#   graph = graph_from_data_frame(data, directed = TRUE)
#   E(graph)$weight <- edge.weights
#   print(plot(graph,  
#              layout=layout_nicely,
#              vertex.size=8,
#              vertex.shape='circle', 
#              vertex.color="yellow",
#              vertex.label=NULL,
#              vertex.label.cex=0.8,
#              vertex.label.color='black',
#              edge.arrow.size=0.5,
#              edge.width = 0.8,
#              edge.label=NA,
#              edge.color="black"))
# }
# 
# plot.network(Type.1.summary[Type.1.summary$pathway == "MK", ])
# summary.table = Type.1.summary
# 
# 
# 
# subtypes.results = list()
# for(i in 2:16) {
#   subtypes.results[[as.character(i)]] =  compare.net(net.list, LUAD.patients.order, LUSC.patients.order, normal.net, normal.net.order, i)
# }
# print(plot.score(subtypes.results, 2:16))
# 
# Type.1.summary = subtypes.results[[8]]$Type.1.summary
# Type.2.summary = subtypes.results[[8]]$Type.2.summary
# Type.1.summary[["source"]] = paste0(Type.1.summary$ligand, "(", Type.1.summary$celltype1, ")")
# Type.1.summary[["target"]] = paste0(Type.1.summary$receptor, "(", Type.1.summary$celltype2, ")")
# Type.1.summary = calc.radio.of.interaction(Type.1.summary, LUAD.patients.order, LUSC.patients.order)
# Type.2.summary = calc.radio.of.interaction(Type.2.summary, LUAD.patients.order, LUSC.patients.order)
# Type.1.summary[["rado.diff"]] = Type.1.summary$Type.1.radio - Type.1.summary$Type.2.radio
# Type.2.summary[["rado.diff"]] = Type.2.summary$Type.2.radio - Type.2.summary$Type.1.radio
# 
# setdiff(Type.1.summary$pathway, Type.2.summary$pathway)
# setdiff(Type.2.summary$pathway, Type.1.summary$pathway)
# 
# Type.1.patients.size = get.patients.celltypes.counts(net.list, N.stage$Y)
# Type.2.patients.size = get.patients.celltypes.counts(net.list, N.stage$N)
# Type.1.patients.order = names(Type.1.patients.size)[order(unlist(Type.1.patients.size), decreasing = T)]
# Type.2.patients.order = names(Type.2.patients.size)[order(unlist(Type.2.patients.size), decreasing = T)]
# 
# N.stage.results = list()
# for(i in 2:16) {
#   N.stage.results[[as.character(i)]] = compare.net(net.list, Type.1.patients.order, Type.2.patients.order, normal.net, normal.net.order, i)
# }
# print(plot.score(N.stage.results, 2:16))
# 
# 
# Type.1.summary = Type.1.summary[Type.1.summary$rado.diff > 0, ]
# Type.2.summary = Type.2.summary[Type.2.summary$rado.diff > 0, ]
# Type.1.summary[["source"]] = paste0(Type.1.summary$ligand, "(", Type.1.summary$celltype1, ")")
# Type.1.summary[["target"]] = paste0(Type.1.summary$receptor, "(", Type.1.summary$celltype2, ")")
# 
# Type.2.summary[["source"]] = paste0(Type.2.summary$ligand, "(", Type.2.summary$celltype1, ")")
# Type.2.summary[["target"]] = paste0(Type.2.summary$receptor, "(", Type.2.summary$celltype2, ")")
# 
# write.xlsx(Type.1.summary, file = "subtype_summary.xlsx", sheetName = "LUAD", row.names = F)
# write.xlsx(Type.2.summary, file = "subtype_summary.xlsx", sheetName = "LUSC", append = T, row.names = F)
# 


# ------------------------------------------------------------------------------

calc.radio.of.interaction2 = function(net.list, template, Type.1.patients, Type.2.patients) {
  template[["Type.1.radio"]] = rep(0, nrow(template))
  template[["Type.2.radio"]] = rep(0, nrow(template))
  template[["Type.1.offset"]] = rep(0, nrow(template))
  template[["Type.2.offset"]] = rep(0, nrow(template))

  patients.networks = list()
  for(patient in c(Type.1.patients, Type.2.patients)) {
    patients.networks[[patient]] = data.frame(Source=paste0(net.list[[patient]]$ligand, "(", net.list[[patient]]$source, ")"),
                                              Target=paste0(net.list[[patient]]$receptor, "(", net.list[[patient]]$target, ")"))
  }
  for(i in 1:nrow(template)) {
    for(patient in Type.1.patients) {
      if(sum(c(template$Celltype1[i], template$Celltype2[i]) %in% patient.celltypes[[patient]]) == 2) {
        if(sum(patients.networks[[patient]]$Source == template$Source[i] & patients.networks[[patient]]$Target == template$Target[i]) > 0) {
          template[["Type.1.radio"]][i] = template[["Type.1.radio"]][i] + 1
        }
      } else {
        template[["Type.1.offset"]][i] = template[["Type.1.offset"]][i] + 1
      }
    }
    for(patient in Type.2.patients) {
      if(sum(c(template$Celltype1[i], template$Celltype2[i]) %in% patient.celltypes[[patient]]) == 2) {
        if(sum(patients.networks[[patient]]$Source == template$Source[i] & patients.networks[[patient]]$Target == template$Target[i]) > 0) {
          template[["Type.2.radio"]][i] = template[["Type.2.radio"]][i] + 1
        }
      } else {
        template[["Type.2.offset"]][i] = template[["Type.2.offset"]][i] + 1
      }
    }
  }
  return(template)
}

idenTypeUsingFeature = function(net.list, patients, Type.1.features, Type.2.features) {
  results = c()
  for(patient in patients) {
    Type.1.score = 0
    Type.2.score = 0
    patient.net = get.testing.net(net.list, patient)
    
    for(i in 1:nrow(Type.1.features)) {
      if(sum(patient.net$Source == Type.1.features$Source[i] & patient.net$Target == Type.1.features$Target[i]) == 1) {
        Type.1.score = Type.1.score + 1
      }
    }
    for(i in 1:nrow(Type.2.features)) {
      if(sum(patient.net$Source == Type.2.features$Source[i] & patient.net$Target == Type.2.features$Target[i]) == 1) {
        Type.2.score = Type.2.score + 1
      }
    }
    Type.1.score = Type.1.score / nrow(Type.1.features)
    Type.2.score = Type.2.score / nrow(Type.2.features)
    # results[[patient]] = c(Type.1.score, Type.2.score)
    results[[patient]] = ifelse(Type.1.score >= Type.2.score, 1, 2)
  }
  return(results)
}

findFeatures = function(net.list, Type.1.patients, Type.2.patients) {
  template.all = get.template.net(net.list, c(Type.1.patients, Type.2.patients))
  template.all = calc.radio.of.interaction2(net.list, template.all, Type.1.patients, Type.2.patients)

  template.all[["radio.diff"]] = template.all$Type.1.radio - template.all$Type.2.radio
  
  Type.1.feature = template.all[template.all[["radio.diff"]] >= quantile(template.all[["radio.diff"]][template.all[["radio.diff"]] > 0.1], 0.5), ]
  Type.2.feature = template.all[template.all[["radio.diff"]] <= quantile(template.all[["radio.diff"]][template.all[["radio.diff"]] < -0.1], 0.5), ]
  # Type.1.feature = template.all[template.all[["radio.diff"]] >= 0.1, ]
  # Type.2.feature = template.all[template.all[["radio.diff"]] <= -0.1, ]
  predict = c(unlist(idenTypeUsingFeature(net.list, Type.1.patients, Type.1.feature, Type.2.feature)), 
              unlist(idenTypeUsingFeature(net.list, Type.2.patients, Type.1.feature, Type.2.feature)))
  return(list(
    F1=F1_Score(predict, c(rep(1, length(Type.1.patients)), rep(2, length(Type.2.patients)))),
    MCC=mcc(predict, c(rep(1, length(Type.1.patients)), rep(2, length(Type.2.patients)))),
    Type.1.feature=Type.1.feature,
    Type.2.feature=Type.2.feature
  ))
}

findFeatures2 = function(net.list, Type.1.patients, Type.2.patients, radio.threshold = 0.1, pvalue.threshold = 0.05) {
  template.all = get.template.net(net.list, c(Type.1.patients, Type.2.patients))
  template.all = calc.radio.of.interaction2(net.list, template.all, Type.1.patients, Type.2.patients)
  
  template.all = template.all[template.all$Type.1.offset < (length(Type.1.patients) * 0.75) & template.all$Type.2.offset < (length(Type.2.patients) * 0.75), ]

  print(paste0("排除前: ", nrow(template.all)))
  
  template.all = template.all[template.all$Type.1.radio > ((length(Type.1.patients) - template.all[["Type.1.offset"]]) * radio.threshold) | template.all$Type.2.radio > ((length(Type.2.patients) - template.all[["Type.2.offset"]]) * radio.threshold), ]

  print(paste0("排除后: ", nrow(template.all)))
  all.interaction = template.all
  template.all[["pvalue"]] = rep(1, nrow(template.all))
  n1 = length(Type.1.patients)
  n2 = length(Type.2.patients)
  for(i in 1:nrow(template.all)) {
    data = matrix(c(template.all$Type.1.radio[i], n1 - template.all[["Type.1.offset"]][i] - template.all$Type.1.radio[i],
                    template.all$Type.2.radio[i], n2 - template.all[["Type.2.offset"]][i] - template.all$Type.2.radio[i]), nrow = 2)
    template.all[["pvalue"]][i] = fisher.test(data)$p.value
  }
  template.all = template.all[template.all$pvalue < pvalue.threshold, ]

  template.all$Type.1.radio = template.all$Type.1.radio / (length(Type.1.patients) - template.all[["Type.1.offset"]])
  template.all$Type.2.radio = template.all$Type.2.radio / (length(Type.2.patients) - template.all[["Type.2.offset"]])
  
  Type.1.feature = template.all[template.all$Type.1.radio > template.all$Type.2.radio, ]
  Type.2.feature = template.all[template.all$Type.1.radio < template.all$Type.2.radio, ]
  
  Type.1.feature = Type.1.feature[order(Type.1.feature$pvalue), ]
  Type.2.feature = Type.2.feature[order(Type.2.feature$pvalue), ]
  
  predict = c(unlist(idenTypeUsingFeature(net.list, Type.1.patients, Type.1.feature, Type.2.feature)), 
              unlist(idenTypeUsingFeature(net.list, Type.2.patients, Type.1.feature, Type.2.feature)))
  return(list(
    F1=F1_Score(predict, c(rep(1, length(Type.1.patients)), rep(2, length(Type.2.patients)))),
    MCC=mcc(predict, c(rep(1, length(Type.1.patients)), rep(2, length(Type.2.patients)))),
    Type.1.feature=Type.1.feature,
    Type.2.feature=Type.2.feature,
    all.interaction=all.interaction
  ))
}

featuresToXlsx = function(features, Type.1.name, Type.2.name, path) {
  features$Type.1.feature[["radio.diff"]] = features$Type.1.feature$Type.1.radio - features$Type.1.feature$Type.2.radio
  features$Type.2.feature[["radio.diff"]] = features$Type.2.feature$Type.1.radio - features$Type.2.feature$Type.2.radio
  Type.1 = data.frame("源"=features$Type.1.feature$Source, 
                      "目标"=features$Type.1.feature$Target,
                      "通路"=features$Type.1.feature$Pathway, 
                      "A" = features$Type.1.feature$Type.1.radio,
                      "B" = features$Type.1.feature$Type.2.radio,
                      "P值" = features$Type.1.feature$pvalue,
                      "Diff"=features$Type.1.feature$radio.diff)
  colnames(Type.1) = c("源", "目标", "通路", Type.1.name, Type.2.name, "P值", "比例差")
  Type.1 = Type.1[order(Type.1[["P值"]], decreasing = F), ]
  # Type.1[[Type.1.name]] = paste0(round(Type.1[[Type.1.name]] * 100, 2), "%")
  # Type.1[[Type.2.name]] = paste0(round(Type.1[[Type.2.name]] * 100, 2), "%")
  
  Type.2 = data.frame("源"=features$Type.2.feature$Source, 
                      "目标"=features$Type.2.feature$Target,
                      "通路"=features$Type.2.feature$Pathway, 
                      "A" = features$Type.2.feature$Type.1.radio,
                      "B" = features$Type.2.feature$Type.2.radio,
                      "P值" = features$Type.2.feature$pvalue,
                      "Diff"=features$Type.2.feature$radio.diff)
  colnames(Type.2) = c("源", "目标", "通路", Type.1.name, Type.2.name, "P值", "比例差")
  Type.2 = Type.2[order(Type.2[["P值"]], decreasing = F), ]
  # Type.2[[Type.1.name]] = paste0(round(Type.2[[Type.1.name]] * 100, 2), "%")
  # Type.2[[Type.2.name]] = paste0(round(Type.2[[Type.2.name]] * 100, 2), "%")
  
  write.xlsx(Type.1, file = path, sheetName = Type.1.name, row.names = F)
  write.xlsx(Type.2, file = path, sheetName = Type.2.name, append = T, row.names = F)
}

features.to.input = function(net.list, features, path) {
  tree.input = data.frame(matrix(0, 
                                 nrow = length(net.list), 
                                 ncol = (nrow(features$Type.1.feature) + nrow(features$Type.2.feature))))
  colnames(tree.input) = c(paste0(features$Type.1.feature$Source, "->", features$Type.1.feature$Target),
                           paste0(features$Type.2.feature$Source, "->", features$Type.2.feature$Target))
  rownames(tree.input) = names(net.list)
  for(patient in names(net.list)) {
    patient.net = get.testing.net(net.list, patient)
    for(i in 1:nrow(features$Type.1.feature)) {
      if(sum(patient.net$Source == features$Type.1.feature$Source[i] & patient.net$Target == features$Type.1.feature$Target[i]) != 0) {
        tree.input[patient, paste0(features$Type.1.feature$Source[i], "->", features$Type.1.feature$Target[i])] = 1
      }
    }
    for(i in 1:nrow(features$Type.2.feature)) {
      if(sum(patient.net$Source == features$Type.2.feature$Source[i] & patient.net$Target == features$Type.2.feature$Target[i]) != 0) {
        tree.input[patient, paste0(features$Type.2.feature$Source[i], "->", features$Type.2.feature$Target[i])] = 1
      }
    }
  }
  write.xlsx(tree.input, file = path)
}

plot.pie = function(label, counts, path, min.percentage = 0.01) {
  png(file=path, width=1200, height=1200, res=120, type='cairo')
  par(mfrow = c(1, 1), mai = c(0, 0.6, 0, 1))
  fitter = (counts / sum(counts)) >= min.percentage
  pie(counts[fitter], labels = paste0(label, " (", round(counts / sum(counts) * 100, 2), "%)")[fitter], cex=1.2)
  dev.off()
}

net.list = readRDS("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/net_list_UCRSI_without_cnv.rds")
patients_info = read_json("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/patients_info.json")
patients_N_stage = list(
  "Normal" = as.character(unlist(patients_info[["N stage"]][["Normal"]])),
  "LUAD_N0" = c(),
  "LUAD_N1|N2|N3" = c(),
  "LUSC_N0" = c(),
  "LUSC_N1|N2|N3" = c()
)
for(datasets in names(patients_info[["N stage"]][["N0"]])){
  for(stage in c("N0", "N1|N2|N3")) {
    patients = patients_info[["N stage"]][[stage]][[datasets]]
    patients_N_stage[[paste0("LUAD_", stage)]] = c(patients_N_stage[[paste0("LUAD_", stage)]],
                                                   as.character(patients[patients %in% patients_info[["Subtypes"]][["LUAD"]][[datasets]]]))
    patients_N_stage[[paste0("LUSC_", stage)]] = c(patients_N_stage[[paste0("LUSC_", stage)]],
                                                   as.character(patients[patients %in% patients_info[["Subtypes"]][["LUSC"]][[datasets]]]))
  }
}

LUAD.Normal.N0 = findFeatures2(net.list, patients_N_stage$Normal, patients_N_stage$LUAD_N0, 0.1, 0.05)
LUAD.N0.N123 = findFeatures2(net.list, patients_N_stage$LUAD_N0, patients_N_stage$`LUAD_N1|N2|N3`, 0.1, 0.05)

LUSC.Normal.N0 = findFeatures2(net.list, patients_N_stage$Normal, patients_N_stage$LUSC_N0, 0.1, 0.05)
LUSC.N0.N123 = findFeatures2(net.list, patients_N_stage$LUSC_N0, patients_N_stage$`LUSC_N1|N2|N3`, 0.1, 0.05)

save(LUAD.Normal.N0, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD.Normal.N0.rds")
save(LUAD.N0.N123, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD.N0.N123.rds")
save(LUSC.Normal.N0, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC.Normal.N0.rds")
save(LUSC.N0.N123, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC.N0.N123.rds")

featuresToXlsx(LUAD.Normal.N0, "Normal", "N0", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD_Normal_N0.xlsx")
featuresToXlsx(LUAD.N0.N123, "N0", "N123", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD_N0_N123.xlsx")
featuresToXlsx(LUSC.Normal.N0, "Normal", "N0", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC_Normal_N0.xlsx")
featuresToXlsx(LUSC.N0.N123, "N0", "N123", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC_N0_N123.xlsx")

LUAD.Normal.N0.pathway = as.data.frame(summary(as.factor(LUAD.Normal.N0$Type.2.feature$Pathway)))
LUAD.N0.N123.pathway.N0 = as.data.frame(summary(as.factor(LUAD.N0.N123$Type.1.feature$Pathway)))
LUAD.N0.N123.pathway.N123 = as.data.frame(summary(as.factor(LUAD.N0.N123$Type.2.feature$Pathway)))

LUSC.Normal.N0.pathway = as.data.frame(summary(as.factor(LUSC.Normal.N0$Type.2.feature$Pathway)))
LUSC.N0.N123.pathway.N0 = as.data.frame(summary(as.factor(LUSC.N0.N123$Type.1.feature$Pathway)))
LUSC.N0.N123.pathway.N123 = as.data.frame(summary(as.factor(LUSC.N0.N123$Type.2.feature$Pathway)))

write.xlsx(LUAD.Normal.N0.pathway, sheetName = "LUAD Normal-N0", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = F)
write.xlsx(LUAD.N0.N123.pathway.N0, sheetName = "LUAD N0-N1|N2|N3 N0特征", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = T)
write.xlsx(LUAD.N0.N123.pathway.N123, sheetName = "LUAD N0-N1|N2|N3 NN1|N2|N3特征", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = T)

write.xlsx(LUSC.Normal.N0.pathway, sheetName = "LUSC Normal-N0", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = T)
write.xlsx(LUSC.N0.N123.pathway.N0, sheetName = "LUSC N0-N1|N2|N3 N0特征", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = T)
write.xlsx(LUSC.N0.N123.pathway.N123, sheetName = "LUSC N0-N1|N2|N3 NN1|N2|N3特征", file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/pathway.xlsx", append = T)



plot.pie(rownames(LUAD.Normal.N0.pathway), LUAD.Normal.N0.pathway[, 1], "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD.Normal.N0.png")
plot.pie(rownames(LUAD.N0.N123.pathway.N123), LUAD.N0.N123.pathway.N123[, 1], "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUAD.N0.N123.png")

plot.pie(rownames(LUSC.Normal.N0.pathway), LUSC.Normal.N0.pathway[, 1], "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC.Normal.N0.png")
plot.pie(rownames(LUSC.N0.N123.pathway.N123), LUSC.N0.N123.pathway.N123[, 1], "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/结果/LUSC.N0.N123.png")



subtypes.features = findFeatures2(net.list, LUAD.patients, LUSC.patients, 0.1, 0.05)
subtypes.features$Type.1.feature = subtypes.features$Type.1.feature[order(subtypes.features$Type.1.feature$pvalue), ]
subtypes.features$Type.2.feature = subtypes.features$Type.2.feature[order(subtypes.features$Type.2.feature$pvalue), ]
featuresToXlsx(subtypes.features, "LUAD", "LUSC", "subtypes_summary.xlsx")
features.to.input(net.list, subtypes.features, "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/features/subtypes_input.xlsx")

N.stage.features = findFeatures2(net.list, N.stage$N, N.stage$Y, 0.1, 0.05)
features.to.input(net.list, N.stage.features, "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/features/N_stage_input.xlsx")

LUAD.N.stage.features = findFeatures2(net.list, LUAD.N.stage$N, LUAD.N.stage$Y, 0.1, 0.05)
LUAD.N.stage.features$Type.1.feature = LUAD.N.stage.features$Type.1.feature[order(LUAD.N.stage.features$Type.1.feature$pvalue), ]
LUAD.N.stage.features$Type.2.feature = LUAD.N.stage.features$Type.2.feature[order(LUAD.N.stage.features$Type.2.feature$pvalue), ]

LUSC.N.stage.features = findFeatures2(net.list, LUSC.N.stage$N, LUSC.N.stage$Y, 0.1, 0.05)
LUSC.N.stage.features$Type.1.feature = LUSC.N.stage.features$Type.1.feature[order(LUSC.N.stage.features$Type.1.feature$pvalue), ]
LUSC.N.stage.features$Type.2.feature = LUSC.N.stage.features$Type.2.feature[order(LUSC.N.stage.features$Type.2.feature$pvalue), ]

summary.N.stage.features = list(
  "Type.1.feature" = rbind(LUAD.Normal.N0$Type.1.feature, LUSC.Normal.N0$Type.1.feature),
  "Type.2.feature" = rbind(LUAD.N.stage.features$Type.2.feature, LUSC.N.stage.features$Type.2.feature)
)

Type.1.feature = rbind(LUAD.N0.N123$Type.1.feature, LUSC.N0.N123$Type.1.feature)
Type.1.feature = Type.1.feature[!duplicated(Type.1.feature$Source.Target), ]
Type.2.feature = rbind(LUAD.N0.N123$Type.2.feature, LUSC.N0.N123$Type.2.feature)
Type.2.feature = Type.2.feature[!duplicated(Type.2.feature$Source.Target), ]
summary.N.stage.features = list(
  "Type.1.feature" = Type.1.feature,
  "Type.2.feature" = Type.2.feature
)

sum(idenTypeUsingFeature(net.list, 
                         c(patients_N_stage$LUAD_N0, patients_N_stage$LUSC_N0), 
                         summary.N.stage.features$Type.1.feature,
                         summary.N.stage.features$Type.2.feature) == 1)
sum(idenTypeUsingFeature(net.list, 
                         c(patients_N_stage$`LUAD_N1|N2|N3`, patients_N_stage$`LUSC_N1|N2|N3`), 
                         summary.N.stage.features$Type.1.feature,
                         summary.N.stage.features$Type.2.feature) == 2)


sum(c(unlist(idenTypeUsingFeature(net.list, N.stage$N, summary.N.stage.features$Type.1.feature, summary.N.stage.features$Type.2.feature)) == 1, 
      unlist(idenTypeUsingFeature(net.list, N.stage$Y, summary.N.stage.features$Type.1.feature, summary.N.stage.features$Type.2.feature)) == 2))

features.to.input(net.list, summary.N.stage.features, "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/features/summary_N_stage_input.xlsx")

featuresToXlsx(LUAD.N.stage.features, "N0", "N1|N2|N3", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/LUAD_N_stage_summary.xlsx")
featuresToXlsx(LUSC.N.stage.features, "N0", "N1|N2|N3", "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/LUSC_N_stage_summary.xlsx")


# featuresToXlsx(N.stage.features, "N0", "N1|N2|N3", "N_stage_summary.xlsx")

# idenTypeUsingFeature(net.list, names(Song), 
#                      N.stage.features$Type.1.feature,
#                      N.stage.features$Type.2.feature)
# 
# 


# 
# 
# 
features.to.input(net.list, subtypes.features, "./FNN/features/subtypes_input.xlsx")
# features.to.input(N.stage.features, "./FNN/features/N_stage_input.xlsx")
# 
# saveRDS(subtypes.features, "./data/subtypes.features.rds")
# saveRDS(N.stage.features, "./data/N.stage.features.rds")


# GSE131907 --------------------------------------------------------------------
# cellchat.obj = readRDS(C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/cellchat.obj.rds")
# datasets.net = list()
# for(patient in names(cellchat.obj)) {
#   if(patient == "LT_S13") {
#     next
#   }
#   datasets.net[[patient]] = subsetCommunication(cellchat.obj[[patient]])
# }
# 
# idenTypeUsingFeature(datasets.net, names(datasets.net), 
#                      subtypes.features$Type.1.feature,
#                      subtypes.features$Type.2.feature)
# idenTypeUsingFeature(datasets.net, names(datasets.net), 
#                      N.stage.features$Type.1.feature,
#                      N.stage.features$Type.2.feature)
