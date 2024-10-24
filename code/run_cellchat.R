seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Zhang et al.2022/Zhang.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Zhang et al.2022/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Leader et al.2021/Leader.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Leader et al.2021/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/GSE148071.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/Maynard.seurat.subset.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/E_MTAB_6149.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/GSE127465.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/GSE117570.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE153935/GSE153935.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE153935/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/GSE131907.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/UKIM-V/UKIM.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/UKIM-V/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/PRJEB31843/PRJEB31843.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/PRJEB31843/UCRSI.cellchat.obj.rds"

seurat.file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/GSE136831.Control.seurat.rds"
net.save.path = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/UCRSI.cellchat.obj.rds"

library(CellChat)
cellchat.paths = c("C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Zhang et al.2022/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Leader et al.2021/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE148071/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/Maynard et al.2020/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/E-MTAB-6149/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE127465/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE117570/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE153935/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE131907/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/UKIM-V/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/PRJEB31843/UCRSI.cellchat.obj.rds",
                   "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/Lung Cancer/GSE136831/UCRSI.cellchat.obj.rds")
net.list = list()
for(path in cellchat.paths) {
  cellchat.list = readRDS(path)
  for(patient in names(cellchat.list)) {
    net.list[[patient]] = subsetCommunication(cellchat.list[[patient]])
  }
}
saveRDS(net.list, file = "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/data/net_list_UCRSI_without_cnv.rds")

for(patient in names(net.list)) {
  network.to.json(get.testing.net(net.list, patient), 
                  "C:/Users/tkken/Nextcloud/研究生/硕士毕业论文实验/FNN/dataset/LungCancer_UCRSI_without_cnv/", patient)
}
