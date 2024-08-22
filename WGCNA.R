rm(list = ls())
workingDir = "C:/Users/"
setwd(workingDir)
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = FALSE)
enableWGCNAThreads()
dataExpr <- read.csv("SSH222-17.csv", header = TRUE, row.names = 1)
dim(dataExpr)

# 检查数据框的列数
ncol_dataExpr <- ncol(dataExpr)
print(ncol_dataExpr)
head(dataExpr)[,1:16]
## 筛选中位绝对偏差(MAD)前50%的基因，且至少MAD大于0.01
##probs=seq(0, 1, 0.25)指定计算0%、25%、50%、75%和100%的分位数
##然后通过索引[2]选择中位数。
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
## 转换为样品在行，基因在列的矩阵
dataExpr0 <- as.data.frame(t(dataExpr))

## 检测缺失值
gsg = goodSamplesGenes(dataExpr0, verbose = 3);
if (!gsg$allOK) {
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes) > 0) {
    print(paste("Removing genes:", 
                paste(names(dataExpr0)[!gsg$goodGenes], collapse = ",")))
  }
  
  if (sum(!gsg$goodSamples) > 0) {
    print(paste("Removing samples:", 
                paste(rownames(dataExpr0)[!gsg$goodSamples], collapse = ",")))
  }
  
  # Remove the offending genes and samples from the data:
  dataExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr0)
nSamples = nrow(dataExpr0)

dim(dataExpr0)
head(dataExpr0)
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
# 使用 hclust 对样品进行层次聚类，然后通过 plot 函数可视化聚类结果。
#这段代码是用于检查是否有离群样品，步骤如下：
#计算样品之间的距离： 使用 dist(dataExpr) 计算基因表达数据中样品之间的欧氏距离。
#进行层次聚类： 使用 hclust 函数对样品进行层次聚类，方法为 "average"。
#这将生成一个样品层次聚类树，其中样品根据它们的表达式模式被分组。
#可视化聚类树： 使用 plot 函数将样品层次聚类树可视化，
#标题为 "Sample clustering to detect outliers"。
#制表符分隔(两种输入方式根据自己的文件类型选择)
myPhenotypeData <- read.table("ssh22-17-2.txt", header = TRUE, sep = "\t", row.names = 1)
#逗号分隔
myPhenotypeData <- read.table("ssh22-17-2.txt", header = TRUE, row.names = 1)  # 如果是制表符分隔，使用 "\t"；如果是逗号分隔，使用 ","
datTraits <- myPhenotypeData
sampleTree2 = hclust(dist(dataExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

powers = c(c(1:10), seq(from = 10, to=30, by=2))
sft = pickSoftThreshold(dataExpr0, powerVector=powers, verbose=5)

par(mfrow = c(1,2))
# 将画布一分为二开始绘图
cex1 = 0.9
# 将文本调整为0.9倍
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# 筛选标准：加一条水平线。R-square=0.85（根据自己的需求选择）
abline(h=0.85,col="red")

# Soft threshold与平均连通性
# 其中 x 轴是软阈值的幂次，y 轴是与软阈值相对应的平均连通性。
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
sft$powerEstimate

power_select=sft[['powerEstimate']]
power=power_select

#有时候计算结果是NA，只需自己合理选择软阈值即可
net = blockwiseModules(dataExpr0, power = 16,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
table(net$colors)
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save.image(file = "workspace.RData")

# 在下次需要时加载工作区
load("workspace.RData")

nGenes = ncol(dataExpr0);
nSamples = nrow(dataExpr0);
# 重新计算带有颜色标签的模块
MEs0 = moduleEigengenes(dataExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
# 通过相关值对每个关联进行颜色编码
sizeGrWindow(10,6)
# 展示模块与表型数据的相关系数和 P值
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# 用热图的形式展示相关系数
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#此处的L为表型变量，根据自己的数据修改
L = as.data.frame(datTraits$L);
names(L) = "L";
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(dataExpr0, L, use = "p"));#和体重性状的关联
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(L), sep="");
names(GSPvalue) = paste("p.GS.", names(L), sep="")

#本次和表型相关的模块
module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;
table(moduleGenes)
red_module<-as.data.frame(dimnames(data.frame(dataExpr0))[[2]][moduleGenes])
write.csv(red_module,"redmodulegene.csv")


sizeGrWindow(7, 7)
par(mfrow = c(1,1))
point_color <- "red"  # 你可以根据需要更改颜色

#绘制相关性散点图
verboseScatterplot(abs(geneTraitSignificance[moduleGenes, 1]),
                   abs(geneModuleMembership[moduleGenes, column]),
                   xlab = "Gene significance for L",
                   ylab = paste("Module Membership in", module, "module"),
                   main = paste("Gene significance vs. module membership\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = point_color,
                   pch = 1,  # 设置点形状为空心圆
                   xlim = c(0, max(abs(geneTraitSignificance[moduleGenes, 1]))),
                   ylim = c(0, max(abs(geneModuleMembership[moduleGenes, column]))))

#可选编辑，非必要
# 在条件处画上表示截断的虚线
# 在条件处画上表示截断的蓝色线
abline(h = 0.8, col = "skyblue", lty = 2, lwd = 2)  # 水平虚线
abline(v = 0.2, col = "skyblue", lty = 2, lwd = 2)  # 垂直虚线

#计算TOM的时间较长
TOM = TOMsimilarityFromExpr(dataExpr0, power = 16)

# Choose modules
modules <- c("blue")
probes <- names(dataExpr0)
inModule <- is.finite(match(moduleColors, modules))
modProbes <- probes[inModule]

# Choose relevant TOM matrix
modTOM <- TOM[inModule, inModule]
dimnames(modTOM) <- list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
#绘制cytoscape互作图（不同于string绘制的蛋白互作分析图）
cyt <- exportNetworkToCytoscape(modTOM,
                                edgeFile = paste("CytoscapeInput-edges-ablue", paste(modules, collapse="-"), ".txt", sep=""),
                                nodeFile = paste("CytoscapeInput-nodes-ablue", paste(modules, collapse="-"), ".txt", sep=""),
                                weighted = TRUE,
                                threshold = 0.02,
                                nodeNames = modProbes,
                                nodeAttr = moduleColors[inModule])

save.image(file = "my_session.RData")

#进一步操作
#筛选模块内的高表达基因进行保存和分析
module = "cyan"
column = match(module, modNames)
moduleGenes = moduleColors == module

# 提取模块内基因的名称
cyan_module <- as.data.frame(dimnames(data.frame(dataExpr0))[[2]][moduleGenes])
names(cyan_module) <- "genename"

# 计算基因模块成员值（MM）和性状显著性（GS）
MM <- abs(geneModuleMembership[moduleGenes, column])
GS <- abs(geneTraitSignificance[moduleGenes, 1])

# 将 MM 和 GS 合并到数据框中，并添加基因名称作为行名
c <- as.data.frame(cbind(MM, GS))
rownames(c) <- cyan_module$genename

# 筛选出 MM > 0.8 且 GS > 0.2 的核心基因
cyan_hub <- c[MM > 0.8 & GS > 0.2, ]

# 输出核心基因名称
hub_genes <- rownames(cyan_hub)
print(hub_genes)

# 保存符合条件的基因
write.csv(hub_genes, "hubgene_MMGS_cyan.csv")

# 如这部分的基因需要注释，可将其重新match到注释表上
# 返回属于颜色模块的基因ID
names(dataExpr0)[moduleColors=="cyan"] 

annot = read.csv(file = "geneanovation.csv");
dim(annot)
names(annot)
probes = names(dataExpr0) # 匹配信息
probes2annot = match(probes, annot$Gene_name);
sum(is.na(probes2annot)) # 检测是否有没有匹配上的ID号，正常来说为0，即全匹配上了。
#输出必要的信息：
geneInfo0 = data.frame(Gene_name = probes,
                       Gene_description = annot$Gene_description[probes2annot],
                       Gene_id = annot$Gene_id[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue);
#按照与体重的显著水平将模块进行排序:
modOrder = order(-abs(cor(MEs, group100, use = "p")))
#添加模块成员的信息：
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.group));  # 排序
geneInfo = geneInfo0[geneOrder, ]
#输出为CSV格式，可用fix(geneInfo)在R中查看：
write.csv(geneInfo, file = "geneInfo.csv")