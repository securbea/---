# 安装并加载 VennDiagram 包
install.packages("VennDiagram")
library(VennDiagram)

# 设置工作目录
workingDir <- "C:/Users/12345"
setwd(workingDir)

# 不同数量的vene图绘制并无大区别，只需保证元素存在，在list中修改参数即可
# 这里举例两个元素和四个元素的绘制

# 两个元素
proflies <- read.csv("COL5A3.csv", header = TRUE, row.names = 1)
print(hubgenes)
# 创建 Venn 图
venn.plot <- venn.diagram(
  x = list(
    "proflies" = rownames(hubgenes),
    "DEGs" = rownames(c)
  ),
  filename = NULL, # 不再保存为文件
  col = "transparent",
  fill = c("cornflowerblue", "green"),
  alpha = 0.50,
  label.col = "red",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0.2
)

# 绘制图形
grid.draw(venn.plot)


# 四个元素
# 读取数据
hubgenes <- read.csv("hubgene.csv", header = TRUE, row.names = 1)
a <- read.csv("90-100.csv", header = TRUE, row.names = 1)
b <- read.csv("90-110.csv", header = TRUE, row.names = 1)
c <- read.csv("90-120.csv", header = TRUE, row.names = 1)
print(hubgenes)
# 创建 Venn 图
venn.plot <- venn.diagram(
  x = list(
    "Hubgenes" = rownames(hubgenes),
    "A/B" = rownames(a),
    "A/C" = rownames(b),
    "A/D" = rownames(c)
  ),
  filename = NULL, # 不再保存为文件
  col = "transparent",
  fill = c("cornflowerblue", "green", "yellow", "darkorchid1"),
  alpha = 0.50,
  label.col = "white",
  cex = 1.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("darkblue", "darkgreen", "orange", "darkorchid4"),
  cat.cex = 1.5,
  cat.pos = 0,
  cat.dist = 0.07,
  cat.fontfamily = "serif",
  rotation.degree = 270,
  margin = 0.2
)

# 绘制图形
grid.draw(venn.plot)

