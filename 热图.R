# 下载包#
install.packages("pheatmap")
install.packages("RColorBrewer")
rm(list = ls())
# 加载包#
library("pheatmap")
library("RColorBrewer")
library(ggplot2)
# 加载绘图数据#
data<-read.table(file='C:/Users/12345/Desktop/cell cycle.txt'
                 ,header=TRUE,row.names= 1)

data<-read.table(file='C:/Users/12345/Desktop/DEGs热图818.csv',
                 header=TRUE,row.names= 1,sep=',')
head(data) #查看数据
data=log2(data+1) #对基因表达量数据处理（所有列）
data <- as.matrix(data) #转变为matrix格式矩阵
head(data)

# 基础热图绘制
pheatmap(data) 
data[is.infinite(data)] <- NA  # 将Inf替换为NA
data[is.nan(data)] <- NA  # 将NaN替换为NA
data[is.na(data)] <- 0  # 将NA替换为0
sum(is.na(data))     # 检查是否还有NA值
sum(is.nan(data))    # 检查是否还有NaN值
sum(is.infinite(data))  # 检查是否还有Inf值

# 个性化设置
pheatmap(data,scale="row",
         
         )

# 个性化设置详细代码，选择需要的复制到上面

pheatmap(data, scale = "row", #表示进行均一化的方向，值为 “row”, “column” 或者"none"
      
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",#表示行、列聚类使用的度量方法，默认为欧式距离“euclidean”， "correlation"表示按照 Pearson correlation方法进行聚类
         
         clustering_method = "complete", #表示聚类方法，包括：‘ward’, ‘ward.D’, ‘ward.D2’, ‘single’, ‘complete’, ‘average’, ‘mcquitty’, ‘median’, ‘centroid’
         
         cluster_rows = T,cluster_cols = T, #cluster_rows表示仅对行聚类，cluster_cols表示仅对列聚类，值为TRUE或FALSE
         
         cutree_rows = NA, cutree_cols = NA, #若进行了行/列聚类，根据行/列聚类数量分隔热图行,cutree_rows=num分割行，cutree_cols=num分割列
         
         treeheight_row = 30, treeheight_col = 30, #若行、列聚类树高度调整
         
         border_color = "grey60", #表示热图每个小的单元格边框的颜色，默认为 "grey60"
         
         cellwidth = 60, cellheight = 0.39,  #表示单个单元格的宽度\高度，默认为 “NA”
         
         display_numbers = F, #表示是否在单元格上显示原始数值或按照特殊条件进行区分标记
         
         fontsize_number = 6, #表示热图上显示数字的字体大小
         
         number_format = "%.2f", #表示热图单元格上显示的数据格式，“%.2f” 表示两位小数,“%.1e”表示科学计数法
         
         number_color = "grey30", #表示热图单元格上显示的数据字体颜色
         
         fontsize =10, fontsize_row = 6, fontsize_col = 10, #热图中字体大小、行、列名字体大小
         
         show_rownames = F, show_colnames = T, #表示是否显示行名、列名
         
         main = "pheatmap of Cell cycle", #表示热图的标题名字
         
         color = colorRampPalette(c("navy","white","firebrick3"))(100), #表示热图颜色,(100)表示100个等级
         
         legend = T, #表示是否显示图例，值为TRUE或FALSE
         
         legend_breaks = NA, #设置图例的范围legend_breaks=c(-2.5,0,2.5)表示图例断点的设置，默认为NA,
         
         legend_labels = NA, #表示图例断点的标签
         
         angle_col = "45", #表示列标签的角度
         
         gaps_row = NULL,  #仅在未进行行聚类时使用，表示在行方向上热图的隔断位置
         
         gaps_col = c(1,2,3,4,5,6),  #仅在未进行列聚类时使用，表示在列方向上热图的隔断位置
         
         labels_row = NULL, #表示使用行标签代替行名
         
         labels_col = c("Sushan2","Sushan3","Sushan1","Huai2","Huai3","Huai1"),  #表示使用列标签代替列名
         
         filename = NA,  #表示保存图片的位置及命名
         
         width = NA, height = NA,
         
         margin = 0.1) #表示输出绘制热图的宽度/高度

ggsave(filename = "plot1.png", width = 50, height = 50)
