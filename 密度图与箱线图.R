#本文档提供两种不同的图像保存方式
#一种是直接画图，在右侧会直接出现预览
#一种是将图像存储在变量中，手动输出图片后查看
#根据个人习惯选择，本库中其他文件默认为直接画图

library(data.table)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)

# 禁止Y轴科学计数法显示，可根据需求删除
options(scipen=200) 

#设置工作路径
setwd("D:/黄小马/Documents/")


#读取数据
len <-  fread("distance.csv")

#整理数据
# 首先将数据框转换为长格式
library(tidyr)
genes_long <- pivot_longer(genes, cols = c(Upstream_neighbor, Downstream_neighbor), 
                           names_to = "gene_type", values_to = "gene_length")
colnames(genes)
names(genes)

# 直接绘制图形
ggplot(data = len_data, aes(x = length, group = RNA, fill = RNA)) +
  geom_density(adjust = 1.5, alpha = .4, linewidth = 0.6) +#调整曲线平滑度、峰透明度、线条宽度
  scale_x_continuous(expand = c(0, 0), limits = c(0, 150000), 
                     breaks = seq(0, 150000, 25000)) + #调整x轴与y轴的显示范围（会将此范围外的数值过滤掉）
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.00003), 
                     breaks = seq(0, 0.00003, 0.00002)) +#如不想过滤数值，需删除
  labs(x = "Transcript_length", y = "density") +
  theme_classic()#如需保留网格线，使用theme_minimal()

#将图像保存至tiff文件
ggsave("density_plot.tiff", width = 10, height = 8, dpi = 300)

#关闭图像设备：使用图像函数后，需关闭图像设备以保证图像能够正常保存
dev.off()


n_exon <- fread("lncRNA-mRNA-exonnumber.txt")

exon_data_lnc <- data.frame(V1 = n_exon$LncRNA[!is.na(n_exon$LncRNA)], V2 =
                              "LncRNA")
exon_data_m <- data.frame(V1 = n_exon$mRNA[!is.na(n_exon$mRNA)], V2 = "mRNA")
exon_data <- rbind(exon_data_lnc, exon_data_m)
names(exon_data) <- c("n_exon", "RNA")

#先将图像储存在变量中，再将变量导出为图像保存
p2 <- ggplot(data = exon_data, aes(x = n_exon, group = RNA, fill = RNA)) +
  geom_density(adjust = 1.5, alpha = .4, size = 0.5) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 30), breaks = seq(0, 30, 5)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2.3), breaks = seq(0, 2.0, 0.5)) +
  labs(x = "Exon_Number", y = "density") +
  theme_classic()

png(filename = "Exon_Number.png", width = 1975, height = 1271, res = 300)
print(p2)
dev.off()


#用箱线图查看异常值
#绘制箱线图
summary(genes_long$gene_length)
range(genes_long$gene_length, na.rm = TRUE)
boxplot(genes_long$gene_length, main = "Boxplot of Gene Lengths")


