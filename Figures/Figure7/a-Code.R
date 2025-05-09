library(data.table)
library(plyr)
library(ggpubr)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(maxLik)
library(Rtsne)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(PerformanceAnalytics)
library(gplots)
library(DESeq2)
library(ellipse)
library(ggbiplot)
library(leaps)
library(plotROC)
library(lme4)
library(psych)
library(ade4)
library(ggplot2)
library(ggrepel)
library(qqman)
library(gridExtra)
library(factoextra)
library(ggpubr)
library(UpSetR)
library(dplyr)
library(factoextra)
library(cluster)
library(cowplot)
library(viridis) 
library(ggpointdensity)
theme_set(theme_bw())
###################################
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#09194F","#DE3C3C","#1f640a","#490a3d")###################################
###################################
library(ggplot2)
library(gggenes)
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/a/chr1.Region.txt",header=F)
colnames(data) <- c("molecule","gene","start","end")
data$gene <- factor(data$gene, levels = data$gene)
ggplot(data, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow() +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  #geom_gene_label(aes(label = gene), nudge_x = 0.1,  vjust = -0.5, hjust = 0.5, size = 4)  +
  scale_fill_brewer(palette = "RdBu") +
  theme_genes() + ylab(NULL) +
  scale_fill_manual(values= c(colors10,colors10)) + theme(legend.position = 'none')


















