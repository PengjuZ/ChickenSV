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
theme_set(theme_bw())
library(ggplot2)
library("ggsci")
library(ggpubr)
library(UpSetR)
library(dplyr)
library(factoextra)
library(cluster)

###################################
colors <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/a/04-stats.txt",header=F)
colnames(data) <- c("Count","Genotye","Strains")
#data <- data[data$Strains != "BRB" & data$Strains != "BRA" & data$Strains != "BL" & data$Strains != "WL" & data$Strains != "RJFt" & data$Strains != "RJFi" & data$Strains != "LinD",]
data <- data[data$Strains != "RJFi",]
data$Genotye <- factor(data$Genotye, levels = c("1|1","0|1","0|0"))
YES <- data[data$Genotye == "1|1",]
data$Strains <- factor(data$Strains, levels = YES[order(YES$Count),]$Strains )
data<- data %>% group_by(Strains) %>% mutate(Percent = Count/sum(Count))

ggbarplot(data, x="Strains", 
          y="Percent", 
          alpha = 0.8,
          fill = "Genotye", 
          color = "black",
          #position = position_dodge(0.9),
          palette = c("#C00000","#ED7D31","#FFC000","#005eaa"),
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          #sort.by.groups=FALSE,
          #x.text.angle=45,
          rotate=TRUE,
          width =1,
          ylab = "Proportion", xlab = "") +
 theme(legend.position = "right") 





