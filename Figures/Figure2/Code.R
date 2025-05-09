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
library(ggpubr)
require(cowplot)
theme_set(theme_bw())
#################################################################
colors <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")
#################################################################
##########
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure2/b/02.Length.txt",header=F)
colnames(data) <- c("Chr","Loc","Number","length")

P1 <- gghistogram(data, x = "Number", fill = "#C00000", color = "black", bins = 30,alpha = 0.8) + 
  labs(x = "Number of assemblies", y = "Count")

dataL <- as.data.frame(data[data$length > -10000 & data$length < 10000,])
dataL[dataL$length > -10000 & dataL$length < 0,]$length <- 0 - dataL[dataL$length > -10000 & dataL$length < 0,]$length
dataL$length <- dataL$length/1000

P2 <- gghistogram(dataL, x = "length", fill = "#005eaa",bins = 30,alpha = 0.8) + 
  labs(x = "Change in sequence length (Kb)", y = "Count")
##########
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure2/c/01-Stats.txt",header=T)
colnames(data) <- c("Types","Type","Value")
P3 <- ggbarplot(data, x="Type", y="Value", fill = "Types", color = "black",
                alpha = 0.8,
          position = position_dodge(0.9),
          palette = colors,
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          sort.by.groups=FALSE,
          #x.text.angle=45,
          #rotate=TRUE,
          ylab = "Count", xlab = "")  +
          scale_y_continuous(sec.axis = sec_axis(~./1000, name = "Length (Mb)"))
##########
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure2/d/Merge.stats.txt", header =F)
colnames(dat) <- c("Ind", "Types", "SVTypes", "Count")
dat <- dat[dat$Count < 20000,]
dat$SVTypes <- factor(dat$SVTypes, levels = c("All SVs","Deletion","Insertion"))
P4 <- ggboxplot(dat,
          x="SVTypes",
          y="Count", 
          #color = "Tissues", 
          color = "Types", 
          #shape = "Tissues",
          #shape = "3",
          alpha = 0.8,
          palette = colors, 
          sort.by.groups=T,
          #rotate=TRUE,
          #facet.by = "Histone", 
          ylab = "Count of large SVs", xlab = "",
          add.params=list(size=1, shape = 4),
          add = "jitter")+ 
          theme(strip.text.x = element_text(size = 12))
##########
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure2/e/Whole.frq.strat.stats.txt",header=F)
colnames(data) <- c("Count", "Strains")
data$Strains <- factor(data$Strains, levels = data$Strains)
P5 <- ggbarplot(data, x="Strains", y="Count", fill = "Strains", color = "black",
          position = position_dodge(0.9),
          alpha = 0.8,
          palette = colors,
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          sort.by.groups=FALSE,
          #x.text.angle=45,
          #rotate=TRUE,
          ylab = "Count", xlab = "Number of strains")  + 
          theme(legend.position = 'none')
##########
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure2/f/POPClass.chr.class.txt", header =F)
colnames(dat) <- c("Count", "Strains", "Chr", "Size", "Types", "SVdensity")
P6 <- ggboxplot(dat,
          x="Strains",
          y="SVdensity", 
          alpha = 0.8,
          #color = "Tissues", 
          color = "Types", 
          #shape = "Tissues",
          #shape = "3",
          palette = colors, 
          sort.by.groups=T,
          #rotate=TRUE,
          #facet.by = "Histone", 
          ylab = "SV density (times/Mb)", xlab = "Number of strains",
          add.params=list(size=2, shape = 4),
          add = "jitter")


M1 <- plot_grid(P1,P3,P4,ncol=3, nrow=1, labels = c('b', 'c', 'd'),
                rel_widths = c(1.1, 1, 1.2) ,
                rel_heights = c(1)
                )
M2 <- plot_grid(P5,P6,ncol=2, nrow=1,labels = c('e', 'f'),
                rel_widths = c(1, 2) ,
                rel_heights = c(1)                
                )

plot_grid(M1,M2,ncol=1, nrow=2)







