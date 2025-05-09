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
library(gggenes)
###################################
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#09194F","#DE3C3C","#1f640a","#490a3d")###################################
##############
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/d/06-plot.txt",header=F)
colnames(data) <- c("Location","P1_P2","P3","fdM")
dat <-aggregate(data$fdM, by=list(data$Location,data$P3),mean)
colnames(dat) <- c("Location","P3","avefdM")
dat$Location <- dat$Location/1000
dat <- dat[dat$Location > 8000 & dat$Location < 12000,]

P1 <- ggscatter(dat, x = "Location", y = "avefdM",
                shape = 21, 
                #rotate=TRUE,
                alpha = 0.8,
                #label = "GeneID",
                #repel = T,
                color = "black",
                palette = c("#C00000","#005eaa"), 
                fill= "P3",
                size=2,
) + theme(legend.position = c(0.7, 0.9), legend.justification = c(0, 1)) +
  xlab(paste0("Chromosome (kb)")) + ylab(paste0("Average fdM")) +
  theme(axis.text.x=element_text(vjust=1,size=11,color = "black")) +
  theme(axis.text.y=element_text(vjust=1,size=11,color = "black"))

datN <- data
datN$Location <- datN$Location/1000
datN <- datN[datN$Location > 9704 & datN$Location < 9736,]
datS <-aggregate(datN$fdM, by=list(datN$P1_P2,datN$P3),mean)

P2 <-ggpaired(datS, x="Group.2", y="x",
         color="Group.2", 
         #width = 0.5,
         line.color="gray",
         line.size=0.4, 
         point.size=1,
         palette = colors10) + 
         stat_compare_means(label = "p.format",paired = TRUE) +
         theme (legend.position="none") + 
         xlab(paste0("")) + 
         ylab(paste0("fdM")) 

S1 <- plot_grid(P1,P2, ncol=2, nrow=1 ,rel_widths = c(2,0.8), align = "vh")
plot_grid(NULL,S1, ncol=1, nrow=2 ,rel_heights = c(0.02,1), align = "vh")


