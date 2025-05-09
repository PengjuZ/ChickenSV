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
##############library(ggplot2)
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure6/b/chr1.Region.txt",header=F)
colnames(data) <- c("molecule","gene","start","end")
data$gene <- factor(data$gene, levels = data$gene)
P1 <- ggplot(data, aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
       geom_gene_arrow() +
       facet_wrap(~ molecule, scales = "free", ncol = 1) +
       #geom_gene_label(aes(label = gene), nudge_x = 0.1,  vjust = -0.5, hjust = 0.5, size = 4)  +
       #scale_fill_brewer(palette = "RdBu") +
       theme_genes() + ylab(NULL) + theme(axis.text.x = element_blank()) +
       scale_fill_manual(values= c(colors10,colors10)) + theme(legend.position = 'none')
#####################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure6/c/04-Pggb7614.gwas.txt",header=F)
colnames(data) <- c("ID","Loc","Pvalue","R2","R2Class","GeneID")
data$P <- -log10(data$Pvalue)
data$Locs <- data$Loc/1000
data <- data[data$Locs > 8800 & data$Locs < 10300,]
P2 <- ggscatter(data, x = "Locs", y = "P",
          shape = 21, 
          #rotate=TRUE,
          #alpha = 0.8,
          #label = "GeneID",
          #repel = T,
          color = "black",
          palette = c("#005eaa","#1E925E","#FFC000","#ED7D31","#C00000"), 
          fill= "R2Class",
          size=1.5,
  ) + theme(legend.position = "none") +
  xlab(paste0("")) + ylab(paste0("-log(p-value)")) +
  theme(axis.text.x=element_text(vjust=1,size=11,color = "black")) +
  theme(axis.text.y=element_text(vjust=1,size=11,color = "black"))

#####################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure6/d/03-nIBD.plot.txt", header =F)
colnames(dat) <- c("chr", "start", "end", "I1", "I2", "N", "C", "IBD", "class")
dat$CLASS <- "Colored feather"
dat[dat$class == "CB_LD" | dat$class == "CB_WB" | dat$class == "LD_WB",]$CLASS <- "White feather"
dat$MID=(dat$start + dat$end)/2000
dat$nIBD <- dat$IBD/1600
dat <- dat[dat$MID > 8800 & dat$MID < 10300,]
P3 <- ggplot()+ geom_line(data=dat,aes(x=MID,y=nIBD,color=CLASS),alpha = 0.8,lwd=0.5)+  
      scale_color_manual(values = c("#005eaa","#C00000"))+   #color对应的颜色
      labs(x='',y='nIBD')+ 
      theme(legend.position = c(0.7, 1.1), legend.justification = c(0, 1))+ 
      theme(axis.line = element_line(color='black'),
        legend.title = element_blank(),    
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
      theme(axis.text.x=element_text(vjust=1,size=11,color = "black")) +
      theme(axis.text.y=element_text(vjust=1,size=11,color = "black"))

#####################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure6/e/07-XP-EHH.window.plot.txt",header=F)
colnames(data) <- c("Chr","Start","End","xpehh")
data$Location <- (data$Start + data$End)/2000
data$Sign <- "Un"
data[data$xpehh > 1,]$Sign <- "Sign"
P4 <- ggscatter(data, x = "Location", y = "xpehh",
          shape = 21, 
          #rotate=TRUE,
          alpha = 0.8,
          #label = "GeneID",
          #repel = T,
          color = "black",
          palette = c("#C00000","#005eaa"), 
          fill= "Sign",
          size=2,
     ) + theme(legend.position = "none") +
      #geom_text_repel(aes(label = GeneID),size=3.8,color="black",direction="both",min.segment.length = 0.05,segment.alpha=0.6,max.overlaps =100,nudge_x = 6,nudge_y=4) +
     xlab(paste0("Chromosome (kb)")) + ylab(paste0("XP-EHH")) +
     theme(axis.text.x=element_text(vjust=1,size=11,color = "black")) +
     theme(axis.text.y=element_text(vjust=1,size=11,color = "black"))

#####################
S1 <-plot_grid(NULL,P1,NULL, ncol=3, nrow=1 ,rel_widths = c(-0.065,1,0.026), align = "vh")
M1 <-plot_grid(S1,NULL,P2,NULL,P3,NULL,P4, ncol=1, nrow=7 ,rel_heights = c(1,-0.3,1.2,-0.1,1,-0.1,1), align = "vh")

ggsave("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure6/myplot.pdf", M1, width = 16, height = 19, units = "cm")
