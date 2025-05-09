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
library(cowplot)
library(forcats)
theme_set(theme_bw())
###################################
colors <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
###################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure1/b/Nx.plot.txt", header =F)
colnames(dat) <- c("Bin", "Length", "Genome")
dat$Genomes <- "Others"
dat$Types <- "Others"
dat[dat$Genome == "G01.C",]$Genomes <- "Huxu"
dat[dat$Genome == "G01.C",]$Types <- "Contigs"
dat[dat$Genome == "G01.S",]$Genomes <- "Huxu"
dat[dat$Genome == "G01.S",]$Types <- "Scaffolds"
dat[dat$Genome == "G02.C",]$Genomes <- "GRCg7b"
dat[dat$Genome == "G02.C",]$Types <- "Contigs"
dat[dat$Genome == "G02.S",]$Genomes <- "GRCg7b"
dat[dat$Genome == "G02.S",]$Types <- "Scaffolds"
dat[dat$Genome == "G03.C",]$Genomes <- "GRCg7w"
dat[dat$Genome == "G03.C",]$Types <- "Contigs"
dat[dat$Genome == "G03.S",]$Genomes <- "GRCg7w"
dat[dat$Genome == "G03.S",]$Types <- "Scaffolds"
dat[dat$Genome == "G04.C",]$Genomes <- "WChap1"
dat[dat$Genome == "G04.C",]$Types <- "Contigs"
dat[dat$Genome == "G04.S",]$Genomes <- "WChap1"
dat[dat$Genome == "G04.S",]$Types <- "Scaffolds"
dat[dat$Genome == "G05.C",]$Genomes <- "WChap2"
dat[dat$Genome == "G05.C",]$Types <- "Contigs"
dat[dat$Genome == "G05.S",]$Genomes <- "WChap2"
dat[dat$Genome == "G05.S",]$Types <- "Scaffolds"
dat$Types <- factor(dat$Types,levels =c("Scaffolds","Contigs"))

P1 <- ggplot()+ geom_line(data=dat,aes(x=Bin,y=Length,color=Genomes,linetype=Types),lwd=0.5)+  
        guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
        scale_color_manual(values = colors)+   #color对应的颜色
        labs(x='N(x) %',y='Length (bp)')+ 
         #geom_vline(xintercept = 50,linetype="dotted") +  
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  #geom_vline(xintercept = 90,linetype="dotted") + 
  scale_y_continuous(
        trans = log10_trans(),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x))) 

###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure1/c/Busco.list.plot.txt",header=F)
colnames(data) <- c("genomes","Methods","MethodsC","ID","Counts","Types")
data <- aggregate(data[,5],data[,c(3,6)],mean)
colnames(data) <- c("genomes","Types","Counts")
data$genomes <- factor(data$genomes, levels = c("Abinitio","RNA-seq","Proteins","EVM","PASA","Final","GRCg7b","GRCg7w"))
P2 <- ggbarplot(data, x = "genomes", 
          y="Counts",
          #color = c("black"),
          fill = "Types",
          #rotate=TRUE,
          alpha = 0.8, 
          x.text.angle=30,
          #position = "stack",
          color = "Types", palette = colors,
          ylab = "Count", xlab = "") 
###################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure1/d/00-Glist.Nx.plot.txt", header =F)
colnames(dat) <- c("Bin", "Length", "Genome", "Annotations")
dat <- aggregate(dat[,2],dat[,c(1,4)],mean)
colnames(dat) <- c("Bin", "Annotations", "Length")
dat$Annotations <- factor(dat$Annotations, levels = c("Abinitio","RNA-seq","Proteins","EVM","PASA","Final","GRCg7b","GRCg7w"))
dat$Genomes <- "WChap"
dat[dat$Annotations == "GRCg7b",]$Genomes <- "GRCg7b"
dat[dat$Annotations == "GRCg7w",]$Genomes <- "GRCg7w"
dat[dat$Annotations == "GRCg7b",]$Annotations <- "Final"
dat[dat$Annotations == "GRCg7w",]$Annotations <- "Final"
dat$Genomes <- factor(dat$Genomes, levels = c("WChap","GRCg7b","GRCg7w"))

P3 <- ggplot()+ geom_line(data=dat,aes(x=Bin,y=Length,color=Annotations,linetype=Genomes),lwd=0.5) + 
     guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2)) +
     scale_color_manual(values = colors)+   #color对应的颜色
     labs(x='N(x) %',y='Length (bp)')+
     theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())+
    scale_y_continuous(
    trans = log10_trans(),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x)))


###################################

M1 <- plot_grid(P2,P1,ncol=2, nrow=1, labels = c('', ''),
                rel_widths = c(2, 1.5) ,
                rel_heights = c(1), scale = c(1, .9, .9, .8)
                
)
M2 <- plot_grid(P2,P3,ncol=2, nrow=1, labels = c('', ''),
                rel_widths = c(2, 1.5) ,
                rel_heights = c(1), scale = c(1, .9, .9, .8)
)
plot_grid(M1,M2,ncol=1, nrow=2)

