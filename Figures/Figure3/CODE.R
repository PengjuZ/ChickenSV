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
library(ggplotify)
###################################
colors <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure3/a/06-mechanisms.plot.txt",header=F)
colnames(data) <- c("CHR","LOC","ID","Frequency","Length","Mechanisms","Types")
P1 <- gghistogram(data, x = "Frequency", 
               #color = "black",
                fill = "Types",
                alpha = 0.8, bins = 50,
                position = "fill",
                color = "Types", palette = colors,
                ylab = "Proportion", xlab = "Frequency of SVs") + theme(legend.position = "right")

#########
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure3/b/03-Permtest.out.txt",header=F)
colnames(dat) <- c("Regions", "SVs", "Pvalue", "Zscore")
dataM <- dat[,-3]
dataM <- spread(dataM , key=Regions, value=Zscore)
row.names(dataM) <- dataM$SVs
dataM <- as.matrix.data.frame(dataM[,-1]) 
dataM <- dataM[,c(1,2,6,8,5,3,4,7)]
dataM <- dataM[c(4,3,1,2,6,5),]

dataP <- dat[,-4]
dataP <- spread(dataP, key=Regions, value=Pvalue)
row.names(dataP) <- dataP$SVs
dataP <- as.matrix.data.frame(dataP[,-1]) 
dataP <- dataP[,c(1,2,6,8,5,3,4,7)]
dataP <- dataP[c(4,3,1,2,5,6),]

if (!is.null(dataP)){
  sssmt <- dataP < 0.001
  dataP[sssmt] <-'***'
  ssmt <- dataP < 0.01
  dataP[ssmt] <-'**'
  smt <-dataP > 0.01 & dataP <0.05
  dataP[smt] <- '*'
  dataP[!sssmt&!ssmt&!smt]<- ''
} else {
  dataP <- F
}

P2 <- pheatmap(dataM,
              color = colorRamp2(c(-30,0,30), c("#005eaa", "white", "#C00000")),
              cellwidth = 25, cellheight = 25,
              #color = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), c("#C00000","#e41749","#ED7D31","#FEC306","white", "white")),
              #clustering_method = "average",
              legend = T,
              display_numbers = dataP,
              number_color = "black",
              border_color = "black",
              fontsize = 11,
              cluster_cols = F,        
              cluster_rows = F, 
              angle_col = c("45"),
              gaps_col = c(2,3,6),
              gaps_row = c(4),
              # legend_labels = c("10", "20", "30", "40", "title\n"),
) 
#########
library("ggdensity")
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure3/a/06-mechanisms.plot.txt",header=F)
colnames(data) <- c("CHR","LOC","ID","Frequency","Length","Mechanisms","Types")
Mechanisms <- data.frame(Mechanisms = data$Mechanisms,
                         Frequency = as.numeric(data$Frequency),
                         Length = log10(as.numeric(data$Length)),
                         stringsAsFactors = FALSE )
Mechanisms$Mechanisms <- factor(Mechanisms$Mechanisms, levels = c("VNTR","NAHR","STEI","MTEI","NHR"))
P3 <- ggplot(Mechanisms, aes(Length, Frequency, fill = Mechanisms)) +
      geom_hdr(xlim = c(0, 5), ylim = c(0, 1)) +
      scale_fill_brewer(palette = "Set1") + # 调色
      theme_bw() + # 主题
      labs(x="SV Length (log10)",y="Frequency") +
      facet_wrap(~Mechanisms,  nrow = 1) +
      theme(panel.grid.major=element_line(colour=NA),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            panel.grid.minor = element_blank()) 
#########
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure3/d/input.txt",header=T)
data$Length <- data$Length/1000000

dataF <- data[data$Type == "Family",]
P4 <- ggbarplot(dataF, x="ID", y="Length", fill = "ID", color = "black",
          position = position_dodge(0.9),
          palette = colors,
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          sort.by.groups=FALSE,
          x.text.angle=45,
          #rotate=TRUE,
          ylab = "Length (Mb)", xlab = "")  + 
  theme(legend.position = 'none') + theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

dataS <- data[data$Type == "Subfamily",]
P5 <- ggbarplot(dataS, x="ID", y="Length", fill = "ID", color = "black",
          position = position_dodge(0.9),
          palette = colors,
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          sort.by.groups=FALSE,
          x.text.angle=45,
          #rotate=TRUE,
          ylab = "Length (Mb)", xlab = "")  + 
  theme(legend.position = 'none') + theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 12)
  )

dataT <- data[data$Type == "Tandem",]
P6 <- ggbarplot(dataT, x="ID", y="Length", fill = "ID", color = "black",
          position = position_dodge(0.9),
          palette = colors,
          #sort.val = "desc", #上升排序,区别于desc，具体看图演示
          sort.by.groups=FALSE,
          x.text.angle=45,
          #rotate=TRUE,
          fontsize = 11,          
          ylab = "Length (Mb)", xlab = "")  + 
          theme(legend.position = 'none') + theme(
            axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 11),
            axis.title = element_text(size = 12)
          )
#########
M1 <- plot_grid(P6,P4,P5,ncol=3, nrow=1, labels = c('', '', ''),
                rel_widths = c(1) ,
                rel_heights = c(1,1,1) ,align = 'vh'
)

plot_grid(P3,M1, ncol=1, nrow=2, labels = c('', ''),
          rel_widths = c(1),
          rel_heights = c(1,1)
)










###################################
M1 <- plot_grid(P1,as.ggplot(P2),ncol=2, nrow=1, labels = c('', ''),
                rel_widths = c(0.8,1) ,
                rel_heights = c(1)
)


M2 <- plot_grid(M1,P3,ncol=1, nrow=2, labels = c('', ''),
                rel_widths = c(1) ,
                rel_heights = c(1,1),align = 'vh'
)


M2 <- plot_grid(P6,P4,P5,ncol=3, nrow=1, labels = c('', '', ''),
                rel_widths = c(1) ,
                rel_heights = c(1,1,1) ,align = 'vh'
)

plot_grid(M1,P3,M2, ncol=1, nrow=3, labels = c('', '', ''),
          rel_widths = c(1),
          rel_heights = c(1,1,1)
          )








