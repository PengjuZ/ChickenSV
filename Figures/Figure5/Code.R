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
library(viridis) # 使用viridis提供的翠绿色标度：scale_fill_viridis()
library(ggpointdensity)
library(ggplotify)
library(cowplot)
library(forcats)
library(tidyr)
theme_set(theme_bw())
###################################
col <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
###################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure5/a/Permtest.peak.out.plot.txt", header =F)
colnames(dat) <- c("Histone", "Tissue", "SV", "Zscore")
dat <- dat[dat$SV!="SV.1" & dat$SV!="SV.2" & dat$SV!="SV.3" & dat$SV!="SV.4" & dat$SV!="SV.5" & dat$SV!="SV.6" & dat$SV!="SV.7" & dat$SV!="SV.8" & dat$SV!="SV.9" & dat$SV!="SV.10",]
dat$SV <- factor(dat$SV,
                 levels = rev(c("SV.withnold","SV.withld","SV.Strains.specific","SV.Strains.shared","SV.rare","SV.low-frequency","SV.common","SV.high-frequency")))
P1 <- ggboxplot(dat,
          x="SV",
          y="Zscore", 
          #color = "Tissues", 
          color = "Histone", 
          #shape = "Tissues",
          #shape = "3",
          palette = col, 
          sort.by.groups=T,
          rotate=TRUE,
          #facet.by = "Histone", 
          ylab = "Zscore", xlab = "",
          add.params=list(size=2, shape = 4),
          add = "jitter")+ 
  geom_hline(yintercept = 0,linetype="dotted") + 
  facet_grid(.~Histone, scales = "free_x") +
  theme (axis.text = element_text (size = 9)) +
  theme(legend.position = 'none') + 
  theme(panel.spacing = unit(0.3,"lines")) +
  theme(strip.text.x = element_text(size = 9)) +
  theme(plot.margin =margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))
####
#datS <- dat[dat$SV == "SV.Strains.specific" & dat$Histone != "ATAC" & dat$Histone != "DNaseSeq" & dat$Histone != "CTCF",]
#datS <- dat[dat$SV == "SV.Strains.specific" & dat$Histone == "CTCF",]
datS <- dat[dat$Histone == "H3K27me3",]
datO <- datS[datS$SV == "SV.common",]
datS$Tissue <- factor(datS$Tissue,levels = datO[order(datO$Zscore),]$Tissue)
datS <- datS[,-1]
Heat <- spread(datS, key=Tissue, value=Zscore)
row.names(Heat) <- Heat$SV
Heat <- as.matrix.data.frame(Heat[,-1]) 
P2 <- pheatmap(Heat,
         color = colorRamp2(c(-10,0,10), c("#005eaa", "white", "#C00000")),
         cellwidth = 15, cellheight = 15,
         #color = colorRamp2(c(1,0.8,0.6,0.4,0.2,0), c("#C00000","#e41749","#ED7D31","#FEC306","white", "white")),
         #clustering_method = "average",
         legend = T,
         # display_numbers = dataP,
         number_color = "black",
         border_color = "black",
         fontsize = 9,
         cluster_cols = T,        
         cluster_rows = F, 
         angle_col = c("45"),
         # gaps_col = c(7,17,23),
         gaps_row = c(4,6,8),
         # legend_labels = c("10", "20", "30", "40", "title\n"),
) 
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure5/a/Permtest.state.out.txt", header =F)
colnames(data) <- c("Tissue", "StateID", "State",  "SV", "Pvalue", "Zscore")
data <- data[data$SV!="SV.1" & data$SV!="SV.2" & data$SV!="SV.3" & data$SV!="SV.4" & data$SV!="SV.5" & data$SV!="SV.6" & data$SV!="SV.7" & data$SV!="SV.8" & data$SV!="SV.9" & data$SV!="SV.10",]
#data <- data[data$SV == "SV.common" | data$SV == "SV.high-frequency" | data$SV == "SV.low-frequency" | data$SV == "SV.rare",]
data$State <- factor(data$State,
                     levels = c("TssA","TssAHet","TxFlnk","TxFlnkWk","TxFlnkHet","EnhA","EnhAMe","EnhAWk","EnhAHet","EnhPois","ATAC_Is","TssBiv","Repr","ReprWk","Qui"))
data$SV <- factor(data$SV,
                  levels = rev(c("SV.withnold","SV.withld","SV.Strains.specific","SV.Strains.shared","SV.rare","SV.low-frequency","SV.common","SV.high-frequency")))
#data <- data[data$Tissue != "Kidney",]
datas <- aggregate(data$Zscore, list(data$SV,data$State), sd)
datas$M <- aggregate(data$Zscore, list(data$SV,data$State), mean)$x
datas$CV <- abs(datas$x/datas$M )
datas <- datas[,c(1,2,5)]
colnames(datas) <- c("SV", "State",  "CV")
P3 <- ggline(datas,x="State",
        plot_type = c("l"),
        palette = col,
        #rotate=TRUE,
        color = "SV",
        #linetype = "Types",
        #shape = "",
        y="CV",
        size=1.2,
        point.size = 0.01) + 
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) +
        geom_hline(yintercept = 1,linetype="dotted") + theme(legend.position = 'right')
######
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure5/a/Permtest.state.out.txt", header =F)
colnames(data) <- c("Tissue", "StateID", "State",  "SV", "Pvalue", "Zscore")
data <- data[data$SV!="SV.1" & data$SV!="SV.2" & data$SV!="SV.3" & data$SV!="SV.4" & data$SV!="SV.5" & data$SV!="SV.6" & data$SV!="SV.7" & data$SV!="SV.8" & data$SV!="SV.9" & data$SV!="SV.10",]
#"SV.withnold","SV.withld","SV.Strains.specific","SV.Strains.shared","SV.rare","SV.low-frequency","SV.common","SV.high-frequency"      

dataP <- data[data$StateID == "E2" & data$SV == "SV.high-frequency",]
P4 <- ggbarplot(dataP, x="Tissue", y="Zscore", fill = "#C00000", color = "black",
          position = position_dodge(0.9),
          palette = col,
          sort.val = "asc", #上升排序,区别于desc，具体看图演示
          #sort.by.groups=FALSE,
          x.text.angle=0,
          rotate=TRUE,
          ylab = "Zscore", xlab = "")  +  theme(legend.position = 'none')

dataP <- data[data$StateID == "E6" & data$SV == "SV.Strains.specific",]
P5 <- ggbarplot(dataP, x="Tissue", y="Zscore", fill = "#005eaa", color = "black",
          position = position_dodge(0.9),
          palette = col,
          sort.val = "asc", #上升排序,区别于desc，具体看图演示
          #sort.by.groups=FALSE,
          x.text.angle=0,
          rotate=TRUE,
          ylab = "Zscore", xlab = "")  +  theme(legend.position = 'none')

dataP <- data[data$State == "E12" & data$SV == "SV.low-frequency",]
P6 <- ggbarplot(dataP, x="Tissue", y="Zscore", fill = "#ED7D31", color = "black",
          position = position_dodge(0.9),
          palette = col,
          sort.val = "asc", #上升排序,区别于desc，具体看图演示
          #sort.by.groups=FALSE,
          x.text.angle=0,
          rotate=TRUE,
          ylab = "Zscore", xlab = "")  +  theme(legend.position = 'none')

M <- plot_grid(P4,P5,ncol =2, align = 'hv', rel_widths =  c(1,1))
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure5/c/07-heatmap.txt",header=F)
row.names(data) <- data$V1
data <- data[,-1]
colnames(data) <- c('LinD','RJFt','BL','WL','WZ','WY','WJ','WB','TN','LD','LA','CU3','CK','CB','BRB','BRA','RJFi')
data <- data[,c(3,4,16,15,1,10,11,12,13,14,9,7,6,5,8,2,17)]
P7 <- densityHeatmap(data,
               title = "",
               ylab = "Frequency",
               cluster_columns = TRUE, 
               clustering_distance_columns = "ks")
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure5/d/03-AF.plotc.txt",header=F)
colnames(data) <- c("Rank", "AF", "Type", "ID", "GeneID")
data <- data[data$Type == "RJFs" | data$Type == "LRs",]
P8 <- ggscatter(data, x = "Rank", y = "AF",
               color = "Type", 
               palette = c("#C00000","#005eaa","#ED7D31","#1E925E","#791E94"), 
               fill= "Type",
               ellipse = TRUE,
               #rotate=TRUE,
               label = "GeneID",
               repel = T,
               font.label = c(7), 
               size=2,
               alpha=0.1,
               ellipse.level = 0,
               ellipse.alpha=0) +
  xlab(paste0("Rank")) + ylab(paste0("AF")) + 
  theme(legend.position = "none") + 
  facet_grid(.~Type, scales = "free_x") 


P1
          
plot_grid(as.ggplot(P2),P3, nrow=2,align = 'vh',rel_heights = c(1,1))




#M1 <-plot_grid(P1,NULL, as.ggplot(P2), nrow=3,align = 'vh',rel_heights = c(1,0.1,1))
M1 <- plot_grid(P1,NULL,P1, nrow=3,align = 'vh',rel_heights = c(1,0.1,1))
M2 <- ggarrange(P4,P5,P6,ncol = 3, common.legend = T)


M2 <- plot_grid(P4,P5,P6,ncol=3, nrow=1,align = 'vh', rel_widths =  c(1,1,1))
M3 <- plot_grid(P3,M2,ncol=1, nrow=2,align = 'vh', rel_heights = c(1.2,2))
      plot_grid(M1,M3,ncol=1, nrow=2,align = 'vh', rel_widths =  c(1,1))




M1 <-plot_grid(NULL,as.ggplot(P2), ncol=2, nrow=1,align = 'vh', rel_widths =  c(0.2,1))





M1 <- plot_grid(P1,P2,P3, ncol=3, nrow=1,align = 'vh', rel_widths =  c(1,1,2))
M2 <- plot_grid(P4,P5, ncol=2, nrow=1,align = 'vh', rel_widths =  c(1,1))

plot_grid(M1,M2, 
          ncol=1, nrow=2, 
          align = 'vh', 
          rel_widths =  c(1,1))


M0 <- plot_grid(NULL,as.ggplot(P2),ncol=2, nrow=1, labels = c('', ''),
                rel_widths = c(0.03,1) ,
                rel_heights = c(1)
)
M1 <- plot_grid(NULL,P1,M0,ncol=1, nrow=3, labels = c('', ''),
                rel_widths = c(1) ,
                rel_heights = c(0.05,1,1)
)
M2 <- plot_grid(NULL,P3,NULL,P4,ncol=1, nrow=4, labels = c('', ''),
                rel_widths = c(1) ,
                rel_heights = c(0.1,1,-0.1,1),align = 'vh'
)


M3 <- plot_grid(M2,NULL,P5,ncol=1, nrow=3, labels = c('', ''),
                rel_widths = c(1) ,
                rel_heights = c(1.3,-0.1,1)
)

plot_grid(M1,NULL,M3, ncol=3, nrow=1 ,rel_widths = c(1.5,0,1))


