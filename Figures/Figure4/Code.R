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
theme_set(theme_bw())
###################################
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure4/c-f/MergeLD.plot.txt",header=T)
colnames(data) <- c("ID","SV","Length","Frequency","Type")
data <- data[data$Length < 10000,]
P1 <- ggboxplot(data, x = "Type", y = "Frequency", width = 0.6, 
          color = "black",#轮廓颜色
          fill="Type",#填充
          #alpha = 0.8,
          palette =c("#C00000","#005eaa"),#分组着色
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          #size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          #legend = "right"
          ) + theme (legend.position="none") +
             stat_compare_means(method = "t.test", label.y = 0.75) + 
             coord_cartesian(ylim = c(0, 0.75)) 


P2 <- ggboxplot(data, x = "Type", y = "Length", width = 0.6, 
          color = "black",#轮廓颜色
          fill="Type",#填充
          #alpha = 0.8,
          palette =c("#C00000","#005eaa"),#分组着色
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          #size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          ylab = "Length (bp)",
          #legend = "right"
          )  + theme (legend.position="none") + 
          stat_compare_means(method = "t.test", label.y = 2000) + 
          coord_cartesian(ylim = c(0, 2000)) 


###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure4/e/SNP.SV.IBS.plot.txt",header=F)
colnames(data) <- c("Breeds","Strains","IBS","NonLD","LD","Ratio")
data <- data[data$Strains != "RJFi",]
P3 <- ggscatter(data, x = "IBS", y = "Ratio",
          color = "black",
          shape = 21, 
          palette = colors10 , 
          fill= "Strains",
          size=3,
          alpha = 1,
          label = "Breeds",
          repel = T,
          add = "reg.line",
          conf.int = TRUE,
          add.params = list(color = "#A6A6A6",fill = "#D9D9D9"),
          xlab="SNP-based IBS",
          ylab="SNP-tagged / independent SVs",         
) + stat_cor(method = "pearson",
             label.x = 0.855, 
             label.y = 1.38)  +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + 
  theme(legend.position = 'none') 
  #theme(plot.margin =margin(t = 0, r = 10, b = 0, l = 0, unit = "pt"))
###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure4/g/03.plot.txt",header=F)
colnames(data) <- c("Strain1","Strain2","SNP","LD","Non-LD")
data <- data[data$Strain1 != "RJFi" & data$Strain2 != "RJFi",]
data$SNP <- 1- data$SNP
data$LD <- 1- data$LD
data[,5] <- 1- data[,5]
#data$Non-LD <- 1- data$Non-LD
data <- reshape2::melt(data,
                       measure.vars = c('LD','Non-LD'),
                       variable.name = "Type",
                       value.name = "IBS")

P4 <- ggscatter(data, x = "SNP", y = "IBS",
          color = "black",
          shape = 21,
          palette = colors10 , 
          fill= "Type",
          #ellipse = TRUE,
          size=2,
          ellipse.level = 0,
          alpha = 0.8,
          xlab="SNP-based IBS",
          ylab="SV-based IBS",         
) +  geom_smooth(aes(color = Type),linetype="dashed", method = lm , se = FALSE, fullrange = TRUE)+
  scale_color_manual(values = c("#C00000","#005eaa")) +
  ggpubr::stat_cor(aes(color = Type), label.x = 0.76, label.y = c(0.82,0.83)) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())  + theme(legend.position = 'none') 
#+ theme(legend.position = 'right',legend.direction = 'vertical',legend.key.size = unit(18,"point"))
#####################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure4/g/03.plot.txt",header=F)
colnames(data) <- c("Strain1","Strain2","SNP","LD","NonLD")
data <- data[data$Strain1 != "RJFi" & data$Strain2 != "RJFi",]
data$SNP <- 1- data$SNP
data$LD <- 1- data$LD
data$NonLD <- 1- data$NonLD
data$Dm <- data$NonLD - data$SNP
dataM <- data[data$Strain1 == "RJFt",]
dataM$Strain <- "WC"
dataM[dataM$Strain2 == "BL",]$Strain <- "BL"
dataM[dataM$Strain2 == "WL",]$Strain <- "WL"
dataM[dataM$Strain2 == "BRA",]$Strain <- "BRA"
dataM[dataM$Strain2 == "BRB",]$Strain <- "BRB"
dataM[dataM$Strain2 == "LinD",]$Strain <- "LinD"

P5 <- ggdotchart(dataM, x = "Strain2", y = "Dm",
           color = "Strain",
          # shape = 21, 
           palette = colors10 , 
           fill= "Strain",
           rotate = TRUE,
           sorting = "descending",
           size=3,
           alpha = 0.8,
           xlab="Population",
           ylab="IBS value with RJFt (SV - Snp)",
) + theme(legend.position = 'none') 


#####################################
plot_grid(P1,P2,P3,P4, ncol=4, nrow=1,align = 'vh', rel_widths =  c(0.84,0.84,1,1))


M1 <- plot_grid(P1,P2, ncol=2, nrow=1,align = 'vh', rel_widths =  c(1,1))
plot_grid(M1,P3,P4, ncol=1, nrow=3,align = 'vh', rel_heights =  c(1,1,1))
plot_grid(NULL,M1,M2, 
          ncol=1, nrow=3, 
          align = 'vh', 
          rel_heights =  c(0.01,1,1))




