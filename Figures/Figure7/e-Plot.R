library(devtools)
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
library(rstatix)
theme_set(theme_bw())
########################################
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#09194F","#DE3C3C","#1f640a","#490a3d")###################################
########################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/e/data.txt")

dat$Color <- "Colored"
#dat[dat$POP == "CB",]$Color <- "CB-White"
#dat[dat$POP == "LD",]$Color <- "LD-White"
#dat[dat$POP == "WB",]$Color <- "WB-White"
dat[dat$POP == "CB" | dat$POP == "LD" | dat$POP == "WB",]$Color <- "White"

P1 <- ggline(dat,x="Date",y="Weights",add = "mean_ci",color = "Color",
       palette = colors10,
       size = 1,
       ylab = "Body weight (g)", xlab = "Weeks of age"
       )


datE <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/e/dataE.txt")
datE <-datE[datE$Tissue != "Crop",]
P2 <- ggbarplot(datE[rev(order(datE$GRM5)),][10:1,], x="Tissue", y="GRM5", fill = "Tissue", color = "black",
                position = position_dodge(0.9),
                palette = rev(colors10),
                sort.val = "asc", #上升排序,区别于desc，具体看图演示
                #sort.by.groups=FALSE,
                x.text.angle=0,
                rotate=TRUE,
                ylab = "Average TPM", xlab = "")  +  
                theme(legend.position = 'none') + theme(axis.text = element_text(size = 10))

P3 <- ggbarplot(datE[rev(order(datE$TYR)),][10:1,], x="Tissue", y="TYR", fill = "Tissue", color = "black",
                position = position_dodge(0.9),
                palette = rev(colors10),
                sort.val = "asc", #上升排序,区别于desc，具体看图演示
                #sort.by.groups=FALSE,
                x.text.angle=0,
                rotate=TRUE,
                ylab = "Average TPM", xlab = "")  +  
                theme(legend.position = 'none') + theme(axis.text = element_text(size = 10))

P4 <- ggbarplot(datE[rev(order(datE$LOC107052320)),][10:1,], x="Tissue", y="LOC107052320", fill = "Tissue", color = "black",
                position = position_dodge(0.9),
                palette = rev(colors10),
                sort.val = "asc", #上升排序,区别于desc，具体看图演示
                #sort.by.groups=FALSE,
                x.text.angle=0,
                rotate=TRUE,
                ylab = "Average TPM", xlab = "")  +  
                theme(legend.position = 'none') + theme(axis.text = element_text(size = 10))

P5 <- ggbarplot(datE[rev(order(datE$NOX4)),][10:1,], x="Tissue", y="NOX4", fill = "Tissue", color = "black",
                position = position_dodge(0.9),
                palette = rev(colors10),
                sort.val = "asc", #上升排序,区别于desc，具体看图演示
                #sort.by.groups=FALSE,
                x.text.angle=0,
                rotate=TRUE,
                ylab = "Average TPM", xlab = "")  +  
                theme(legend.position = 'none') + theme(axis.text = element_text(size = 10))





S1 <- plot_grid(P2,P3, ncol=2, nrow=1 ,rel_widths = c(1,1), align = "vh")
S2 <- plot_grid(P4,P5, ncol=2, nrow=1 ,rel_widths = c(1,1), align = "vh")
S3 <- plot_grid(S1,S2, ncol=1, nrow=2 ,rel_heights = c(1,1), align = "vh")
S4 <- plot_grid(P1,NULL, ncol=1, nrow=2 ,rel_heights = c(1,0.05), align = "vh")
plot_grid(S4,S3, ncol=2, nrow=1 ,rel_widths = c(1.2,3), align = "vh")




plot_grid(P1,NULL,P2,NULL,P3,NULL,P4,NULL,P5 ,ncol=9, nrow=1 ,rel_widths = c(1.5,-0.1,1,-0.03,1,-0.05,1,-0.05,1), align = "vh")


















ggboxplot(dat,
                x="Pops",
                y="P17", 
                #color = "Tissues", 
                color = "Feather", 
                #shape = "Tissues",
                #shape = "6",
                #palette = col, 
                sort.by.groups=T,
                #rotate=TRUE,
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





























pheatmap(data, 
         color = c("white","white","#f9ae08","#f9802d","#C00000" ),
         #color = colorRampPalette(brewer.pal(n = 5, name = "YlOrRd"))(5), 
         # color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(10),
         #color = colorRampPalette(brewer.pal(n = 9, name = "OrRd"))(10),
         clustering_method = "average",
         border_color = "#EEEEEE",
         fontsize = 11,
         cluster_cols = T,        
         cluster_rows = T, 
         display_numbers = TRUE,  
         number_format = "%.3f",
         show_rownames = TRUE,
         show_colnames = TRUE,
         #annotation_col = annotation_col, 
         #annotation_colors = ann_colors,
         #gaps_col = c(1, 1, 6, 47),
         #gaps_row = c(6,22,27),
         angle_col = c("0")
)
########################################Pi
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/ab/DataPi.txt")
row.names(dat) <- dat$Breeds
ggbarplot(dat, x="Pops", y="Pi", 
          orientation = "horiz", 
                  fill = "Pops", 
                  #color = "Pops",
                  palette = rep(my36colors,4), #杂志nature的配色
                  sort.val = "desc", #下降排序
                  sort.by.groups= F, #不按组排序
                  x.text.angle=0,
                  ylab = "Nucleotide diversity (π)", xlab = "") +   theme(text = element_text(size=11)) + theme(legend.position = "none")
########################################
library(rstatix)
options(scipen = 999)
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/ab/Pdata.txt")
Sex1 <- dat %>% t_test(P1 ~ Sex, p.adjust.method = "bonferroni")
Sex2 <- dat %>% t_test(P2 ~ Sex, p.adjust.method = "bonferroni")
Sex3 <- dat %>% t_test(P3 ~ Sex, p.adjust.method = "bonferroni")
Sex4 <- dat %>% t_test(P4 ~ Sex, p.adjust.method = "bonferroni")
Sex5 <- dat %>% t_test(P5 ~ Sex, p.adjust.method = "bonferroni")
Sex6 <- dat %>% t_test(P6 ~ Sex, p.adjust.method = "bonferroni")
Sex7 <- dat %>% t_test(P7 ~ Sex, p.adjust.method = "bonferroni")
Sex8 <- dat %>% t_test(P8 ~ Sex, p.adjust.method = "bonferroni")
Sex9 <- dat %>% t_test(P9 ~ Sex, p.adjust.method = "bonferroni")
Sex10 <- dat %>% t_test(P10 ~ Sex, p.adjust.method = "bonferroni")
Sex11 <- dat %>% t_test(P11 ~ Sex, p.adjust.method = "bonferroni")
Sex12 <- dat %>% t_test(P12 ~ Sex, p.adjust.method = "bonferroni")
rbind(Sex1,Sex2,Sex3,Sex4,Sex5,Sex6,Sex7,Sex8,Sex9,Sex10,Sex11,Sex12)
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures/10-Feather/Pdata.txt")
dat <- dat[dat$Sex == "F",]
stat.test <- dat %>% t_test(P1 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P1 <- pheatmap(-log10(PLOTF), 
         color = c("white","white","#f9ae08","#f9802d","#C00000" ),
         clustering_method = "average", border_color = "gray",
         cellwidth = 10, cellheight = 10, fontsize = 10,
         cluster_cols = F, cluster_rows = F, 
         show_rownames = TRUE, show_colnames = TRUE,
         angle_col = c("90"), column_title = "P1",
         name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P2 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P2 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P2",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P3 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P3 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P3",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P4 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P4 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P4",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )



P1+P2+P3+P4
P5+P6+P7+P8
P9+P10+P11+P12




##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/PdataF.txt")
stat.test <- dat %>% t_test(P17 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P17 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P17",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/PdataF.txt")
stat.test <- dat %>% t_test(P18 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P18 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P18",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/PdataF.txt")
stat.test <- dat %>% t_test(P19 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P19 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P19",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/PdataF.txt")
stat.test <- dat %>% t_test(P20 ~ Pops, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("CB", "CB", 1)
row2<- c("WZ", "WZ", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:11]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=10,ncol=10))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P20 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P20",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )

P13+P14+P15+P16
P17+P18+P19+P20



##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/Pdata.txt")
dat <- dat[dat$Sex == "F",]
stat.test <- dat %>% t_test(P1 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P1 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P1",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/Pdata.txt")
dat <- dat[dat$Sex == "M",]
stat.test <- dat %>% t_test(P2 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P2 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P2",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P3 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P3 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P3",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P4 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P4 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P4",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################

P1+P2+P3+P4
P5+P6+P7+P8
P9+P10+P11+P12
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures_tmp/10-Feather/PdataF.txt")
stat.test <- dat %>% t_test(P17 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P17 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P17",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P18 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P18 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P18",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P19 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P19 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P19",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )
##################################################################################
stat.test <- dat %>% t_test(P20 ~ Feather, p.adjust.method = "bonferroni")
PLOT <- data.frame(G1 = stat.test$group1,G2 = stat.test$group2,P = as.numeric(stat.test$p.adj))
row1<- c("Black-G", "Black-G", 1)
row2<- c("Yellowish", "Yellowish", 1)
PLOT <- data.frame(rbind(PLOT,row1,row2))
PLOT <- spread(PLOT, key =G1 , value = P, fill = 1)
row.names(PLOT) <- PLOT$G2
PLOT <-PLOT[,2:7]
PLOT[upper.tri(PLOT)]  <- t(PLOT)[upper.tri(PLOT)]
PLOTF <- data.frame(matrix(as.numeric(unlist(PLOT)),nrow=6,ncol=6))
rownames(PLOTF) <- rownames(PLOT)
colnames(PLOTF) <- rownames(PLOT)
PLOTF[PLOTF <= 0.00001] <- 0.00001
P20 <- pheatmap(-log10(PLOTF), 
               color = c("white","white","#f9ae08","#f9802d","#C00000" ),
               clustering_method = "average", border_color = "gray",
               cellwidth = 10, cellheight = 10, fontsize = 10,
               cluster_cols = F, cluster_rows = F, 
               show_rownames = TRUE, show_colnames = TRUE,
               angle_col = c("90"), column_title = "P20",
               name = "-log10(p)", legend_breaks=seq(1,5,1),breaks=c(seq(1,5,by=1)) )


P13+P14+P15+P16
P17+P18+P19+P20

##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20230706-初稿/Figures/10-Feather/Fst.merge.txt", header =F)

Res <- data.frame(SNP = dat$V1,
                  Chromosome = dat$V2,
                  Position = dat$V3,
                  trait1 = dat$V4,
                  trait2 = dat$V5,
                  trait3 = dat$V6,                                   
                  trait4 = dat$V7,
                  trait5 = dat$V8             
                  )

datapca <- Res[,4:8]
data.pca <- princomp(datapca)
summary(data.pca)
ggbiplot(data.pca, obs.scale = 1, var.scale = 1, 
            varname.abbrev = TRUE,
            #ellipse = TRUE,
            labels.size = 6, 
            varname.size = 6,
            #groups = datapca$type,
            #circle= TRUE,
            alpha=0.8,
) + theme(panel.grid=element_blank()) + scale_color_d3() +
  theme(legend.direction = 'horizontal',
        legend.position = 'top',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) 


library(tsne)

tsne_out <- Rtsne(
  Res[,4:8],
  dims = 2,
  pca = FALSE,
  perplexity = 0.5,
  theta = 0.0,
  max_iter = 1000,
  check_duplicates=FALSE,
)

tsne_res <- as.data.frame(tsne_out$Y)
colnames(tsne_res) <- c("tSNE1","tSNE2")


ggplot(tsne_res,aes(tSNE1,tSNE2)) + 
  geom_point() + theme_bw() + 
  geom_hline(yintercept = 0,lty=2,col="red") + 
  geom_vline(xintercept = 0,lty=2,col="blue",lwd=1) +
  theme(plot.title = element_text(hjust = 0.5)) + 
  labs(title = "tSNE plot",color="Species")





p$layers <- c(p$layers[[2]],p$layers[[3]],p$layers[[1]])
p$layers[[3]]$aes_params$colour <- 'black'
p$layers[[2]]$aes_params$colour <- 'black'


Res$trait1 <- scale(Res$trait1) 
Res$trait2 <- scale(Res$trait2) 
Res$trait3 <- scale(Res$trait3) 
Res$trait4 <- scale(Res$trait4) 
Res$trait5 <- scale(Res$trait5) 

SNPs <-  Res[
    Res$trait1 > sort(Res$trait1,decreasing = T)[96] |
    Res$trait2 > sort(Res$trait2,decreasing = T)[96] |
    Res$trait3 > sort(Res$trait3,decreasing = T)[96] |
    Res$trait4 > sort(Res$trait4,decreasing = T)[96] |      
    Res$trait5 > sort(Res$trait5,decreasing = T)[96], 1]


library('CMplot')
CMplot(Res,type="p",plot.type="m",
       LOG10=FALSE,
       #highlight=SNPs,highlight.type="l",
       threshold=5,
       threshold.col="black",threshold.lty=1,
       col=c("grey60","#4197d8"),
       signal.cex=1.2, signal.col="red", highlight.col="grey", highlight.cex=0.7,
       #file="jpg",dpi=300,
       #file.output=FALSE,
       ylab = "", 
       #height=3, width = 6,
       #pch =".", 
       amplify=FALSE,
       verbose=TRUE,
       multracks=TRUE)





c(0.3,0.1,0.3,0.3,0.3)


CMplot(Res,plot.type="m",
       LOG10=F,
       cex=0.5,
      #ylim=c(0.1,0.4),
       file.output=F)


CMplot(Res, plot.type="m", LOG10=F,
       type="p",
       chr.den.col=NULL, 
       #col = color_set,
       #threshold = 0.3, threshold.col = "red", threshold.lwd= 2, threshold.lty =1,
       amplify = T,
       file.output=F, height=4, width = 6,
       ylab = "XP-CLR",
       pch =".", 
       cex =2
       )



