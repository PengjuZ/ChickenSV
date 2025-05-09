library("pheatmap")
library('RColorBrewer')
library('gplots')
suppressMessages(library('DESeq2'))
library(ComplexHeatmap)
library(circlize)
library("WGCNA")
library("stringi")
library('ggplot2')
library(Rtsne)
theme_set(theme_bw())

data<-read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/k/gene-Os.tsv",row.names = 1)
data<-data[,c(12:23,27:29)]

keep <- rowSums(data) >= 50
data <- data[keep,]
group_list <- factor(c(rep("Colored",3),rep("White",3),rep("Colored",3),rep("Colored",3),rep("White",3)))
#group_list <- factor(c(rep("SW",3),rep("WL",3),rep("SL",3),rep("TB",3),rep("WC",3)))
#group_list <- factor(c(rep("WL",6),rep("SD",3),rep("WL",3),rep("SL",3),rep("TB",3),rep("GP",3),rep("WCB",3),rep("WCY",3)))
colData <- data.frame(row.names=colnames(data),group_list=group_list)

##标准化数据
dds <- DESeqDataSetFromMatrix(countData = data,colData = colData, design = ~ group_list)
dds2 <- estimateSizeFactors(dds)
counts <- counts(dds2, normalized=F)
normalized_counts <- counts(dds2, normalized=TRUE)
rlds <- data.frame(normalized_counts)
write.csv(rlds,file = "G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/k/normalized_counts.csv")

##标准化后展示
rld <- rlogTransformation(dds2) 
vsd<-varianceStabilizingTransformation(dds2, blind=TRUE)
data_new=assay(rld)
n.sample=ncol(data)
cols <- rainbow(n.sample*1.2)
boxplot(data, col = cols,main="expression value",las=2)
boxplot(data_new, col = cols,main="expression value",las=2)

##组间聚类
hmcol<- colorRampPalette(brewer.pal(9, 'GnBu'))(100)
distsRL <- dist(t(rlds))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),paste(group_list,colnames(data) , sep=' : '))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),symm=TRUE, trace='none',col = rev(hmcol), margin=c(13, 13))


##差异分析
dds <- DESeqDataSetFromMatrix(countData = data,colData = colData, design = ~ group_list)
suppressMessages(dds2 <- DESeq(dds))
res <- results(dds2,contrast = c("group_list","White","Colored"))
resAB <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
summary(resAB)
resOrdered <- resAB[order(resAB$padj),]
resOrdered=as.data.frame(resOrdered)
head(resOrdered)
resdata <- merge(resOrdered, as.data.frame(counts(dds2, normalized=TRUE)),by="row.names",sort=FALSE)
write.csv(resdata,file = "G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/k/res.csv")


theme_set(theme_bw())
resr0 <-subset(res,padj > 0.05)
resr1 <-subset(res,padj < 0.05 & (log2FoldChange < 1 & log2FoldChange > -1))
resr2 <-subset(res,padj < 0.05 & (log2FoldChange > 1))
resr3 <-subset(res,padj < 0.05 & (log2FoldChange < -1))

resr0$significant<-"no"
resr1$significant<-"no"
resr2$significant<-"up"
resr3$significant<-"down"

resr<-as.data.frame(rbind(resr0,resr1,resr2,resr3))
r = ggplot(resr,aes(log2FoldChange,-1*log10(padj)))+ geom_point(aes(color=significant)) + labs(title="",x="log2(FC)",y="-log10(padj)")
r+geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
scale_color_manual(values = c("#00ba38","#619cff","#f8766d"))
#r03 + geom_point()+ xlim(-4,4) + ylim(0,30)

##单个基因
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#DE3C3C","#09194F","#1f640a","#490a3d")###################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/k/001-ID.list.txt",row.names = 1)
data$Strains <- factor(data$breed, levels = c("Silkie","Tibetan","Sweden","White Wenchang","White Leghorn"))
P1<-ggboxplot(data, x = "Strains", y = "LOC107052320", width = 0.6, 
                color = "black",#轮廓颜色
                fill="Strains",#填充
                #alpha = 0.8,
                palette = colors10,#分组着色
                bxp.errorbar=T,#显示误差条
                bxp.errorbar.width=0.5, #误差条大小
                #size=1, #箱型图边线的粗细
                x.text.angle=30,
                outlier.shape=NA, #不显示outlier
                #legend = "right"
                #rotate=TRUE,
                xlab = "", 
                ylab = "Normalize counts"
     ) + theme (legend.position="none")

data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/k/002-ID.list.txt")
P2<-ggboxplot(data, x = "color", y = "value", 
          #width = 0.6, 
          color = "black",#轮廓颜色
          fill="color",#填充
          #alpha = 0.8,
          palette =c("#C00000","#005eaa"),#分组着色
          xlab = F, #不显示x轴的标签
          bxp.errorbar=T,#显示误差条
          bxp.errorbar.width=0.5, #误差条大小
          #size=1, #箱型图边线的粗细
          outlier.shape=NA, #不显示outlier
          ylab = '2^-(∆∆Ct)',
          x.text.angle=30,
          #legend = "right"
)  + theme (legend.position="none") + 
  stat_compare_means(method = "t.test")


plot_grid(P1,P2, ncol=2, nrow=1,align = 'h', rel_widths = c(3,2))




##############




#############################################
data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/d/exon-O.ok.tsv",row.names = 1)

data <- data[,6:23]

#设置分组信息以及构建dds对象
group_list <- factor(c(rep("White",6),rep("Colored",3),rep("White",3),rep("Colored",6)))
colData <- data.frame(row.names=colnames(data),group_list=group_list)

##标准化数据
dds <- DESeqDataSetFromMatrix(countData = data,colData = colData, design = ~ group_list)
dds2 <- estimateSizeFactors(dds)
counts <- counts(dds2, normalized=F)
normalized_counts <- counts(dds2, normalized=TRUE)
rlds <- data.frame(normalized_counts)
write.csv(rlds,file = "G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/d/normalized_counts.exon.csv")

##单个基因
datS <- read.csv("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure8/d/normalized_counts.exon.sl.csv",row.names = 1,header = F)
data_scale <- as.data.frame(apply(datS,1,scale))
rownames(data_scale) <- c("Leghorn1","Leghorn2","Leghorn3","Leghorn4","Leghorn5","Leghorn6","Sweden1","Sweden2","Sweden3","Leghorn7","Leghorn8","Leghorn9","Silkie1","Silkie2","Silkie3","Tibetan1","Tibetan2","Tibetan3")

Heatmap(data_scale,name = "heatmap",
        column_names_side = "bottom",
        col = colorRamp2(c(-1,0,1), c("#005eaa", "#FFC000", "#C00000")),
        cluster_columns = F,
        cluster_rows = T,
        row_dend_side = "left",
        show_row_names = TRUE,
        row_names_max_width = unit(6, "cm"),
        row_names_gp = gpar(fontsize = 9),
        clustering_method_rows = "complete")








