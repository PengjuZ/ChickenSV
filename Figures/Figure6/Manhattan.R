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
library(cowplot)
library(forcats)
library(dplyr)
library(ggsci)
library(fitdistrplus)
library('CMplot')
##################################################################################
colors10 <-c("#C00000","#005eaa","#ED7D31","#FFC000","#1E925E","#791E94","#09194F","#DE3C3C","#1f640a","#490a3d")
##################################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure7/01-Manhattan/Merge.5.list", header =F)
dat$ID <- paste(dat$V1, dat$V2, sep = "_")
colnames(dat) <- c('Chr','Start','End','FST','XPCLR','XPEHH','ID') 
dat[dat$FST > quantile(dat$FST,c(0.999)) & dat$XPCLR > quantile(dat$XPCLR,c(0.999)) & dat$XPEHH > quantile(dat$XPEHH,c(0.999)),]

#dat<-dat[dat$Chr == "chr1" | dat$Chr == "chr3" | dat$Chr == "chr4" | dat$Chr == "chr7" | dat$Chr == "chrZ",]
dat<-dat[dat$Chr == "chr1" | dat$Chr == "chrZ",]
dat$FSTr <- rank(dat$FST)/max(rank(dat$FST))
dat$XPCLRr <- rank(dat$XPCLR)/max(rank(dat$XPCLR))
dat$XPEHHr <- rank(dat$XPEHH)/max(rank(dat$XPEHH))

dat$FSTs <- dat$FSTr/(1.001-dat$FSTr)
dat$XPCLRs <- dat$XPCLRr/(1.001-dat$XPCLRr)
dat$XPEHHs <- dat$XPEHHr/(1.001-dat$XPEHHr)
dat$Rank <- (dat$FSTs+dat$XPCLRs+dat$XPEHHs)/3

Res <- data.frame(SNP = dat$ID,
                  Chromosome = dat$Chr,
                  Position = dat$Start,
                  trait1 = dat$Rank)

Res[Res$trait1 > quantile(Res$trait1,c(0.999)),]

CMplot(Res, type="p",
       plot.type="m",
       LOG10=FALSE,
       threshold=c(quantile(Res$trait1,c(0.999)) ,quantile(Res$trait1,c(0.999)) ),
       signal.cex=c(0.8,0.8),
       threshold.lty=2, 
       threshold.lwd=1,  
       col=c("#C00000","#005eaa"),
       #signal.col="red", 
       threshold.col=c("grey"),
       file="jpg",file.name="",dpi=300,
       file.output=TRUE,verbose=TRUE,
       width=10,height=3,
       cex=0.7,
       ylab="Rank value",
       #chr.border=TRUE,
       #chr.labels.angle=45
       )

###################################################################
dat <- read.delim("G:/桌面/投稿文章/Wenchang/20231107-二次稿/Figures/Figure7/01-Manhattan/logistic.results.plot.txt", header =F)
dat$ID <- paste(dat$V1, dat$V2, sep = "_")
colnames(dat) <- c('Chr','Start','Pvalue','ID') 
#dat<-dat[dat$Chr == "1" | dat$Chr == "3" | dat$Chr == "4" | dat$Chr == "7" | dat$Chr == "chrZ",]
dat<-dat[dat$Chr == "1" | dat$Chr == "chrZ",]
Res <- data.frame(SNP = dat$ID,
                  Chromosome = dat$Chr,
                  Position = dat$Start,
                  trait1 = dat$Pvalue
                  
)

CMplot(Res, type="p",
       plot.type="m",
       LOG10=TRUE,
       threshold=c(5.39942E-08,5.39942E-08),
       signal.cex=c(0.8,0.8),
       threshold.lty=2, 
       threshold.lwd=1,  
       col=c("#C00000","#005eaa"),
       #signal.col="red", 
       threshold.col=c("grey"),
       file="jpg",file.name="",dpi=300,
       file.output=T,verbose=TRUE,
       width=10,height=3,
       cex=0.8,
       #ylab="",
       #chr.border=TRUE,
       #chr.labels.angle=45
)




















FST<-fitdist(dat$XPEHH, "norm")

fitdist(dat$XPCLR, "unif")


"norm", "lnorm", "pois", "exp", "gamma", "nbinom", "geom", "beta", "unif",




dat$FSTz <- scale(dat$FST)
dat$XPCLRz <- scale(dat$XPCLR)
dat$XPEHHz <- scale(dat$XPEHH)


ggplot(dat, aes(x = Rank))  +   geom_density()








#dat <- dat[,c(4:6)]
#dat6 <- dat[dat$FST > quantile(dat$FST,c(0.99)) & dat$XPCLR > quantile(dat$XPCLR,c(0.99)) & dat$XPEHH > quantile(dat$XPEHH,c(0.99)),]





CMplot(Res, plot.type="m", 
       multracks=TRUE, 
       LOG10=FALSE,
       threshold=c(1e-6,1e-4),
       threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), 
       threshold.col=c("black","grey"),
       amplify=TRUE, bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"),
       signal.col=c("red","green"),
       signal.cex=c(0.5,1),
       file="jpg",dpi=300,file.output=FALSE,verbose=TRUE)










Res <- data.frame(SNP = res$ID,
                  Chromosome = res$Chr,
                  Position = res$Start,
                  trait1 = res$FST,
                  trait2 = res$XPCLR,
                  trait3 = res$XPEHH                                   
)

CMplot(Res,plot.type="c",
       LOG10=FALSE,
       r=0.4,cir.legend=TRUE,outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",dpi=300,file.output=TRUE,verbose=TRUE)


CMplot(Res,plot.type="c",
       LOG10=FALSE,
       outward=TRUE,col=matrix(c("#4DAF4A",NA,NA,"dodgerblue4","deepskyblue",NA,"dodgerblue1", "olivedrab3", "darkgoldenrod1"), nrow=3, byrow=TRUE),
       threshold=NULL,r=1.2,cir.chr.h=1.5,cir.legend.cex=0.5,cir.band=1,file="jpg",dpi=300,chr.den.col="black",file.output=TRUE,verbose=TRUE)




CMplot(Res, plot.type="m", 
       multracks=TRUE, 
       LOG10=FALSE,
       threshold=c(1e-6,1e-4),threshold.lty=c(1,2), 
       threshold.lwd=c(1,1), threshold.col=c("black","grey"), amplify=TRUE,bin.size=1e6,
       chr.den.col=c("darkgreen", "yellow", "red"), signal.col=c("red","green"),signal.cex=c(1,1),
       file="jpg",dpi=300,file.output=TRUE,verbose=TRUE)





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
       amplify=TRUE,
       verbose=TRUE,
       multracks=TRUE)



data.pca <- princomp(dat)

ggbiplot(data.pca, 
         #obs.scale = 1, var.scale = 1, 
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











dat$FSTp <- pnorm (dat$FST, mean=0.02125450, sd=0.01916994, log = FALSE)



fit <- fitdist(dat$XPCLR, "pois", method = "mme")
plot(fit)


ppois (dat$XPCLR, 4.812729e-08) 

plot(density(dat$XPCLR))

ggplot(dat, aes(x = FST))  +   geom_density()


ecdf(test)




quantile(dat$V6,c(0.01))

ggscatter(dat, x = "V4", y = "V6",
          #color = "Frequency",
          #fill= "Frequency",
          palette = c("#C00000","#ED7D31","#FEC306","#FAE7B3"),
          size=1,
          alpha = 0.8) + 
  #color = "Class", palette = my36colors, fill= "Class",ellipse = TRUE,size=2,ellipse.level = 0,ellipse.alpha=0) + 
  xlab(paste0("UMAP_1")) + ylab(paste0("UMAP_2")) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + theme(legend.position = 'none') +
  ggtitle ("CK") + theme(plot.title = element_text(hjust = 0.5))














datapca <- dat[,4:6]
data.pca <- princomp(datapca)
summary(data.pca)
ggbiplot(data.pca, 
         #obs.scale = 1, var.scale = 1, 
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


















Res <- data.frame(SNP = dat$V1,
                  Chromosome = dat$V2,
                  Position = dat$V3,
                  trait1 = dat$V4,
                  trait2 = dat$V5,
                  trait3 = dat$V6,                                   
                  trait4 = dat$V7,
                  trait5 = dat$V8             
                  )





library(uwot)
Pumap <- uwot::umap(Res[,4:8],n_neighbors = 5, min_dist = 0)
Pumap <- data.frame(Pumap)
All_res <- data.frame(Pumap,Res[,4:8])
colnames(Pumap) <- c('X','Y') 
ggscatter(Pumap, x = "X", y = "Y",
                 #color = "Frequency",
                 #fill= "Frequency",
                 palette = c("#C00000","#ED7D31","#FEC306","#FAE7B3"),
                 size=1,
                 alpha = 0.1) + 
  #color = "Class", palette = my36colors, fill= "Class",ellipse = TRUE,size=2,ellipse.level = 0,ellipse.alpha=0) + 
  xlab(paste0("UMAP_1")) + ylab(paste0("UMAP_2")) +
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank()) + theme(legend.position = 'none') +
  ggtitle ("CK") + theme(plot.title = element_text(hjust = 0.5))



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



