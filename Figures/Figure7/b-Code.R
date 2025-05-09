my36colors <-c("#DE3C3C","#005eaa","#e41749","#FF5A09","#F7B32D","#FFEB28","#CABD37","#41D3BD","#2BDE73","#00818A","#1f640a","#41b6e6","#368cbf","#001871","#09194F","#6600CC","#6a60a9","#791E94","#490a3d","#000000","#e41749","#FF5A09","#F7B32D","#FFEB28","#CABD37","#41D3BD","#2BDE73","#00818A","#1f640a","#41b6e6","#368cbf","#001871","#09194F","#6600CC","#6a60a9","#791E94","#490a3d","#000000")

data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/b/GRM5.plot.txt",header=F,row.names = 1)
#25,45,65,85,105,125,145,165,185,205,225,249,274,294,314
P1 <-pheatmap(data,
         color = colorRamp2(c(-1,0,1,2), c("#09194F", "#FCF9DE", "#C00000", "#C00000")),
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F, 
         show_rownames = F,
         #gaps_row = c(165,225,249,274,294,314),
         #gaps_col = c(4),
         gaps_row = c(25,165,225,249,274,294,314),
         legend = F,
         #border_color = "gray",        
)

data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/b/TYR.plot.txt",header=F,row.names = 1)
#25,45,65,85,105,125,145,165,185,205,225,249,274,294,314
P2 <-pheatmap(data,
         color = colorRamp2(c(-1,0,1,2), c("#09194F", "#FCF9DE", "#C00000", "#C00000")),
         cluster_rows = F,
         cluster_cols = F,
         show_colnames = F, 
         show_rownames = F,
         #gaps_row = c(165,225,249,274,294,314),
         #gaps_col = c(4),
         gaps_row = c(25,165,225,249,274,294,314),
         legend = F,
         #border_color = "gray",        
)

data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/b/LOC107052320.plot.txt",header=F,row.names = 1)
#25,45,65,85,105,125,145,165,185,205,225,249,274,294,314
P3 <-pheatmap(data,
              color = colorRamp2(c(-1,0,1,2), c("#09194F", "#FCF9DE", "#C00000", "#C00000")),
              cluster_rows = F,
              cluster_cols = F,
              show_colnames = F, 
              show_rownames = F,
              #gaps_row = c(165,225,249,274,294,314),
              #gaps_col = c(4),
              gaps_row = c(25,165,225,249,274,294,314),
              legend = F,
              #border_color = "gray",        
)

data <- read.delim("G:/桌面/投稿文章/Wenchang/20240311/Figures/Figure7/b/NOX4.plot.txt",header=F,row.names = 1)
#25,45,65,85,105,125,145,165,185,205,225,249,274,294,314
P4 <-pheatmap(data,
              color = colorRamp2(c(-1,0,1,2), c("#09194F", "#FCF9DE", "#C00000", "#C00000")),
              cluster_rows = F,
              cluster_cols = F,
              show_colnames = F, 
              show_rownames = F,
              #gaps_row = c(165,225,249,274,294,314),
              #gaps_col = c(4),
              gaps_row = c(25,165,225,249,274,294,314),
              legend = F,
              #border_color = "gray",        
)




 plot_grid(as.ggplot(P1), as.ggplot(P2), as.ggplot(P3), as.ggplot(P4),
           ncol=4, nrow=1, 
           labels = c('', '', '', '')
                #rel_widths = c(0.03,1) ,
               # rel_heights = c(1)
           )








