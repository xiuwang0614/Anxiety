setwd('I:/lda/0611/figure/E/E2')
library(pheatmap)
data1<- read.csv("I:/lda/0611/thresholded/0.025/reF1/F1.csv",sep = ",",header = F,row.names = 1)
head(data1)
data1<-as.matrix(data1)
colnames(data1)<-c(1:416)

data2<- read.csv("F2.csv",sep = ",",header = F,row.names = 1)
head(data2)
data3<-as.matrix(data2)
colnames(data2)<-c(1:416)

data3<- read.csv("I:/lda/0611/thresholded/0.025/reF1/F3.csv",sep = ",",header = F,row.names = 1)
head(data3)
data3<-as.matrix(data3)
colnames(data3)<-c(1:416)

data4<- read.csv("F4.csv",sep = ",",header = F,row.names = 1)
head(data4)
data4<-as.matrix(data4)
colnames(data4)<-c(1:416)

datarow<- read.csv('networks0624.csv',sep = ",",header = F,row.names =1)
colnames(datarow)<-"networks"
datarow2<- read.csv('I:/lda/0611/figure/E/E1/reF1/network0625.csv',sep = ",",header = F,row.names =1)
colnames(datarow2)<-"networks"
datacol=datarow
datacol2=datarow2
colors<-list("networks" = c(Cont= "#D3B9B9",
                            SalVentAttn = "#784C4D",
                            Limbic = "#D92523",
                            #TempPar = "#9B3782",
                            Default = "#62498A",
                            DorsAttn = '#0098A3',
                            Vis='#00945E',
                            SomMot='#797B00',
                            Subcortical='#FFE5DB'))

F1<- pheatmap(data1,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 3,
         cellheight = 3,
         annotation_row = datarow2,
         annotation_col = datacol2,
         annotation_colors = colors,
         gaps_row = c(77,129,176,222,283,374,400),
         gaps_col = c(77,129,176,222,283,374,400),
         fontsize_row=20, fontsize_col =10, fontsize=9,
         #F3
         #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#fe2701","#FE2701", "#FA300D" ,"#F73C1B","#F04E32" ,"#F05135" ,"#f9c3bf"))(101),
         #F1
         color = colorRampPalette(c("#bfe4ee","#6E94EA","#163EE7","#2632CD","black","black","black","black","#E0281D","#fe2701","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
         #F2
         #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#FA300D","#F04E32","#f9c3bf"))(101),
         #f4
         #color = colorRampPalette(c("#bfe4ee","#2632CD","black","black","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
         border = F,
         show_rownames = F, show_colnames = F,
         column_names_side = c("top"),
         angle_col = c("90"))
         #main = "Factor 1") # a title for our heatmap)
        # filename = "0624k4.tiff")

F2<- pheatmap(data2,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              cellwidth = 3,
              cellheight = 3,
              annotation_row = datarow,
              annotation_col = datacol,
              annotation_colors = colors,
              gaps_row = c(52,99,125,216,262,323,400),
              gaps_col = c(52,99,125,216,262,323,400),
              fontsize_row=20, fontsize_col =10, fontsize=9,
              #F3
              #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#fe2701","#FE2701", "#FA300D" ,"#F73C1B","#F04E32" ,"#F05135" ,"#f9c3bf"))(101),
              #F1
              #color = colorRampPalette(c("#bfe4ee","#6E94EA","#163EE7","#2632CD","black","black","black","black","#E0281D","#fe2701","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              #F2
              color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#FA300D","#F04E32","#f9c3bf"))(101),
              #f4
              #color = colorRampPalette(c("#bfe4ee","#2632CD","black","black","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              border = F,
              show_rownames = F, show_colnames = F,
              column_names_side = c("top"),
              angle_col = c("90"))
              #main = "Factor 2") # a title for our heatmap)
              #filename = "0624k4.tiff")

F3<- pheatmap(data3,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              cellwidth = 3,
              cellheight = 3,
              annotation_row = datarow2,
              annotation_col = datacol2,
              annotation_colors = colors,
              gaps_row = c(77,129,176,222,283,374,400),
              gaps_col = c(77,129,176,222,283,374,400),
              fontsize_row=40, fontsize_col =40, fontsize=9,
              #F3
              color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#fe2701","#FE2701", "#FA300D" ,"#F73C1B","#F04E32" ,"#F05135" ,"#f9c3bf"))(101),
              #F1
              #color = colorRampPalette(c("#bfe4ee","#6E94EA","#163EE7","#2632CD","black","black","black","black","#E0281D","#fe2701","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              #F2
              #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#FA300D","#F04E32","#f9c3bf"))(101),
              #f4
             # color = colorRampPalette(c("#bfe4ee","#2632CD","black","black","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              border = F,
              show_rownames = F, show_colnames = F,
              column_names_side = c("top"),
              angle_col = c("90"))
             #main = "Factor 3")

F4<- pheatmap(data4,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              cellwidth = 3,
              cellheight = 3,
              annotation_row = datarow,
              annotation_col = datacol,
              annotation_colors = colors,
              gaps_row = c(52,99,125,216,262,323,400),
              gaps_col = c(52,99,125,216,262,323,400),
              fontsize_row=20, fontsize_col =10, fontsize=9,
              #F3
              #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#fe2701","#FE2701", "#FA300D" ,"#F73C1B","#F04E32" ,"#F05135" ,"#f9c3bf"))(101),
              #F1
              #color = colorRampPalette(c("#bfe4ee","#6E94EA","#163EE7","#2632CD","black","black","black","black","#E0281D","#fe2701","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              #F2
              #color = colorRampPalette(c("#bfe4ee","#0b34e7","black","black","#fe2701","#FA300D","#F04E32","#f9c3bf"))(101),
              #f4
              color = colorRampPalette(c("#bfe4ee","#2632CD","black","black","#FA300D","#F73A19","#F73A19" ,"#F04E32","#f9c3bf"))(101),
              border = F,
              show_rownames = F, show_colnames = F,
              column_names_side = c("top"),
              angle_col = c("90"))
              #main = "Factor 1") # a title for our heatmap)
              #filename = "0624k4.tiff")

library(gridExtra)
library(grid)
combined <- grid.arrange(grobs = list(F1[[4]], F2[[4]], F3[[4]], F4[[4]]), ncol = 2, nrow = 2)
# Save combined image to file
tiff(filename = "combined_heatmaps2.tiff", width = 3000, height = 3000)
grid.draw(combined)

# Add titles with customized sizes
grid.text("Factor 1", x = 0.25, y = 0.99, gp = gpar(fontsize = 50))
grid.text("Factor 2", x = 0.75, y = 0.99, gp = gpar(fontsize = 50))
grid.text("Factor 3", x = 0.25, y = 0.49, gp = gpar(fontsize = 50)) # Customized title size for Factor 3
grid.text("Factor 4", x = 0.75, y = 0.49, gp = gpar(fontsize = 50))

dev.off()