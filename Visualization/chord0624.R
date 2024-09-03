library(tidyverse)
setwd('I:/lda/0611/figure/E/E2')
EC<- read.csv('F4.csv',sep = ",",header = F,row.names = 1)
ECmatrix<-as.matrix(EC)
# Define the ranges
ranges <- list(c(1,52), c(53,99), c(100,125), c(126,216),
               c(217,262),c(263,323),c(324,400),
               c(401,416))

# Initialize an empty 4x4 matrix to store the results
result <- matrix(nrow=8, ncol=8)

# Calculate the mean for each specific point in the matrix
for(i in 1:length(ranges)){
  for(j in 1:length(ranges)){
    result[i,j] <- mean(ECmatrix[ranges[[i]][1]:ranges[[i]][2], ranges[[j]][1]:ranges[[j]][2]])
  }
}

colnames(result)<- c("Cont","SAN","LN","DMN","DAN","Vis","SMN","Sub")
rownames(result)<- c("Cont","SAN","LN","DMN","DAN","Vis","SMN","Sub")

# Specify the full path to the folder where you want to save the result data
output_folder <- "I:/lda/0611/figure/E/E2"
# Update the output file paths
ec_output_file <- file.path(output_folder, "chordF4.csv")
# Write the matrices to CSV files
write.csv(result, file = ec_output_file, row.names = TRUE)


# figure the connectome figure
library(reshape2)
library(circlize)
setwd('I:/lda/0611/figure/E/E2')
ECmatrix <- read.csv('chordF4.csv',sep = ",",header = T, row.names = 1)
# Convert the matrix to a data frame
df <- as.data.frame(ECmatrix)
# Add row names to the data frame
#df$row_name <- rownames(df)
df<-as.matrix(df)

# 提取矩阵对角线数据
diag(df)
df[!upper.tri(df, diag = TRUE)] <- NA

dim(df)
#df<-df[,1:9]
# Use melt to convert the data frame to a long format
long_df<- melt(df)
head(long_df)
long_df<- long_df %>% filter(value != "NA")
colnames(long_df)<-c("from","to","value")
dat<- long_df

# figure the conenctome fiugur

folder_path <- "I:/lda/0611/figure/E/E2"
tiff_file <- file.path(folder_path, "chordF4.tiff")

tiff(filename = tiff_file, width = 960, height = 960,units = "px", pointsize = 12,
     #compression = c("none", "rle", "lzw", "jpeg", "zip", "lzw+p", "zip+p"),
     res = 130)


circos.clear()
circos.par(start.degree = 90, 
           gap.degree = 8, 
           points.overflow.warning = FALSE,
           gap.after=c("Cont"=3,"SAN"=3,"LN"=3,
                       "DMN"=3,"DAN"=3,"Vis"=3,"SMN"=3, "Sub"=3))

#线段的颜色范围取决于数据的最小值与最大值
cols <- colorRamp2(c(-8.611055e-06,0,1.155609e-05), c("#2632CD","white","#FF3F3F"))


grid.col<-c("Cont"="#ebce6a","SAN"="#f5a889","LN"="#ea71ae",
            "DMN"="#b22222","DAN"="#383838","Vis"="#000080","SMN"="#5ca7c7",
            "Sub"="#acd6ec")
#group<-c("1"="A","2"="A","3"="A","4"="A",
#         "5"="B","6"="B","7"="B","8"="B",
#         "9"="C","10"="C","11"="C","12"="C",
#         "13"="D","14"="D","15"="D","16"="D","17"="D","18"="D")
#col_fun = function(x) ifelse(x < 0.001, "#00000000", "#FF000080")
chordDiagram(dat,
             grid.col = grid.col,
             #col = rand_color(nrow(dat)),
             #col=col_fun,
             #group = group,
             transparency = 0.2,
             #directional = 1,
             #direction.type = c("diffHeight", "arrows"), 
             #diffHeight  = -0.04,
             annotationTrack = "grid", 
             annotationTrackHeight = c(0.08, 0.1),
             link.arr.type = "big.arrow", 
             link.sort = TRUE, 
             link.decreasing = TRUE,
             link.largest.ontop = TRUE,
             preAllocateTracks = list(
               track.height = 0.1,
               track.margin = c(0.01, 0.02)),
             col = cols)

circos.trackPlotRegion(
  track.index = 2, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    #sector.index = get.cell.meta.data("sector.index")
    sector.index = gsub("[a-z]+_", "", get.cell.meta.data("sector.index"))
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 2, 
      col = "black",
      labels = sector.index, 
      facing = "bending", 
      cex = 1.4,
      niceFacing = TRUE
    )
  }
)

dev.off()