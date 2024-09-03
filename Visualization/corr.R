# Load the necessary library
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggdensity)
# Load the CSV files
A <- fread("I:/lda/0611/D/0.03mask.csv", header = FALSE)
B <- fread("I:/lda/0611/figure/group_difference/t.csv", header = FALSE)
C <- fread("I:/lda/0611/D/Dc.csv", header = FALSE)
p <- fread("I:/lda/0611/figure/group_difference/p.csv", header = FALSE)
# Ensure the data is in matrix format
#A <- as.matrix(A)
B <- as.matrix(B)
C <- as.matrix(C)
p <- as.matrix(p)
p.adjust(p,method="fdr")
# Create a sample 416x416 matrix

B[!lower.tri(B)] <- NA
C[!lower.tri(C)] <- NA
p[!lower.tri(p)] <- NA

# Identify the locations in A that contain '4'
locations <- which(A == 0, arr.ind = TRUE)

# Screen B at these locations
B[locations] <- NA
C[locations] <- NA
p[locations] <- NA
#C[C>0.0004]<-NA

#locations2 <- which(p > 0.05, arr.ind = TRUE)

#B[locations2] <- NA
#C[locations2] <- NA
#p[locations2] <- NA

# Convert the matrix to a column vector
B_vector <- as.vector(B)

# Optionally, convert the vector to a data frame for saving as CSV
B_vector_df <- data.frame(B_vector)
B_vector_df<-abs(B_vector_df)

# Convert the matrix to a column vector
C_vector <- as.vector(C)

# Optionally, convert the vector to a data frame for saving as CSV
C_vector_df <- data.frame(C_vector)

p_vector <- as.vector(p)

# Optionally, convert the vector to a data frame for saving as CSV
p_vector_df <- data.frame(p_vector)

data<-cbind(B_vector_df,C_vector_df,p_vector_df)

data_matrix <- as.matrix(data)
colnames(data_matrix)<-c("t","LDA",'p')
data_matrix<-na.omit(data_matrix)
data_matrix<-data.frame(data_matrix)
#conjuction
p<- ggplot(data_matrix,aes(x=t, y=LDA))+
  stat_cor(method = "pearson", 
           #label.x.npc = "left", 
           label.y.npc = 0.5,
           label.x.npc = 0.03)+
  #xlim(c(1.5,5))+
  ylim(c(0,0.000025))+
  labs(x = "t statistic", y = "D (ROI pairs|Patients)")+
  scale_x_continuous(breaks = seq(0,4.5, by =2),expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  #scale_y_continuous(breaks = seq(0,0.00055, by =0.0002),expand = c(0,0))+
  expand_limits(x = 0)+
  theme_classic()+
  theme(axis.line.x=element_line(linetype=1,color="black",
                                 size=1),
        axis.line.y=element_line(linetype=1,color="black",
                                 size=1))+
  theme(axis.text=element_text(size=12))+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14))+
  theme(legend.position = "none")+
  stat_density_2d(geom = "polygon",
                  aes(alpha = ..level..,fill = ..level..),
                  bins = 45)+
  scale_fill_distiller(palette="YlOrRd",direction = 1) +
  geom_smooth(aes(x=t, y=LDA), method = 'lm',
              level =0.95, se= F, linewidth= 0.6,color="#403990")


p

