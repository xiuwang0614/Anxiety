library(reshape2)
library(MASS)
library(sjPlot)
setwd('I:/lda/demo/demography')
data<- read.csv('data_GLM.csv',header = TRUE)
head(data)
#row_condition <- rowSums(!is.na(data[,14:31])) > 0
# Step 2: Subset the data frame based on the condition
#data <- data[row_condition, ]
df <- as.data.frame(data)
# Add row names to the data frame
#df$row_name <- rownames(df)
df1=df[df$diagnoses=='FoH',]
data1<-as.matrix(df1[,1:4])

long_df1<- melt(data1)
head(long_df1)
levels(long_df1$Var2) = paste0('F',1:4)

contrasts(long_df1$Var2) = contr.treatment(4,base=1)
#需要手动修改base值，不知道对不对???
ModelFoHb1 = lm(data = long_df1, value~Var2)
a1<-confint(ModelFoHb1, level = 0.95)
b1<-round(coef(summary(ModelFoHb1)),digits = 3)
b1[,4]<-round(p.adjust(b1[,4], "BH"), 3)
data1= cbind(a1,b1)
rownames(data1) <- c('Interpret','F2-F1','F3-F1','F4-F1')

contrasts(long_df1$Var2) = contr.sdif(4)
#与第i个因子的前一个相比
ModelFoHsd = lm(data = long_df1, value~Var2)
a2<-confint(ModelFoHsd, level = 0.95)
b2<-round(coef(summary(ModelFoHsd)),digits = 3)
b2[,4]<-round(p.adjust(b2[,4], "BH"), 3)
data2= cbind(a2,b2)
rownames(data2) <- c('Interpret','F2-F1','F3-F2','F4-F3')

contrasts(long_df1$Var2) = contr.treatment(4,base=4)
#需要手动修改base值，不知道对不对???
ModelFoHb4 = lm(data = long_df1, value~Var2)
plot_model(ModelFoHb4,vline.color="red")
a3<-confint(ModelFoHb4, level = 0.95)
b3<-round(coef(summary(ModelFoHb4)),digits = 3)
#FDR矫正
b3[,4]<-round(p.adjust(b3[,4], "BH"), 3)
data3= cbind(a3,b3)
rownames(data3) <- c('Interpret','F1-F4','F2-F4','F3-F4')

data = rbind(data1,data2,data3)
write.table(data,"FoH.csv",row.names=TRUE,col.names=TRUE,sep=",")


#######PD#########
df1=df[df$diagnoses=='PD',]
data1<-as.matrix(df1[,1:4])

long_df1<- melt(data1)
head(long_df1)
levels(long_df1$Var2) = paste0('F',1:4)
contrasts(long_df1$Var2) = contr.treatment(4,base=1)
#需要手动修改base值，不知道对不对???
ModelFoHb1 = lm(data = long_df1, value~Var2)
a1<-confint(ModelFoHb1, level = 0.95)
b1<-round(coef(summary(ModelFoHb1)),digits = 3)
b1[,4]<-round(p.adjust(b1[,4], "BH"), 3)
data1= cbind(a1,b1)
rownames(data1) <- c('Interpret','F2-F1','F3-F1','F4-F1')

contrasts(long_df1$Var2) = contr.sdif(4)
#与第i个因子的前一个相比
ModelFoHsd = lm(data = long_df1, value~Var2)
a2<-confint(ModelFoHsd, level = 0.95)
b2<-round(coef(summary(ModelFoHsd)),digits = 3)
b2[,4]<-round(p.adjust(b2[,4], "BH"), 3)
data2= cbind(a2,b2)
rownames(data2) <- c('Interpret','F2-F1','F3-F2','F4-F3')

contrasts(long_df1$Var2) = contr.treatment(4,base=4)
#需要手动修改base值，不知道对不对???
ModelFoHb4 = lm(data = long_df1, value~Var2)
plot_model(ModelFoHb4,vline.color="red")
a3<-confint(ModelFoHb4, level = 0.95)
b3<-round(coef(summary(ModelFoHb4)),digits = 3)
b3[,4]<-round(p.adjust(b3[,4], "BH"), 3)
data3= cbind(a3,b3)
rownames(data3) <- c('Interpret','F1-F4','F2-F4','F3-F4')

data = rbind(data1,data2,data3)
write.table(data,"PD.csv",row.names=TRUE,col.names=TRUE,sep=",")
###################GAD########################
df1=df[df$diagnoses=='GAD',]
data1<-as.matrix(df1[,1:4])

long_df1<- melt(data1)
head(long_df1)
levels(long_df1$Var2) = paste0('F',1:4)

contrasts(long_df1$Var2) = contr.treatment(4,base=1)
#需要手动修改base值，不知道对不对???
ModelFoHb1 = lm(data = long_df1, value~Var2)
a1<-confint(ModelFoHb1, level = 0.95)
b1<-round(coef(summary(ModelFoHb1)),digits = 3)
b1[,4]<-round(p.adjust(b1[,4], "BH"), 3)
data1= cbind(a1,b1)
rownames(data1) <- c('Interpret','F2-F1','F3-F1','F4-F1')

contrasts(long_df1$Var2) = contr.sdif(4)
#与第i个因子的前一个相比
ModelFoHsd = lm(data = long_df1, value~Var2)
a2<-confint(ModelFoHsd, level = 0.95)
b2<-round(coef(summary(ModelFoHsd)),digits = 3)
b2[,4]<-round(p.adjust(b2[,4], "BH"), 3)
data2= cbind(a2,b2)
rownames(data2) <- c('Interpret','F2-F1','F3-F2','F4-F3')

contrasts(long_df1$Var2) = contr.treatment(4,base=4)
#需要手动修改base值，不知道对不对???
ModelFoHb4 = lm(data = long_df1, value~Var2)
plot_model(ModelFoHb4,vline.color="red")
a3<-confint(ModelFoHb4, level = 0.95)
b3<-round(coef(summary(ModelFoHb4)),digits = 3)
b3[,4]<-round(p.adjust(b3[,4], "BH"), 3)
data3= cbind(a3,b3)
rownames(data3) <- c('Interpret','F1-F4','F2-F4','F3-F4')

data = rbind(data1,data2,data3)
write.table(data,"GAD.csv",row.names=TRUE,col.names=TRUE,sep=",")
###################SAD########################
df1=df[df$diagnoses=='SAD',]
data1<-as.matrix(df1[,1:4])

long_df1<- melt(data1)
head(long_df1)
levels(long_df1$Var2) = paste0('F',1:4)

contrasts(long_df1$Var2) = contr.treatment(4,base=1)
#需要手动修改base值，不知道对不对???
ModelFoHb1 = lm(data = long_df1, value~Var2)
a1<-confint(ModelFoHb1, level = 0.95)
b1<-round(coef(summary(ModelFoHb1)),digits = 3)
b1[,4]<-round(p.adjust(b1[,4], "BH"), 3)
#b1[,4]<-round(b1[,4], 3)
data1= cbind(a1,b1)
rownames(data1) <- c('Interpret','F2-F1','F3-F1','F4-F1')

contrasts(long_df1$Var2) = contr.sdif(4)
#与第i个因子的前一个相比
ModelFoHsd = lm(data = long_df1, value~Var2)
a2<-confint(ModelFoHsd, level = 0.95)
b2<-round(coef(summary(ModelFoHsd)),digits = 3)
b2[,4]<-round(p.adjust(b2[,4], "BH"), 3)
#b2[,4]<-round(b2[,4], 3)
data2= cbind(a2,b2)
rownames(data2) <- c('Interpret','F2-F1','F3-F2','F4-F3')

contrasts(long_df1$Var2) = contr.treatment(4,base=4)
#需要手动修改base值，不知道对不对???
ModelFoHb4 = lm(data = long_df1, value~Var2)
plot_model(ModelFoHb4,vline.color="red")
a3<-confint(ModelFoHb4, level = 0.95)
b3<-round(coef(summary(ModelFoHb4)),digits = 3)
b3[,4]<-round(p.adjust(b3[,4], "BH"), 3)
#b3[,4]<-round(b3[,4], 3)
data3= cbind(a3,b3)
rownames(data3) <- c('Interpret','F1-F4','F2-F4','F3-F4')

data = rbind(data1,data2,data3)
write.table(data,"SAD.csv",row.names=TRUE,col.names=TRUE,sep=",")
