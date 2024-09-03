#install.packages("devtools")
#devtools::install_github("NightingaleHealth/ggforestplot")
library(ggplot2)
library(tidyverse)
library(forestploter)
install.packages("grid")
library(grid)
setwd('I:/lda/demo/demography')
data<- read.csv('FoH.csv',header = TRUE)
colnames(data)<-c('Δ','CI1','CI2','Estimate','SE','t','p')
head(data)
data$ID<-factor(data$Δ, levels = data$Δ) #将类别转换为因子类型，由于R语音默认为按首字母进行排序，因此我们手动设置顺序，即可自定义Y轴顺序
data$` `<-paste(rep(" ",20),collapse = " ")
tm<-forest_theme(base_size = 10,
                 reline_col = "red4")

p<-forest(data[,c(1,9,7)],
          est = data$Estimate,
          lower = data$CI1,
          upper = data$CI2,
          #size = data[,5],
          ci =2,
          ref_line = 0,
          theme = tm)
pp<- edit_plot(p,col = c(1,3),part = "header", gp = gpar(fontface = "italic"))
pp<-edit_plot(pp,row = c(3,5,6),gp = gpar(col = "red4",fontface = "bold"))
pp<-insert_text(pp,
                text  = "(A) Acrophobia",
                col =1,
                part = "header",
                gp = gpar(fontface="bold"))
pp
################SAD#####################
data<- read.csv('SAD.csv',header = TRUE)
colnames(data)<-c('Δ','CI1','CI2','Estimate','SE','t','p')
head(data)
data$ID<-factor(data$Δ, levels = data$Δ) #将类别转换为因子类型，由于R语音默认为按首字母进行排序，因此我们手动设置顺序，即可自定义Y轴顺序
data$` `<-paste(rep(" ",20),collapse = " ")
tm<-forest_theme(base_size = 10,
                 reline_col = "red4")

p1<-forest(data[,c(1,9,7)],
          est = data$Estimate,
          lower = data$CI1,
          upper = data$CI2,
          #size = data[,5],
          ci =2,
          ref_line = 0,
          theme = tm)
pp1<- edit_plot(p1,col = c(1,3),part = "header", gp = gpar(fontface = "italic"))
pp1<-edit_plot(pp1,row = c(5,6),gp = gpar(col = "red4",fontface = "bold"))
pp1<-insert_text(pp1,
                text  = "(B) SAD",
                col = 1,
                part = "header",
                gp = gpar(fontface="bold"))
pp1
########################GAD#################
data<- read.csv('GAD.csv',header = TRUE)
colnames(data)<-c('Δ','CI1','CI2','Estimate','SE','t','p')
head(data)
data$ID<-factor(data$Δ, levels = data$Δ) #将类别转换为因子类型，由于R语音默认为按首字母进行排序，因此我们手动设置顺序，即可自定义Y轴顺序
data$` `<-paste(rep(" ",20),collapse = " ")
tm<-forest_theme(base_size = 10,
                 reline_col = "red4")

p2<-forest(data[,c(1,9,7)],
           est = data$Estimate,
           lower = data$CI1,
           upper = data$CI2,
           #size = data[,5],
           ci =2,
           ref_line = 0,
           theme = tm)
pp2<- edit_plot(p2,col = c(1,3),part = "header", gp = gpar(fontface = "italic"))
pp2<-edit_plot(pp2,row = c(1:4),gp = gpar(col = "red4",fontface = "bold"))
pp2<-insert_text(pp2,
                 text  = "(C) GAD",
                 col = 1,
                 part = "header",
                 gp = gpar(fontface="bold"))
pp2
######################PD#########################

data<- read.csv('PD.csv',header = TRUE)
colnames(data)<-c('Δ','CI1','CI2','Estimate','SE','t','p')
head(data)
data$ID<-factor(data$Δ, levels = data$Δ) #将类别转换为因子类型，由于R语音默认为按首字母进行排序，因此我们手动设置顺序，即可自定义Y轴顺序
data$` `<-paste(rep(" ",20),collapse = " ")
tm<-forest_theme(base_size = 10,
                 reline_col = "red4")

p3<-forest(data[,c(1,9,7)],
           est = data$Estimate,
           lower = data$CI1,
           upper = data$CI2,
           #size = data[,5],
           ci =2,
           ref_line = 0,
           theme = tm)
pp3<- edit_plot(p3,col = c(1,3),part = "header", gp = gpar(fontface = "italic"))
pp3<-edit_plot(pp3,row = c(1:4,6),gp = gpar(col = "red4",fontface = "bold"))
pp3<-insert_text(pp3,
                 text  = "(D) PD",
                 col = 1,
                 part = "header",
                 gp = gpar(fontface="bold"))
pp3

library(gridExtra)
merged_plot<- grid.arrange(pp,pp1,pp2,pp3,
                           ncol=2,nrow=2)
merged_plot
