setwd('I:/lda/0803')
data<- read.csv('zhongjie_data.csv',header = TRUE)
head(data)

install.packages('mediation')
library(mediation)
b<- lm(ZDE0.025 ~ CTQSUM + DIG+age+sex,data = data)
c<- lm(HAMASUM ~ CTQSUM + ZDE0.025 + DIG + age + sex, data = data)

contcont<- mediate(b,c,sims= 5000,treat = "CTQSUM",mediator = "ZDE0.025")
summary(contcont)
plot(contcont)


