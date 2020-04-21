## NIR and MIR spectrophotomwtry 
## Pablo Salazar 
## Universidad de Piura
## 21/01/2019

## Libraries----------
library(readxl) ## To open excel files
library(pls) ## To calculate PCR and LS
library(caret) ## it split the data randomly
library(tidyverse) ## For data manipulation
library(devtools) ## for ggbiplot
library(ggbiplot) ## for the PCA
library(factoextra) ## To get eigenvalues in a PCA
library(corrplot) ## for corrplot graphs
library(Hmisc) ## for rcorr function
library(reshape2) ## To use melt 
library(chillR) ## To calculate RPD
library(ggpubr) ## To use ggarrange
library(neuralnet) ## To use an artifitial neural network to predict
library(NeuralNetTools) ## To improve the graphical result
library(Metrics) ## to calculate rmse in the NN
library(qpcR) ## to cbind unequal colummns
library(dplyr) ## To run SummarySE
## Summarize function --------
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

## loading and analyzing  --------------

NIR <- read_excel("Datos/Atrib y Refl PERU.xlsx")
names(NIR)
colnames(NIR)[1]<-"number" ## to put a real name to v1

# statistical data distribution before analisis 

means <- NIR[,1:8] %>%
  pivot_longer(cols = -number, names_to = "nutrient" , values_to = "cc") %>%
  dplyr::group_by(nutrient)%>%
  dplyr::summarise(Mean = mean(cc, na.rm = TRUE),
            SD = sd(cc, na.rm = TRUE), 
            Min = min(cc, na.rm = TRUE), 
            Max = max(cc, na.rm = TRUE))

write.csv(means, "Results/Nutrient means.csv")

# wavelength analysis 

NIR[, 9:178]<-log10(100/NIR[, 9:178]) ## To transform the data to absorbance

NIR.pca <- NIR[, 9:178] 
numbers <- NIR[, 1] ## just renamed the categorical variable
ir.pca <- prcomp(NIR.pca,
                 center = FALSE,
                 scale = TRUE
) ## PCA operation with center and scale active to reduce skewness effect
print(ir.pca)
plot(ir.pca,type="barplot")
summary(ir.pca)

test2 <- get_pca_ind(ir.pca) ## Individual contribution to the PCA
test2$contrib
test2$coord
test2$cos2
coords<-cbind(numbers,test2$coord)
histogram(coords$Dim.2)
subset(coords,coords$Dim.2< -2) ## here we notice we have to erase 77 and 78
## PCA to explore the data 
pca <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
              ellipse = TRUE, 
              circle = FALSE, repel =TRUE)
pca <-pca+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            axis.text=element_text(size=16),
            legend.text=element_text(size=16),
            legend.direction = 'horizontal',
            legend.position="top",
            axis.title.x=element_text(size=16),
            axis.title.y=element_text(size=16))
pca

# Now we know "number" 77 and 78 have to go

NIR<-NIR[-c(77,78),,drop=F] ## So we erase it
## To split the data in training and test set -----------

set.seed(123) ## To fix the random matrix

training.samples <- NIR$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- NIR[training.samples, ] ## To make the train data
test.data <- NIR[-training.samples, ] ## To make the test data
train.data1<- as.matrix(train.data[,-c(1:8),drop=F]) ## To subset NIR variable
test.data1<- as.matrix(test.data[,-c(1:8),drop=F]) ## To subset NIR variable

## To inspect the data distribution ---------------

train.data2<-melt(train.data[,-c(2:8),drop=F],id="number")
train.data2$variable<-as.numeric(paste(train.data2$variable))
means<-summarySE(train.data2, measurevar="value",groupvars="variable")

Temas<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             panel.border =element_rect(size = 0.5, linetype="solid",fill=NA, colour ="black"),
             axis.text =element_text(size=16),
             axis.title=element_text(size=16,face="bold"),
             plot.title=element_text(size=20,face="bold", hjust = 0.5),
             legend.text=element_text(size=20),
             legend.title=element_text(size=16),
             legend.position="top")
scaleNIR<-scale_x_continuous(breaks = seq(400,2100,by=400),
                             limits = c(300,2200) ,
                             labels=paste0(seq(400,2100,by=400)),
                             expand=c(0.01,0))
a<-ggplot(data=train.data2)+
  aes(x=variable)+
  geom_line(aes(y=value, group=factor(number)))+
  labs(y="Absorbance", x= "nm")+
  Temas + scaleNIR
a

aa<-ggplot(data=means)+
  aes(x=variable)+
  geom_line(aes(y=value),size=1)+
  geom_line(aes(y=value-sd),linetype="dashed")+
  geom_line(aes(y=value+sd),linetype="dashed")+
  labs(y="Absorbance", x= "Wavelength (nm)")+
  Temas + scaleNIR
aa

## Phosphorous analysis ----------
## run train model with plsr 

model <- plsr(P ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
               segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

b<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Phosphorus", x= "nm")+
  Temas + scaleNIR
b

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$P, pred))
rpd_cal<-RPD(pred,train.data$P) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$P),
  Rsquare = caret::R2(pred, train.data$P))
## calibration graph

c<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
c

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$P, prediction))
rpd_test<-RPD(prediction,test.data$P) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$P),
  Rsquare = caret::R2(prediction, test.data$P))

## prediction accuracy graph with external data 

cc<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
cc

# combine vector of results
Phosphorus <- c("P","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

NIR_P_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Nitrogen analysis --------------
## run train model with plsr 

model <- plsr(N ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

d<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Nitrogen", x= "nm")+
  Temas + scaleNIR
d

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$N, pred))
rpd_cal<-RPD(pred,train.data$N) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$N, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$N, na.rm = TRUE))
## calibration graph

e<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
e

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$N, prediction))
rpd_test<-RPD(prediction,test.data$N) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$N),
  Rsquare = caret::R2(prediction, test.data$N))

## prediction accuracy graph with external data 

ee<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
ee
# combine vector of results 
Nitrogen <- c("N","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
NIR_N_reg_pls<-qpcR:::cbind.na(calib,ex_calib)

## Carbon analysis --------------
## run train model with plsr 

model <- plsr(C ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

f<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Carbon", x= "nm")+
  Temas + scaleNIR
f

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$C, pred))
rpd_cal<-RPD(pred,train.data$C) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$C, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$C, na.rm = TRUE))
## calibration graph

g<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
g

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$C, prediction))
rpd_test<-RPD(prediction,test.data$C) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$C),
  Rsquare = caret::R2(prediction, test.data$C))

## prediction accuracy graph with external data 

gg<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,4)) + scale_y_continuous(limits=c(0,4))
gg
# combine vector of results 

Carbon <- c("C","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
NIR_C_reg_pls<-qpcR:::cbind.na(calib,ex_calib)

## Iron analysis --------------
## run train model with plsr 

model <- plsr(Fe ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:5, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:5, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

h<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Iron", x= "nm")+
  Temas + scaleNIR
h

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Fe, pred))
rpd_cal<-RPD(pred,train.data$Fe) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Fe, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Fe, na.rm = TRUE))
## calibration graph

i<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
i

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Fe, prediction))
rpd_test<-RPD(prediction,test.data$Fe) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Fe),
  Rsquare = caret::R2(prediction, test.data$Fe))

## prediction accuracy graph with external data 

ii<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
ii
# combined vector of results 

Iron <- c("Fe","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

## Manganese analysis --------------
## run train model with plsr 

model <- plsr(Mn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=1:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

j<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Manganese", x= "nm")+
  Temas + scaleNIR
j

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Mn, pred))
rpd_cal<-RPD(pred,train.data$Mn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Mn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Mn, na.rm = TRUE))
## calibration graph

k<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
k

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Mn, prediction))
rpd_test<-RPD(prediction,test.data$Mn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Mn),
  Rsquare = caret::R2(prediction, test.data$Mn))

## prediction accuracy graph with external data 

kk<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
kk
# combine vector of results 

Manganese <- c("Mn","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
NIR_Mn_reg_pls<-qpcR:::cbind.na(calib,ex_calib)

## Cupper analysis --------------
## run train model with plsr 

model <- plsr(Cu ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)

vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

l<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Cupper", x= "nm")+
  Temas + scaleNIR
l

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Cu, pred))
rpd_cal<-RPD(pred,train.data$Cu) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Cu, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Cu, na.rm = TRUE))
## calibration graph

m<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
m

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Cu, prediction))
rpd_test<-RPD(prediction,test.data$Cu) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Cu),
  Rsquare = caret::R2(prediction, test.data$Cu))

## prediction accuracy graph with external data 

mm<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
mm
# combine vector of results 

Cupper <- c("Cu","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)


## Zinc analysis --------------
## run train model with plsr 

model <- plsr(Zn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

n<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Zinc", x= "nm")+
  Temas + scaleNIR
n

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Zn, pred))
rpd_cal<-RPD(pred,train.data$Zn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Zn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Zn, na.rm = TRUE))
## calibration graph

o<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
o

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Zn, prediction))
rpd_test<-RPD(prediction,test.data$Zn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Zn),
  Rsquare = caret::R2(prediction, test.data$Zn))

## prediction accuracy graph with external data 

oo<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
oo
# combine vector of results 

Zinc <- c("Zn","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)


## Vis-NIR final figures and tables only run once ------------

Vis_NIR<-as.data.frame(rbind(Phosphorus,Nitrogen,Carbon,Iron,Manganese,Cupper,Zinc))
colnames(Vis_NIR)[1:8]<-c("Nutrient","N? of Comp","RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
write.csv(Vis_NIR,"results/Vis NIR Results.csv")

VisNIR<- ggarrange(f,d,b,j,h,l,n, ncol=2, nrow=4, align="v",
              labels = c("A)","B)","C)","D)","E)","F)","G)"),
              font.label = list(size = 22, face = "bold"))
tiff("results/VIP_NIR.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
VisNIR
dev.off()

VisNIR_cal<- ggarrange(g,e,c,k,i,m,o, ncol=2, nrow=4, align="v",
                   labels = c("A)","B)","C)","D)","E)","F)","G)"),
                   font.label = list(size = 22, face = "bold"))
tiff("results/VisNIR_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
VisNIR_cal
dev.off()

VisNIR_val<- ggarrange(gg,ee,cc,kk,ii,mm,oo, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("results/VisNIR_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
VisNIR_val
dev.off()

## Neural network data setting for VisNIR --------

NIR1<-NIR[,-1,drop=F]
colnames(NIR1)[8:177] <- paste("V", colnames(NIR1[8:177]), sep = "")

max = apply(NIR1 , 2 , max,na.rm=TRUE)
min = apply(NIR1, 2 , min,na.rm=TRUE)
scaleNIR1 = as.data.frame(scale(NIR1, center = min, scale = max - min))

set.seed(123) ## To fix the random matrix

training.samples <- scaleNIR1$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- scaleNIR1[training.samples, ] ## To make the train data
test.data <- scaleNIR1[-training.samples, ] ## To make the test data
train.data1  <- NIR1[training.samples, ] ## To make the train data without scale
test.data1<- NIR1[-training.samples, ] ## To make the test data without scale

## buidling the formula for the phosphorus Neural network --------------
name<-names(train.data[,-c(2:7),drop=F])
formule <- as.formula(paste("P ~", paste(name[!name %in% "P"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_P<-b + geom_vline(xintercept=lines,lty="dashed",size=1)

b1<-VIP_nn_P + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"P_predicted"
P_predicted= (nn_result$P_predicted * (max(NIR1$P) - min(NIR1$P))) + min(NIR1$P)

nn_train_P<-as.data.frame(cbind(P_predicted,train.data1$P))
r_P_val<-cor(nn_result$P_predicted,train.data1$P)
r2cal<-r_P_val*r_P_val
nn_c<-ggplot(data=nn_train_P)+
  aes(x=V2)+
  geom_point(aes(y=P_predicted))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,P_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_val*r_P_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_c

nn_RMSE_P_cal<-rmse(nn_result$P_predicted,train.data1$P)
nn_rpd_cal<-RPD(nn_result$P_predicted,train.data1$P)

## external validation of the neural network with testdata
nn_val_P<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_P= (nn_val_P * (max(NIR1$P) - min(NIR1$P))) + min(NIR1$P)

nn_val_P<-as.data.frame(cbind(nn_val_P,test.data1$P))
r_P_cal<-cor(nn_val_P$V1,nn_val_P$V2)
r2val<-r_P_cal*r_P_cal

nn_cc<-ggplot(data=nn_val_P)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_cal*r_P_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_cc

nn_RMSE_P_val<-rmse(nn_val_P$V1,nn_val_P$V2)
nn_rpd_val<-RPD(nn_val_P$V1,nn_val_P$V2)

nn_phosphorus<-c(nn_RMSE_P_cal,r2cal,nn_rpd_cal,nn_RMSE_P_val,r2val,nn_rpd_val)
NIR_P_reg_NN<-qpcR:::cbind.na(nn_train_P,nn_val_P)
NIR_P_reg<-qpcR:::cbind.na(NIR_P_reg_pls,NIR_P_reg_NN)
colnames(NIR_P_reg)<-c("train_P", "pls_P_train","test_P","pls_P_test","NN_P_Train","train_P2","NN_P_Test", "test_P2")

## buidling the formula for the Nitrogen Neural network --------------
name<-names(train.data[,-c(1,3:7),drop=F])
train.data.NC <- train.data %>% drop_na()
train.data1.NC<-train.data1 %>% drop_na()
formule <- as.formula(paste("N ~", paste(name[!name %in% "N"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data.NC,hidden=c(5,3),linear.output = T,rep=110)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_N<-d + geom_vline(xintercept=lines,lty="dashed",size=1)

d1<-VIP_nn_N + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"N_predicted"
N_predicted= (nn_result$N_predicted * (max(NIR1$N,na.rm=TRUE) - min(NIR1$N,na.rm=TRUE))) + min(NIR1$N,na.rm=TRUE)

nn_train_N<-as.data.frame(cbind(N_predicted,train.data1.NC$N))
r_N_val<-cor(nn_result$N_predicted,train.data1.NC$N)
r2cal<-r_N_val*r_N_val
nn_e<-ggplot(data=nn_train_N)+
  aes(x=V2)+
  geom_point(aes(y=N_predicted))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,N_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_val*r_N_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_e

nn_RMSE_N_cal<-rmse(nn_result$N_predicted,train.data1.NC$N)
nn_rpd_cal<-RPD(nn_result$N_predicted,train.data1.NC$N)


## external validation of the neural network with testdata
nn_val_N<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_N= (nn_val_N * (max(NIR1$N,na.rm=TRUE) - min(NIR1$N,na.rm=TRUE))) + min(NIR1$N,na.rm=TRUE)

nn_val_N<-as.data.frame(cbind(nn_val_N,test.data1$N))
r_N_cal<-cor(nn_val_N$V1,nn_val_N$V2)
r2val<-r_N_cal*r_N_cal

nn_ee<-ggplot(data=nn_val_N)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_cal*r_N_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_ee

nn_RMSE_N_val<-rmse(nn_val_N$V1,nn_val_N$V2)
nn_rpd_val<-RPD(nn_val_N$V1,nn_val_N$V2)

nn_Nitrogen<-c(nn_RMSE_N_cal,r2cal,nn_rpd_cal,nn_RMSE_N_val,r2val,nn_rpd_val)
NIR_N_reg_NN<-qpcR:::cbind.na(nn_train_N,nn_val_N)
NIR_N_reg<-qpcR:::cbind.na(NIR_N_reg_pls,NIR_N_reg_NN)
colnames(NIR_N_reg)<-c("train_N","pls_N_train","test_N","pls_N_test","NN_N_Train","train_N2","NN_N_Test", "test_N2")

## buidling the formula for the Carbon Neural network --------------
name<-names(train.data[,-c(1:2,4:7),drop=F])
formule <- as.formula(paste("C ~", paste(name[!name %in% "C"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data.NC,hidden=c(5,3),linear.output = T,rep=110)
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_C<-f + geom_vline(xintercept=lines,lty="dashed",size=1)

f1<-VIP_nn_C + theme(axis.title=element_blank(), axis.text.x=element_blank())
## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"C_predicted"
C_predicted= (nn_result$C_predicted * (max(NIR1$C,na.rm=TRUE) - min(NIR1$C,na.rm=TRUE))) + min(NIR1$C,na.rm=TRUE)
nn_train_C<-as.data.frame(cbind(C_predicted,train.data1.NC$C))
r_C_val<-cor(nn_result$C_predicted,train.data1.NC$C)
r2cal<-r_C_val*r_C_val
nn_g<-ggplot(data=nn_train_C)+
  aes(x=V2)+
  geom_point(aes(y=C_predicted))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,C_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_val*r_C_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_g

nn_RMSE_C_cal<-rmse(nn_result$C_predicted,train.data1.NC$C)
nn_rpd_cal<-RPD(nn_result$C_predicted,train.data1.NC$C)


## external validation of the neural network with testdata
nn_val_C<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_C= (nn_val_C * (max(NIR1$C,na.rm=TRUE) - min(NIR1$C,na.rm=TRUE))) + min(NIR1$C,na.rm=TRUE)

nn_val_C<-as.data.frame(cbind(nn_val_C,test.data1$C))
r_C_cal<-cor(nn_val_C$V1,nn_val_C$V2)
r2val<-r_C_cal*r_C_cal

nn_gg<-ggplot(data=nn_val_C)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_cal*r_C_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_gg

nn_RMSE_C_val<-rmse(nn_val_C$V1,nn_val_C$V2)
nn_rpd_val<-RPD(nn_val_C$V1,nn_val_C$V2)

nn_Carbon<-c(nn_RMSE_C_cal,r2cal,nn_rpd_cal,nn_RMSE_C_val,r2val,nn_rpd_val)
NIR_C_reg_NN<-qpcR:::cbind.na(nn_train_C,nn_val_C)
NIR_C_reg<-qpcR:::cbind.na(NIR_C_reg_pls,NIR_C_reg_NN)
colnames(NIR_C_reg)<-c("train_C", "pls_C_train","test_C","pls_C_test","NN_C_Train","train_C2","NN_C_Test", "test_C2")

## buidling the formula for the Iron Neural network --------------
name<-names(train.data[,-c(1:3,5:7),drop=F])
formule <- as.formula(paste("Fe ~", paste(name[!name %in% "Fe"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Fe<-h + geom_vline(xintercept=lines,lty="dashed",size=1)

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Fe_predicted"
Fe_predicted= (nn_result$Fe_predicted * (max(NIR1$Fe) - min(NIR1$Fe))) + min(NIR1$Fe)

nn_train_Fe<-as.data.frame(cbind(Fe_predicted,train.data1$Fe))
r_Fe_val<-cor(nn_result$Fe_predicted,train.data1$Fe)
r2cal<-r_Fe_val*r_Fe_val
nn_i<-ggplot(data=nn_train_Fe)+
  aes(x=V2)+
  geom_point(aes(y=Fe_predicted))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,Fe_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_val*r_Fe_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_i

nn_RMSE_Fe_cal<-rmse(nn_result$Fe_predicted,train.data1$Fe)
nn_rpd_cal<-RPD(nn_result$Fe_predicted,train.data1$Fe)

## external validation of the neural network with testdata
nn_val_Fe<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_Fe= (nn_val_Fe * (max(NIR1$Fe) - min(NIR1$Fe))) + min(NIR1$Fe)

nn_val_Fe<-as.data.frame(cbind(nn_val_Fe,test.data1$Fe))
r_Fe_cal<-cor(nn_val_Fe$V1,nn_val_Fe$V2)
r2val<-r_Fe_cal*r_Fe_cal

nn_ii<-ggplot(data=nn_val_Fe)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_cal*r_Fe_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_ii

nn_RMSE_Fe_val<-rmse(nn_val_Fe$V1,nn_val_Fe$V2)
nn_rpd_val<-RPD(nn_val_Fe$V1,nn_val_Fe$V2)

nn_iron<-c(nn_RMSE_Fe_cal,r2cal,nn_rpd_cal,nn_RMSE_Fe_val,r2val,nn_rpd_val)

## buidling the formula for the Manganese Neural network --------------
name<-names(train.data[,-c(1:4,6:7),drop=F])
formule <- as.formula(paste("Mn ~", paste(name[!name %in% "Mn"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Mn<-j + geom_vline(xintercept=lines,lty="dashed",size=1)

j1<-VIP_nn_Mn + theme(axis.title=element_blank())
## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Mn_predicted"
Mn_predicted= (nn_result$Mn_predicted * (max(NIR1$Mn) - min(NIR1$Mn))) + min(NIR1$Mn)

nn_train_Mn<-as.data.frame(cbind(Mn_predicted,train.data1$Mn))
r_Mn_val<-cor(nn_result$Mn_predicted,train.data1$Mn)
r2cal<-r_Mn_val*r_Mn_val
nn_k<-ggplot(data=nn_train_Mn)+
  aes(x=V2)+
  geom_point(aes(y=Mn_predicted))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,Mn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_val*r_Mn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_k

nn_RMSE_Mn_cal<-rmse(nn_result$Mn_predicted,train.data1$Mn)
nn_rpd_cal<-RPD(nn_result$Mn_predicted,train.data1$Mn)

## external validation of the neural network with testdata
nn_val_Mn<-predict(NIR_nn,test.data, rep =which.min(NIR_nn$result.matrix[1,]))
nn_val_Mn= (nn_val_Mn * (max(NIR1$Mn) - min(NIR1$Mn))) + min(NIR1$Mn)

nn_val_Mn<-as.data.frame(cbind(nn_val_Mn,test.data1$Mn))
r_Mn_cal<-cor(nn_val_Mn$V1,nn_val_Mn$V2)
r2val<-r_Mn_cal*r_Mn_cal

nn_kk<-ggplot(data=nn_val_Mn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_cal*r_Mn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_kk

nn_RMSE_Mn_val<-rmse(nn_val_Mn$V1,nn_val_Mn$V2)
nn_rpd_val<-RPD(nn_val_Mn$V1,nn_val_Mn$V2)

nn_manganese<-c(nn_RMSE_Mn_cal,r2cal,nn_rpd_cal,nn_RMSE_Mn_val,r2val,nn_rpd_val)
NIR_Mn_reg_NN<-qpcR:::cbind.na(nn_train_Mn,nn_val_Mn)
NIR_Mn_reg<-qpcR:::cbind.na(NIR_Mn_reg_pls,NIR_Mn_reg_NN)
colnames(NIR_Mn_reg)<-c("train_Mn","pls_Mn_train","test_Mn","pls_Mn_test","NN_Mn_Train","train_Mn2","NN_Mn_Test", "test_Mn2")

## buidling the formula for the cupper Neural network --------------
name<-names(train.data[,-c(1:5,7),drop=F])
formule <- as.formula(paste("Cu ~", paste(name[!name %in% "Cu"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Cu<-l + geom_vline(xintercept=lines,lty="dashed",size=1)

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Cu_predicted"
Cu_predicted= (nn_result$Cu_predicted * (max(NIR1$Cu) - min(NIR1$Cu))) + min(NIR1$Cu)

nn_train_Cu<-as.data.frame(cbind(Cu_predicted,train.data1$Cu))
r_Cu_val<-cor(nn_result$Cu_predicted,train.data1$Cu)
r2cal<-r_Cu_val*r_Cu_val
nn_m<-ggplot(data=nn_train_Cu)+
  aes(x=V2)+
  geom_point(aes(y=Cu_predicted))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,Cu_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_val*r_Cu_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_m

nn_RMSE_Cu_cal<-rmse(nn_result$Cu_predicted,train.data1$Cu)
nn_rpd_cal<-RPD(nn_result$Cu_predicted,train.data1$Cu)

## external validation of the neural network with testdata
nn_val_Cu<-predict(NIR_nn,test.data, rep= which.min(NIR_nn$result.matrix[1,]))
nn_val_Cu= (nn_val_Cu * (max(NIR1$Cu) - min(NIR1$Cu))) + min(NIR1$Cu)

nn_val_Cu<-as.data.frame(cbind(nn_val_Cu,test.data1$Cu))
r_Cu_cal<-cor(nn_val_Cu$V1,nn_val_Cu$V2)
r2val<-r_Cu_cal*r_Cu_cal

nn_mm<-ggplot(data=nn_val_Cu)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_cal*r_Cu_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_mm

nn_RMSE_Cu_val<-rmse(nn_val_Cu$V1,nn_val_Cu$V2)
nn_rpd_val<-RPD(nn_val_Cu$V1,nn_val_Cu$V2)

nn_cupper<-c(nn_RMSE_Cu_cal,r2cal,nn_rpd_cal,nn_RMSE_Cu_val,r2val,nn_rpd_val)

## buidling the formula for the zinc Neural network --------------

name<-names(train.data[,-c(1:6),drop=F])
formule <- as.formula(paste("Zn ~", paste(name[!name %in% "Zn"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=80)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,166:170),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Zn<-n + geom_vline(xintercept=lines,lty="dashed",size=1)

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Zn_predicted"
Zn_predicted= (nn_result$Zn_predicted * (max(NIR1$Zn) - min(NIR1$Zn))) + min(NIR1$Zn)

nn_train_Zn<-as.data.frame(cbind(Zn_predicted,train.data1$Zn))
r_Zn_val<-cor(nn_result$Zn_predicted,train.data1$Zn)
r2cal<-r_Zn_val*r_Zn_val
nn_o<-ggplot(data=nn_train_Zn)+
  aes(x=V2)+
  geom_point(aes(y=Zn_predicted))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,Zn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_val*r_Zn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_o

nn_RMSE_Zn_cal<-rmse(nn_result$Zn_predicted,train.data1$Zn)
nn_rpd_cal<-RPD(nn_result$Zn_predicted,train.data1$Zn)

## external validation of the neural network with testdata
nn_val_Zn<-predict(NIR_nn,test.data, rep= which.min(NIR_nn$result.matrix[1,]))
nn_val_Zn= (nn_val_Zn * (max(NIR1$Zn) - min(NIR1$Zn))) + min(NIR1$Zn)

nn_val_Zn<-as.data.frame(cbind(nn_val_Zn,test.data1$Zn))
r_Zn_cal<-cor(nn_val_Zn$V1,nn_val_Zn$V2)
r2val<-r_Zn_cal*r_Zn_cal

nn_oo<-ggplot(data=nn_val_Zn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_cal*r_Zn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_oo

nn_RMSE_Zn_val<-rmse(nn_val_Zn$V1,nn_val_Zn$V2)
nn_rpd_val<-RPD(nn_val_Zn$V1,nn_val_Zn$V2)

nn_zinc<-c(nn_RMSE_Zn_cal,r2cal,nn_rpd_cal,nn_RMSE_Zn_val,r2val,nn_rpd_val)

## NN final figure and tables only run once ----------------

NN_VisNIR<-as.data.frame(rbind(nn_phosphorus,nn_Nitrogen,nn_Carbon,nn_iron,nn_manganese,nn_cupper,nn_zinc))
colnames(NN_VisNIR)<-c("RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
VisNIR2<-cbind(Vis_NIR,NN_VisNIR)
write.csv(VisNIR2,"results/PLS and NN VisNIR Results.csv")

NN_VisNIR_plot<- ggarrange(VIP_nn_C,VIP_nn_N,VIP_nn_P,VIP_nn_Mn,VIP_nn_Fe,VIP_nn_Cu,VIP_nn_Zn, ncol=2, nrow=4, align="v",
                   labels = c("A)","B)","C)","D)","E)","F)","G)"),
                   font.label = list(size = 22, face = "bold"))
tiff("results/NN_VIP_NIR.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_VisNIR_plot
dev.off()

NN_VisNIR_cal_plot<- ggarrange(nn_g,nn_e,nn_c,nn_k,nn_i,nn_m,nn_o, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("results/NN_VisNIR_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_VisNIR_cal_plot
dev.off()

NN_VisNIR_val_plot<- ggarrange(nn_gg,nn_ee,nn_cc,nn_kk,nn_ii,nn_mm,nn_oo, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("results/NN_VisNIR_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_VisNIR_val_plot
dev.off()
## Only Vis setting --------

set.seed(123) ## To fix the random matrix

training.samples <- NIR$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- NIR[training.samples, ] ## To make the train data
test.data <- NIR[-training.samples, ] ## To make the test data
train.data1<- as.matrix(train.data[,-c(1:137),drop=F]) ## To subset VIS variable
test.data1<- as.matrix(test.data[,-c(1:137),drop=F]) ## To subset VIS variable
names(test.data)

scaleVis<- scale_x_continuous(breaks = seq(400,800,by=100),
                   limits = c(350,850) ,
                   labels=paste0(seq(400,800,by=100)),
                   expand=c(0.01,0))
## Phosphorous analysis ----------
## run train model with plsr 

model <- plsr(P ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

b<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Phosphorus", x= "nm")+
  Temas + scaleVis
b

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$P, pred))
rpd_cal<-RPD(pred,train.data$P) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$P),
  Rsquare = caret::R2(pred, train.data$P))
## calibration graph

c<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
c

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$P, prediction))
rpd_test<-RPD(prediction,test.data$P) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$P),
  Rsquare = caret::R2(prediction, test.data$P))

## prediction accuracy graph with external data 

cc<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
cc

# combine vector of results
Phosphorus <- c("P","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
VIS_P_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Nitrogen analysis --------------
## run train model with plsr 

model <- plsr(N ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

d<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Nitrogen", x= "nm")+
  Temas + scaleVis
d

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$N, pred))
rpd_cal<-RPD(pred,train.data$N) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$N, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$N, na.rm = TRUE))
## calibration graph

e<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
e

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$N, prediction))
rpd_test<-RPD(prediction,test.data$N) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$N),
  Rsquare = caret::R2(prediction, test.data$N))

## prediction accuracy graph with external data 

ee<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
ee
# combine vector of results 
Nitrogen <- c("N","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
VIS_N_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Carbon analysis --------------
## run train model with plsr 

model <- plsr(C ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

f<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Carbon", x= "nm")+
  Temas + scaleVis
f

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$C, pred))
rpd_cal<-RPD(pred,train.data$C) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$C, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$C, na.rm = TRUE))
## calibration graph

g<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
g

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$C, prediction))
rpd_test<-RPD(prediction,test.data$C) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$C),
  Rsquare = caret::R2(prediction, test.data$C))

## prediction accuracy graph with external data 

gg<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,4)) + scale_y_continuous(limits=c(0,4))
gg
# combine vector of results 

Carbon <- c("C","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
VIS_C_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Iron analysis --------------
## run train model with plsr 

model <- plsr(Fe ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=3,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:3, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=1:3, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:2,4:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

h<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Iron", x= "nm")+
  Temas + scaleVis
h

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=3)
calib<-as.data.frame(cbind(train.data$Fe, pred))
rpd_cal<-RPD(pred,train.data$Fe) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Fe, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Fe, na.rm = TRUE))
## calibration graph

i<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
i

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=3)
ex_calib<-as.data.frame(cbind(test.data$Fe, prediction))
rpd_test<-RPD(prediction,test.data$Fe) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Fe),
  Rsquare = caret::R2(prediction, test.data$Fe))

## prediction accuracy graph with external data 

ii<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
ii
# combine vector of results 

Iron <- c("Fe","3",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

## Manganese analysis --------------
## run train model with plsr 

model <- plsr(Mn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

j<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Manganese", x= "nm")+
  Temas + scaleVis
j

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Mn, pred))
rpd_cal<-RPD(pred,train.data$Mn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Mn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Mn, na.rm = TRUE))
## calibration graph

k<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
k

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Mn, prediction))
rpd_test<-RPD(prediction,test.data$Mn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Mn),
  Rsquare = caret::R2(prediction, test.data$Mn))

## prediction accuracy graph with external data 

kk<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
kk
# combine vector of results 

Manganese <- c("Mn","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)
VIS_Mn_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Cupper analysis --------------
## run train model with plsr 

model <- plsr(Cu ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

l<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Cupper", x= "nm")+
  Temas + scaleVis
l

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Cu, pred))
rpd_cal<-RPD(pred,train.data$Cu) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Cu, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Cu, na.rm = TRUE))
## calibration graph

m<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
m

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Cu, prediction))
rpd_test<-RPD(prediction,test.data$Cu) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Cu),
  Rsquare = caret::R2(prediction, test.data$Cu))

## prediction accuracy graph with external data 

mm<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
mm
# combine vector of results 

Cupper <- c("Cu","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)


## Zinc analysis --------------
## run train model with plsr 

model <- plsr(Zn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

n<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Zinc", x= "nm")+
  Temas + scaleVis
 n

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Zn, pred))
rpd_cal<-RPD(pred,train.data$Zn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Zn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Zn, na.rm = TRUE))
## calibration graph

o<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
o

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Zn, prediction))
rpd_test<-RPD(prediction,test.data$Zn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Zn),
  Rsquare = caret::R2(prediction, test.data$Zn))

## prediction accuracy graph with external data 

oo<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
oo
# combine vector of results 

Zinc <- c("Zn","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)


## Vis final figures and tables only run once ------------

Vis<-as.data.frame(rbind(Phosphorus,Nitrogen,Carbon,Iron,Manganese,Cupper,Zinc))
colnames(Vis)[1:8]<-c("Nutrient","N? of Comp","RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
write.csv(Vis,"Results/Vis Results.csv")

Vis_plot<- ggarrange(f,d,b,j,h,l,n, ncol=2, nrow=4, align="v",
              labels = c("A)","B)","C)","D)","E)","F)","G)"),
              font.label = list(size = 22, face = "bold"))
tiff("Results/VIP_vis.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
Vis_plot
dev.off()

Vis_cal<- ggarrange(g,e,c,k,i,m,o, ncol=2, nrow=4, align="v",
                   labels = c("A)","B)","C)","D)","E)","F)","G)"),
                   font.label = list(size = 22, face = "bold"))
tiff("Results/Vis_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
Vis_cal
dev.off()

Vis_val<- ggarrange(gg,ee,cc,kk,ii,mm,oo, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("Results/Vis_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
Vis_val
dev.off()

## Neural network data setting for Vis --------

names(NIR1)
NIR1<-NIR[,-1,drop=F]
colnames(NIR1)[8:177] <- paste("V", colnames(NIR1[8:177]), sep = "")

max = apply(NIR1 , 2 , max,na.rm=TRUE)
min = apply(NIR1, 2 , min,na.rm=TRUE)
scaleNIR1 = as.data.frame(scale(NIR1, center = min, scale = max - min))

set.seed(123) ## To fix the random matrix

training.samples <- scaleNIR1$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- scaleNIR1[training.samples, ] ## To make the train data
test.data <- scaleNIR1[-training.samples, ] ## To make the test data
train.data1  <- NIR1[training.samples, ] ## To make the train data without scale
test.data1<- NIR1[-training.samples, ] ## To make the test data without scale

## buidling the formula for the phosphorus Neural network --------------
name<-names(train.data[,-c(2:136),drop=F])
formule <- as.formula(paste("P ~", paste(name[!name %in% "P"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=60)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_P<-b + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_P.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_P
dev.off()

b2<-VIP_nn_P + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"P_predicted"
P_predicted= (nn_result$P_predicted * (max(NIR1$P) - min(NIR1$P))) + min(NIR1$P)

nn_train_P<-as.data.frame(cbind(P_predicted,train.data1$P))
r_P_val<-cor(nn_result$P_predicted,train.data1$P)
r2cal<-r_P_val*r_P_val
nn_c<-ggplot(data=nn_train_P)+
  aes(x=V2)+
  geom_point(aes(y=P_predicted))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,P_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_val*r_P_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_c

tiff("VIS_nn_P_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_c
dev.off()

nn_RMSE_P_cal<-rmse(nn_result$P_predicted,train.data1$P)
nn_rpd_cal<-RPD(nn_result$P_predicted,train.data1$P)

## external validation of the neural network with testdata
nn_val_P<-predict(NIR_nn,test.data, rep= which.min(NIR_nn$result.matrix[1,]))
nn_val_P= (nn_val_P * (max(NIR1$P) - min(NIR1$P))) + min(NIR1$P)

nn_val_P<-as.data.frame(cbind(nn_val_P,test.data1$P))
r_P_cal<-cor(nn_val_P$V1,nn_val_P$V2)
r2val<-r_P_cal*r_P_cal

nn_cc<-ggplot(data=nn_val_P)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_cal*r_P_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_cc

tiff("VIS_nn_P_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_cc
dev.off()

nn_RMSE_P_val<-rmse(nn_val_P$V1,nn_val_P$V2)
nn_rpd_val<-RPD(nn_val_P$V1,nn_val_P$V2)

nn_phosphorus<-c(nn_RMSE_P_cal,r2cal,nn_rpd_cal,nn_RMSE_P_val,r2val,nn_rpd_val)
VIS_P_reg_NN<-qpcR:::cbind.na(nn_train_P,nn_val_P)
VIS_P_reg<-qpcR:::cbind.na(VIS_P_reg_pls,VIS_P_reg_NN)
colnames(VIS_P_reg)<-c("train_P", "pls_P_train","test_P","pls_P_test","NN_P_Train","train_P2","NN_P_Test", "test_P2")

## buidling the formula for the Nitrogen Neural network --------------
name<-names(train.data[,-c(1,3:136),drop=F])
train.data.NC<-train.data[-125,,drop=F]
train.data1.NC<-train.data1[-125,,drop=F]
formule <- as.formula(paste("N ~", paste(name[!name %in% "N"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data.NC,hidden=c(5,3),linear.output = T,rep=60)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_N<-d + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_N.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_N
dev.off()

d2<-VIP_nn_N + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"N_predicted"
N_predicted= (nn_result$N_predicted * (max(NIR1$N,na.rm=TRUE) - min(NIR1$N,na.rm=TRUE))) + min(NIR1$N,na.rm=TRUE)

nn_train_N<-as.data.frame(cbind(N_predicted,train.data1.NC$N))
r_N_val<-cor(nn_result$N_predicted,train.data1.NC$N)
r2cal<-r_N_val*r_N_val
nn_e<-ggplot(data=nn_train_N)+
  aes(x=V2)+
  geom_point(aes(y=N_predicted))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,N_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_val*r_N_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_e

tiff("VIS_nn_N_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_e
dev.off()

nn_RMSE_N_cal<-rmse(nn_result$N_predicted,train.data1.NC$N)
nn_rpd_cal<-RPD(nn_result$N_predicted,train.data1.NC$N)

## external validation of the neural network with testdata
nn_val_N<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_N= (nn_val_N * (max(NIR1$N,na.rm=TRUE) - min(NIR1$N,na.rm=TRUE))) + min(NIR1$N,na.rm=TRUE)

nn_val_N<-as.data.frame(cbind(nn_val_N,test.data1$N))
r_N_cal<-cor(nn_val_N$V1,nn_val_N$V2)
r2val<-r_N_cal*r_N_cal

nn_ee<-ggplot(data=nn_val_N)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_cal*r_N_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_ee

tiff("VIS_nn_N_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_ee
dev.off()

nn_RMSE_N_val<-rmse(nn_val_N$V1,nn_val_N$V2)
nn_rpd_val<-RPD(nn_val_N$V1,nn_val_N$V2)

nn_Nitrogen<-c(nn_RMSE_N_cal,r2cal,nn_rpd_cal,nn_RMSE_N_val,r2val,nn_rpd_val)
VIS_N_reg_NN<-qpcR:::cbind.na(nn_train_N,nn_val_N)
VIS_N_reg<-qpcR:::cbind.na(VIS_N_reg_pls,VIS_N_reg_NN)
colnames(VIS_N_reg)<-c("train_N","pls_N_train","test_N","pls_N_test","NN_N_Train","train_N2","NN_N_Test", "test_N2")

## buidling the formula for the Carbon Neural network --------------
name<-names(train.data[,-c(1:2,4:136),drop=F])
formule <- as.formula(paste("C ~", paste(name[!name %in% "C"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data.NC,hidden=c(5,3),linear.output = T,rep=60)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_C<-f + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_C.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_C
dev.off()

f2<-VIP_nn_C + theme(axis.title=element_blank(), axis.text.x=element_blank())
## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"C_predicted"
C_predicted= (nn_result$C_predicted * (max(NIR1$C,na.rm=TRUE) - min(NIR1$C,na.rm=TRUE))) + min(NIR1$C,na.rm=TRUE)
nn_train_C<-as.data.frame(cbind(C_predicted,train.data1.NC$C))
r_C_val<-cor(nn_result$C_predicted,train.data1.NC$C)
r2cal<-r_C_val*r_C_val
nn_g<-ggplot(data=nn_train_C)+
  aes(x=V2)+
  geom_point(aes(y=C_predicted))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,C_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_val*r_C_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_g

tiff("VIS_nn_C_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_g
dev.off()

nn_RMSE_C_cal<-rmse(nn_result$C_predicted,train.data1.NC$C)
nn_rpd_cal<-RPD(nn_result$C_predicted,train.data1.NC$C)


## external validation of the neural network with testdata
nn_val_C<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_C= (nn_val_C * (max(NIR1$C,na.rm=TRUE) - min(NIR1$C,na.rm=TRUE))) + min(NIR1$C,na.rm=TRUE)

nn_val_C<-as.data.frame(cbind(nn_val_C,test.data1$C))
r_C_cal<-cor(nn_val_C$V1,nn_val_C$V2)
r2val<-r_C_cal*r_C_cal

nn_gg<-ggplot(data=nn_val_C)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_cal*r_C_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_gg

tiff("VIS_nn_C_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_gg
dev.off()

nn_RMSE_C_val<-rmse(nn_val_C$V1,nn_val_C$V2)
nn_rpd_val<-RPD(nn_val_C$V1,nn_val_C$V2)

nn_Carbon<-c(nn_RMSE_C_cal,r2cal,nn_rpd_cal,nn_RMSE_C_val,r2val,nn_rpd_val)
VIS_C_reg_NN<-qpcR:::cbind.na(nn_train_C,nn_val_C)
VIS_C_reg<-qpcR:::cbind.na(VIS_C_reg_pls,VIS_C_reg_NN)
colnames(VIS_C_reg)<-c("train_C", "pls_C_train","test_C","pls_C_test","NN_C_Train","train_C2","NN_C_Test", "test_C2")

## buidling the formula for the Iron Neural network --------------
name<-names(train.data[,-c(1:3,5:136),drop=F])
formule <- as.formula(paste("Fe ~", paste(name[!name %in% "Fe"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=5)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Fe<-h + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_Fe.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Fe
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Fe_predicted"
Fe_predicted= (nn_result$Fe_predicted * (max(NIR1$Fe) - min(NIR1$Fe))) + min(NIR1$Fe)

nn_train_Fe<-as.data.frame(cbind(Fe_predicted,train.data1$Fe))
r_Fe_val<-cor(nn_result$Fe_predicted,train.data1$Fe)
r2cal<-r_Fe_val*r_Fe_val
nn_i<-ggplot(data=nn_train_Fe)+
  aes(x=V2)+
  geom_point(aes(y=Fe_predicted))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,Fe_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_val*r_Fe_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_i

tiff("VIS_nn_Fe_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_i
dev.off()

nn_RMSE_Fe_cal<-rmse(nn_result$Fe_predicted,train.data1$Fe)
nn_rpd_cal<-RPD(nn_result$Fe_predicted,train.data1$Fe)

## external validation of the neural network with testdata
nn_val_Fe<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_Fe= (nn_val_Fe * (max(NIR1$Fe) - min(NIR1$Fe))) + min(NIR1$Fe)

nn_val_Fe<-as.data.frame(cbind(nn_val_Fe,test.data1$Fe))
r_Fe_cal<-cor(nn_val_Fe$V1,nn_val_Fe$V2)
r2val<-r_Fe_cal*r_Fe_cal

nn_ii<-ggplot(data=nn_val_Fe)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_cal*r_Fe_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_ii

tiff("VIS_nn_Fe_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_ii
dev.off()

nn_RMSE_Fe_val<-rmse(nn_val_Fe$V1,nn_val_Fe$V2)
nn_rpd_val<-RPD(nn_val_Fe$V1,nn_val_Fe$V2)

nn_iron<-c(nn_RMSE_Fe_cal,r2cal,nn_rpd_cal,nn_RMSE_Fe_val,r2val,nn_rpd_val)

## buidling the formula for the Manganese Neural network --------------
name<-names(train.data[,-c(1:4,6:136),drop=F])
formule <- as.formula(paste("Mn ~", paste(name[!name %in% "Mn"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=60)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Mn<-j + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_Mn.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Mn
dev.off()

j2<-VIP_nn_Mn + theme(axis.title=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Mn_predicted"
Mn_predicted= (nn_result$Mn_predicted * (max(NIR1$Mn) - min(NIR1$Mn))) + min(NIR1$Mn)

nn_train_Mn<-as.data.frame(cbind(Mn_predicted,train.data1$Mn))
r_Mn_val<-cor(nn_result$Mn_predicted,train.data1$Mn)
r2cal<-r_Mn_val*r_Mn_val
nn_k<-ggplot(data=nn_train_Mn)+
  aes(x=V2)+
  geom_point(aes(y=Mn_predicted))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,Mn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_val*r_Mn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_k

tiff("VIS_nn_Mn_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_k
dev.off()

nn_RMSE_Mn_cal<-rmse(nn_result$Mn_predicted,train.data1$Mn)
nn_rpd_cal<-RPD(nn_result$Mn_predicted,train.data1$Mn)

## external validation of the neural network with testdata
nn_val_Mn<-predict(NIR_nn,test.data, rep= which.min(NIR_nn$result.matrix[1,]))
nn_val_Mn= (nn_val_Mn * (max(NIR1$Mn) - min(NIR1$Mn))) + min(NIR1$Mn)

nn_val_Mn<-as.data.frame(cbind(nn_val_Mn,test.data1$Mn))
r_Mn_cal<-cor(nn_val_Mn$V1,nn_val_Mn$V2)
r2val<-r_Mn_cal*r_Mn_cal

nn_kk<-ggplot(data=nn_val_Mn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_cal*r_Mn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_kk

tiff("VIS_nn_Mn_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_kk
dev.off()

nn_RMSE_Mn_val<-rmse(nn_val_Mn$V1,nn_val_Mn$V2)
nn_rpd_val<-RPD(nn_val_Mn$V1,nn_val_Mn$V2)

nn_manganese<-c(nn_RMSE_Mn_cal,r2cal,nn_rpd_cal,nn_RMSE_Mn_val,r2val,nn_rpd_val)
VIS_Mn_reg_NN<-qpcR:::cbind.na(nn_train_Mn,nn_val_Mn)
VIS_Mn_reg<-qpcR:::cbind.na(VIS_Mn_reg_pls,VIS_Mn_reg_NN)
colnames(VIS_Mn_reg)<-c("train_Mn","pls_Mn_train","test_Mn","pls_Mn_test","NN_Mn_Train","train_Mn2","NN_Mn_Test", "test_Mn2")

## buidling the formula for the cupper Neural network --------------
name<-names(train.data[,-c(1:5,7:136),drop=F])
formule <- as.formula(paste("Cu ~", paste(name[!name %in% "Cu"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=60)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Cu<-l + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_Cu.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Cu
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Cu_predicted"
Cu_predicted= (nn_result$Cu_predicted * (max(NIR1$Cu) - min(NIR1$Cu))) + min(NIR1$Cu)

nn_train_Cu<-as.data.frame(cbind(Cu_predicted,train.data1$Cu))
r_Cu_val<-cor(nn_result$Cu_predicted,train.data1$Cu)
r2cal<-r_Cu_val*r_Cu_val
nn_m<-ggplot(data=nn_train_Cu)+
  aes(x=V2)+
  geom_point(aes(y=Cu_predicted))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,Cu_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_val*r_Cu_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_m

tiff("VIS_nn_Cu_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_m
dev.off()

nn_RMSE_Cu_cal<-rmse(nn_result$Cu_predicted,train.data1$Cu)
nn_rpd_cal<-RPD(nn_result$Cu_predicted,train.data1$Cu)

## external validation of the neural network with testdata
nn_val_Cu<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_Cu= (nn_val_Cu * (max(NIR1$Cu) - min(NIR1$Cu))) + min(NIR1$Cu)

nn_val_Cu<-as.data.frame(cbind(nn_val_Cu,test.data1$Cu))
r_Cu_cal<-cor(nn_val_Cu$V1,nn_val_Cu$V2)
r2val<-r_Cu_cal*r_Cu_cal

nn_mm<-ggplot(data=nn_val_Cu)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_cal*r_Cu_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_mm

tiff("VIS_nn_Cu_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_mm
dev.off()

nn_RMSE_Cu_val<-rmse(nn_val_Cu$V1,nn_val_Cu$V2)
nn_rpd_val<-RPD(nn_val_Cu$V1,nn_val_Cu$V2)

nn_cupper<-c(nn_RMSE_Cu_cal,r2cal,nn_rpd_cal,nn_RMSE_Cu_val,r2val,nn_rpd_val)

## buidling the formula for the zinc Neural network --------------

name<-names(train.data[,-c(1:6,8:136),drop=F])
formule <- as.formula(paste("Zn ~", paste(name[!name %in% "Zn"], collapse = " + ")))
## Running the nnet
NIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=5)
NIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(NIR_nn,rep="best")
nn_importance<-olden(NIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,37:41),,drop=F]
lines<-substr(rownames(nn_importance),2,5)
lines<-as.numeric(lines)
VIP_nn_Zn<-n + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("VIS_RI_nn_Zn.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Zn
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(NIR_nn$net.result[which.min(NIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Zn_predicted"
Zn_predicted= (nn_result$Zn_predicted * (max(NIR1$Zn) - min(NIR1$Zn))) + min(NIR1$Zn)

nn_train_Zn<-as.data.frame(cbind(Zn_predicted,train.data1$Zn))
r_Zn_val<-cor(nn_result$Zn_predicted,train.data1$Zn)
r2cal<-r_Zn_val*r_Zn_val
nn_o<-ggplot(data=nn_train_Zn)+
  aes(x=V2)+
  geom_point(aes(y=Zn_predicted))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,Zn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_val*r_Zn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_o

tiff("VIS_nn_Zn_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_o
dev.off()

nn_RMSE_Zn_cal<-rmse(nn_result$Zn_predicted,train.data1$Zn)
nn_rpd_cal<-RPD(nn_result$Zn_predicted,train.data1$Zn)

## external validation of the neural network with testdata
nn_val_Zn<-predict(NIR_nn,test.data, rep = which.min(NIR_nn$result.matrix[1,]))
nn_val_Zn= (nn_val_Zn * (max(NIR1$Zn) - min(NIR1$Zn))) + min(NIR1$Zn)

nn_val_Zn<-as.data.frame(cbind(nn_val_Zn,test.data1$Zn))
r_Zn_cal<-cor(nn_val_Zn$V1,nn_val_Zn$V2)
r2val<-r_Zn_cal*r_Zn_cal

nn_oo<-ggplot(data=nn_val_Zn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_cal*r_Zn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_oo

tiff("VIS_nn_Zn_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_oo
dev.off()

nn_RMSE_Zn_val<-rmse(nn_val_Zn$V1,nn_val_Zn$V2)
nn_rpd_val<-RPD(nn_val_Zn$V1,nn_val_Zn$V2)

nn_zinc<-c(nn_RMSE_Zn_cal,r2cal,nn_rpd_cal,nn_RMSE_Zn_val,r2val,nn_rpd_val)

## NN final figure and tables for VIS only run once ----------------

NN_Vis<-as.data.frame(rbind(nn_phosphorus,nn_Nitrogen,nn_Carbon,nn_iron,nn_manganese,nn_cupper,nn_zinc))
colnames(NN_Vis)<-c("RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
VisNIR2<-cbind(Vis,NN_Vis)
write.csv(VisNIR2,"Results/PLS and NN Vis Results.csv")

NN_Vis_plot<- ggarrange(VIP_nn_C,VIP_nn_N,VIP_nn_P,VIP_nn_Mn,VIP_nn_Fe,VIP_nn_Cu,VIP_nn_Zn, ncol=2, nrow=4, align="v",
                           labels = c("A)","B)","C)","D)","E)","F)","G)"),
                           font.label = list(size = 22, face = "bold"))
tiff("Results/NN_VIP_VIS.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_Vis_plot
dev.off()

NN_Vis_cal_plot<- ggarrange(nn_g,nn_e,nn_c,nn_k,nn_i,nn_m,nn_o, ncol=2, nrow=4, align="v",
                               labels = c("A)","B)","C)","D)","E)","F)","G)"),
                               font.label = list(size = 22, face = "bold"))
tiff("Results/NN_Vis_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_Vis_cal_plot
dev.off()

NN_Vis_val_plot<- ggarrange(nn_gg,nn_ee,nn_cc,nn_kk,nn_ii,nn_mm,nn_oo, ncol=2, nrow=4, align="v",
                               labels = c("A)","B)","C)","D)","E)","F)","G)"),
                               font.label = list(size = 22, face = "bold"))
tiff("Results/NN_Vis_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_Vis_val_plot
dev.off()
## MIR loading and setting --------

MIR <- read_excel("E:/Dropbox/Fincyt/NIR/Datos/Mario_MIR_absorbancia.xlsx")
names(MIR)
wn<-colnames(MIR)[9:1951]
wn<-as.numeric(wn)
wn<-1/wn*10000000
colnames(MIR)[9:1951]<-wn
MIR<-MIR[,c(1:8,352:1873),drop=F]
MIR<-MIR[-c(171,183),,drop=F] ## To eliminate a really weird sample and a sample without N data
colnames(MIR)[1]<-"number" ## to put a real name to v1
MIR.pca <- MIR[, 9:1530] 
ir.pca <- prcomp(MIR.pca,
                 center = FALSE,
                 scale = TRUE
) ## PCA operation with center and scale active to reduce skewness effect
print(ir.pca)
numbers <- MIR[, 1] ## just renamed the categorical variable
plot(ir.pca,type="barplot")
summary(ir.pca)

test2 <- get_pca_ind(ir.pca) ## Individual contribution to the PCA
test2$contrib
test2$coord
test2$cos2
coords<-cbind(numbers,test2$coord)
histogram(coords$Dim.7)
subset(coords,coords$Dim.7< -1) ## we also exclude 78 just to keep 
## PCA to explore the data 
pca_MIR <- ggbiplot(ir.pca, obs.scale = 1, var.scale = 1, 
                ellipse = TRUE, 
                circle = FALSE, repel =TRUE)
pca_MIR <-pca_MIR+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              axis.text=element_text(size=16),
              legend.text=element_text(size=16),
              legend.direction = 'horizontal',
              legend.position="top",
              axis.title.x=element_text(size=16),
              axis.title.y=element_text(size=16))
pca_MIR

## MIR - To split the data in training and test set ------------------

set.seed(123) ## To fix the random matrix

training.samples <- MIR$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- MIR[training.samples, ] ## To make the train data
test.data <- MIR[-training.samples, ] ## To make the test data
train.data1<- as.matrix(train.data[,-c(1:8),drop=F]) ## To subset MIR variable
test.data1<- as.matrix(test.data[,-c(1:8),drop=F]) ## To subset MIR variable
names(test.data)

## To inspect the data distribution ---------------

train.data2<-melt(train.data[,-c(2:8),drop=F],id="number")
train.data2$variable<-as.numeric(paste(train.data2$variable))
means<-summarySE(train.data2, measurevar="value",groupvars="variable")

Temas<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  ## delete background
             panel.background = element_blank(), axis.line = element_line(colour = "black"),
             panel.border =element_rect(size = 0.5, linetype="solid",fill=NA, colour ="black"),
             axis.text =element_text(size=16),
             axis.title=element_text(size=16,face="bold"),
             plot.title=element_text(size=20,face="bold", hjust = 0.5),
             legend.text=element_text(size=20),
             legend.title=element_text(size=16),
             legend.position="top")
scaleMIR <- scale_x_continuous(breaks = seq(3000,25000,by=4000),
                               limits = c(2500,26000) ,
                               labels=paste0(seq(3000,25000,by=4000)),
                               expand=c(0.01,0))
p<-ggplot(data=train.data2)+
  aes(x=variable)+
  geom_line(aes(y=value, group=factor(number)))+
  labs(y="Absorbance", x= "nm")+
  Temas + scaleMIR
p

pp<-ggplot(data=means)+
  aes(x=variable)+
  geom_line(aes(y=value),size=1)+
  geom_line(aes(y=value-sd),linetype="dashed")+
  geom_line(aes(y=value+sd),linetype="dashed")+
  labs(y="Absorbance", x= "Wavelength (nm)")+
  Temas+ scaleMIR
  
pp

## plotting Absorbance of VisNIR and MIR


absorbance<- ggarrange(aa,pp, ncol=2, nrow=1, align="h",
                labels = c("A)","B)"),
                font.label = list(size = 22, face = "bold"))
tiff("absorbance.tiff", width = 12, height = 4.5, units = 'in', res = 300, compression = 'lzw')
absorbance
dev.off()

## Phosphorous analysis ----------
## run train model with plsr 

model <- plsr(P ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=4:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

b<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Phosphorus", x= "nm")+
  Temas + scaleMIR
b

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$P, pred))
rpd_cal<-RPD(pred,train.data$P) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$P),
  Rsquare = caret::R2(pred, train.data$P))
## calibration graph

c<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
c

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$P, prediction))
rpd_test<-RPD(prediction,test.data$P) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$P),
  Rsquare = caret::R2(prediction, test.data$P))

## prediction accuracy graph with external data 

cc<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed P")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
cc

# combine vector of results
Phosphorus <- c("P","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

MIR_P_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Nitrogen analysis --------------
## run train model with plsr 

model <- plsr(N ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

d<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Nitrogen", x= "nm")+
  Temas + scaleMIR
d

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$N, pred))
rpd_cal<-RPD(pred,train.data$N) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$N, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$N, na.rm = TRUE))
## calibration graph

e<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
e

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$N, prediction))
rpd_test<-RPD(prediction,test.data$N) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$N),
  Rsquare = caret::R2(prediction, test.data$N))

## prediction accuracy graph with external data 

ee<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed N")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
ee
# combine vector of results 
Nitrogen <- c("N","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

MIR_N_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Carbon analysis --------------
## run train model with plsr 

model <- plsr(C ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

f<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Carbon", x= "nm")+
  Temas + scaleMIR
f

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$C, pred))
rpd_cal<-RPD(pred,train.data$C) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$C, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$C, na.rm = TRUE))
## calibration graph

g<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
g

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$C, prediction))
rpd_test<-RPD(prediction,test.data$C) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$C),
  Rsquare = caret::R2(prediction, test.data$C))

## prediction accuracy graph with external data 

gg<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed C")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,4)) + scale_y_continuous(limits=c(0,4))
gg
# combine vector of results 

Carbon <- c("C","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

MIR_C_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Iron analysis --------------
## run train model with plsr 

model <- plsr(Fe ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:5, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

h<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Iron", x= "nm")+
  Temas + scaleMIR
h

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Fe, pred))
rpd_cal<-RPD(pred,train.data$Fe) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Fe, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Fe, na.rm = TRUE))
## calibration graph

i<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
i

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Fe, prediction))
rpd_test<-RPD(prediction,test.data$Fe) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Fe),
  Rsquare = caret::R2(prediction, test.data$Fe))

## prediction accuracy graph with external data 

ii<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Fe")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
ii
# combine vector of results 

Iron <- c("Fe","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

## Manganese analysis --------------
## run train model with plsr 

model <- plsr(Mn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=7,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:6,8:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

j<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Manganese", x= "nm")+
  Temas + scaleMIR
j

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=7)
calib<-as.data.frame(cbind(train.data$Mn, pred))
rpd_cal<-RPD(pred,train.data$Mn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Mn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Mn, na.rm = TRUE))
## calibration graph

k<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
k

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=7)
ex_calib<-as.data.frame(cbind(test.data$Mn, prediction))
rpd_test<-RPD(prediction,test.data$Mn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Mn),
  Rsquare = caret::R2(prediction, test.data$Mn))

## prediction accuracy graph with external data 

kk<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Mn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=15, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
kk
# combine vector of results 

Manganese <- c("Mn","7",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

MIR_Mn_reg_pls<-qpcR:::cbind.na(calib,ex_calib)
## Cupper analysis --------------
## run train model with plsr 

model <- plsr(Cu ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=5,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:5, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=4:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:4,6:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

l<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Cupper", x= "nm")+
  Temas + scaleMIR
l

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=5)
calib<-as.data.frame(cbind(train.data$Cu, pred))
rpd_cal<-RPD(pred,train.data$Cu) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Cu, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Cu, na.rm = TRUE))
## calibration graph

m<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
m

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=5)
ex_calib<-as.data.frame(cbind(test.data$Cu, prediction))
rpd_test<-RPD(prediction,test.data$Cu) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Cu),
  Rsquare = caret::R2(prediction, test.data$Cu))

## prediction accuracy graph with external data 

mm<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Cu")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
mm
# combine vector of results 

Cupper <- c("Cu","5",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)

## Zinc analysis --------------
## run train model with plsr 

model <- plsr(Zn ~ train.data1, ncomp = 10, data = train.data, scale = TRUE,
              segments= 60, validation = "CV",method="oscorespls")
ncomp.onesigma <- selectNcomp(model, method = "onesigma", plot = TRUE) #choosing 
#the model with fewest components that is still less than one standard error away 
#from the overall best model
ncomp.onesigma <- selectNcomp(model, method = "randomization", plot = TRUE) #tests 
#whether adding a new component is beneficial at all

plot(model,ncomp=6,asp = 1 , line =TRUE)
plot(model, "loadings", comps = 1:7, legendpos = "topleft",labels = "numbers", xlab = "nm")
plot(model, plottype = "coef", ncomp=3:7, legendpos = "bottomleft", labels = "numbers", xlab = "nm")

# To see the PLS correlation 
model$loadings
model_cor <- rcorr(as.matrix(model$loadings))
model_cor
r<-model_cor$r
p<-model_cor$P

corrplot(model_cor$r, method = "ellipse",p.mat=model_cor$P, insig ="blank", type = "lower"
         , diag=FALSE, cl.cex=1, tl.cex=0.7)
## To calculate the variable importance for projection (VIP)
vip<-VIP(model)
summary(vip)
vip<- vip[-c(1:5,7:10),,drop=F]
vip2<-melt(vip)
vip2$Var2<-as.numeric(paste(vip2$Var2))

n<-ggplot(data=vip2)+
  aes(x=Var2)+
  geom_line(aes(y=value))+
  labs(y="VIP Zinc", x= "nm")+
  Temas + scaleMIR
n

## To check the accuracy of the model with the train data
pred = predict(model, train.data1, ncomp=6)
calib<-as.data.frame(cbind(train.data$Zn, pred))
rpd_cal<-RPD(pred,train.data$Zn) ## To calculate the residual prediction deviation (RDP)

perfomance_cal <- data.frame(
  RMSE = caret::RMSE(pred, train.data$Zn, na.rm = TRUE),
  Rsquare = caret::R2(pred, train.data$Zn, na.rm = TRUE))
## calibration graph

o<-ggplot(data=calib)+
  aes(x=V1)+
  geom_point(aes(y=pred))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,pred), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_cal[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
o

# Make predictions with external data
prediction = predict(model, newdata=test.data1, ncomp=6)
ex_calib<-as.data.frame(cbind(test.data$Zn, prediction))
rpd_test<-RPD(prediction,test.data$Zn) ## To calculate the residual prediction deviation (RDP)

# Model performance metrics fo external calibration
perfomance_test <- data.frame(
  RMSE = caret::RMSE(prediction, test.data$Zn),
  Rsquare = caret::R2(prediction, test.data$Zn))

## prediction accuracy graph with external data 

oo<-ggplot(data=ex_calib)+
  aes(x=V1)+
  geom_point(aes(y=prediction))+
  labs(y="PLS Predicted", x= "Observed Zn")+
  geom_smooth(aes(V1,prediction), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(perfomance_test[1,2],2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
oo
# combine vector of results 

Zinc <- c("Zn","6",perfomance_cal[1,1],perfomance_cal[1,2],rpd_cal,perfomance_test[1,1],perfomance_test[1,2],rpd_test)


## MIR final figures and tables only run once ------------

MIR_R<-as.data.frame(rbind(Phosphorus,Nitrogen,Carbon,Iron,Manganese,Cupper,Zinc))
colnames(MIR_R)[1:8]<-c("Nutrient","N? of Comp","RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
write.csv(MIR_R,"Results/MIR Results.csv")

MIR_VIP<- ggarrange(f,d,b,j,h,l,n, ncol=2, nrow=4, align="v",
                   labels = c("A)","B)","C)","D)","E)","F)","G)"),
                   font.label = list(size = 22, face = "bold"))
tiff("Results/MIR_VIP.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
MIR_VIP
dev.off()

MIR_cal<- ggarrange(g,e,c,k,i,m,o, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("Results/MIR_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
MIR_cal
dev.off()

MIR_val<- ggarrange(gg,ee,cc,kk,ii,mm,oo, ncol=2, nrow=4, align="v",
                       labels = c("A)","B)","C)","D)","E)","F)","G)"),
                       font.label = list(size = 22, face = "bold"))
tiff("Results/MIR_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
MIR_val
dev.off()

### Neural network data setting for MIR --------
names(MIR)
MIR1<-MIR[,-1,drop=F]
colnames(MIR1)[8:1529] <- paste("V", colnames(MIR1[8:1529]), sep = "")

max = apply(MIR1 , 2 , max,na.rm=TRUE)
min = apply(MIR1, 2 , min,na.rm=TRUE)
scaleMIR1 = as.data.frame(scale(MIR1, center = min, scale = max - min))

set.seed(123) ## To fix the random matrix

training.samples <- scaleMIR1$P %>% 
  createDataPartition(p = 0.7, list = FALSE)
train.data  <- scaleMIR1[training.samples, ] ## To make the train data
test.data <- scaleMIR1[-training.samples, ] ## To make the test data
train.data1  <- MIR1[training.samples, ] ## To make the train data without scale
test.data1<- MIR1[-training.samples, ] ## To make the test data without scale

## buidling the formula for the phosphorus Neural network --------------
name<-names(train.data[,-c(2:7),drop=F])
formule <- as.formula(paste("P ~", paste(name[!name %in% "P"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_P<-b + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_P.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_P
dev.off()

b3<-VIP_nn_P + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"P_predicted"
P_predicted= (nn_result$P_predicted * (max(MIR1$P) - min(MIR1$P))) + min(MIR1$P)

nn_train_P<-as.data.frame(cbind(P_predicted,train.data1$P))
r_P_val<-cor(nn_result$P_predicted,train.data1$P)
r2cal<-r_P_val*r_P_val
nn_c<-ggplot(data=nn_train_P)+
  aes(x=V2)+
  geom_point(aes(y=P_predicted))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,P_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_val*r_P_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_c

tiff("MIR_nn_P_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_c
dev.off()

nn_RMSE_P_cal<-rmse(nn_result$P_predicted,train.data1$P)
nn_rpd_cal<-RPD(nn_result$P_predicted,train.data1$P)

## external validation of the neural network with testdata
nn_val_P<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_P= (nn_val_P * (max(MIR1$P) - min(MIR1$P))) + min(MIR1$P)

nn_val_P<-as.data.frame(cbind(nn_val_P,test.data1$P))
r_P_cal<-cor(nn_val_P$V1,nn_val_P$V2)
r2val<-r_P_cal*r_P_cal

nn_cc<-ggplot(data=nn_val_P)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_P_cal*r_P_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_cc

tiff("MIR_nn_P_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_cc
dev.off()

nn_RMSE_P_val<-rmse(nn_val_P$V1,nn_val_P$V2)
nn_rpd_val<-RPD(nn_val_P$V1,nn_val_P$V2)

nn_phosphorus<-c(nn_RMSE_P_cal,r2cal,nn_rpd_cal,nn_RMSE_P_val,r2val,nn_rpd_val)
MIR_P_reg_NN<-qpcR:::cbind.na(nn_train_P,nn_val_P)
MIR_P_reg<-qpcR:::cbind.na(MIR_P_reg_pls,MIR_P_reg_NN)
colnames(MIR_P_reg)<-c("train_P", "pls_P_train","test_P","pls_P_test","NN_P_Train","train_P2","NN_P_Test", "test_P2")
## buidling the formula for the Nitrogen Neural network --------------
name<-names(train.data[,-c(1,3:7),drop=F])
formule <- as.formula(paste("N ~", paste(name[!name %in% "N"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_N<-d + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_N.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_N
dev.off()

d3<-VIP_nn_N + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"N_predicted"
N_predicted= (nn_result$N_predicted * (max(MIR1$N,na.rm=TRUE) - min(MIR1$N,na.rm=TRUE))) + min(MIR1$N,na.rm=TRUE)

nn_train_N<-as.data.frame(cbind(N_predicted,train.data1$N))
r_N_val<-cor(nn_result$N_predicted,train.data1$N)
r2cal<-r_N_val*r_N_val
nn_e<-ggplot(data=nn_train_N)+
  aes(x=V2)+
  geom_point(aes(y=N_predicted))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,N_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_val*r_N_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_e

tiff("MIR_nn_N_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_e
dev.off()

nn_RMSE_N_cal<-rmse(nn_result$N_predicted,train.data1$N)
nn_rpd_cal<-RPD(nn_result$N_predicted,train.data1$N)


## external validation of the neural network with testdata
nn_val_N<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_N= (nn_val_N * (max(MIR1$N,na.rm=TRUE) - min(MIR1$N,na.rm=TRUE))) + min(MIR1$N,na.rm=TRUE)

nn_val_N<-as.data.frame(cbind(nn_val_N,test.data1$N))
r_N_cal<-cor(nn_val_N$V1,nn_val_N$V2)
r2val<-r_N_cal*r_N_cal

nn_ee<-ggplot(data=nn_val_N)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=0.4, y=0.1), label=paste("R^2 ==",round(r_N_cal*r_N_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
nn_ee

tiff("MIR_nn_N_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_ee
dev.off()

nn_RMSE_N_val<-rmse(nn_val_N$V1,nn_val_N$V2)
nn_rpd_val<-RPD(nn_val_N$V1,nn_val_N$V2)

nn_Nitrogen<-c(nn_RMSE_N_cal,r2cal,nn_rpd_cal,nn_RMSE_N_val,r2val,nn_rpd_val)
MIR_N_reg_NN<-qpcR:::cbind.na(nn_train_N,nn_val_N)
MIR_N_reg<-qpcR:::cbind.na(MIR_N_reg_pls,MIR_N_reg_NN)
colnames(MIR_N_reg)<-c("train_N","pls_N_train","test_N","pls_N_test","NN_N_Train","train_N2","NN_N_Test", "test_N2")

## buidling the formula for the Carbon Neural network --------------
name<-names(train.data[,-c(1:2,4:7),drop=F])
formule <- as.formula(paste("C ~", paste(name[!name %in% "C"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(NIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_C<-f + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_C.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_C
dev.off()

f3<-VIP_nn_C + theme(axis.title=element_blank(), axis.text.x=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"C_predicted"
C_predicted= (nn_result$C_predicted * (max(MIR1$C,na.rm=TRUE) - min(MIR1$C,na.rm=TRUE))) + min(MIR1$C,na.rm=TRUE)
nn_train_C<-as.data.frame(cbind(C_predicted,train.data1$C))
r_C_val<-cor(nn_result$C_predicted,train.data1$C)
r2cal<-r_C_val*r_C_val
nn_g<-ggplot(data=nn_train_C)+
  aes(x=V2)+
  geom_point(aes(y=C_predicted))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,C_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_val*r_C_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_g

tiff("MIR_nn_C_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_g
dev.off()

nn_RMSE_C_cal<-rmse(nn_result$C_predicted,train.data1$C)
nn_rpd_cal<-RPD(nn_result$C_predicted,train.data1$C)


## external validation of the neural network with testdata
nn_val_C<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_C= (nn_val_C * (max(MIR1$C,na.rm=TRUE) - min(MIR1$C,na.rm=TRUE))) + min(MIR1$C,na.rm=TRUE)

nn_val_C<-as.data.frame(cbind(nn_val_C,test.data1$C))
r_C_cal<-cor(nn_val_C$V1,nn_val_C$V2)
r2val<-r_C_cal*r_C_cal

nn_gg<-ggplot(data=nn_val_C)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=4, y=1), label=paste("R^2 ==",round(r_C_cal*r_C_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
nn_gg

tiff("MIR_nn_C_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_gg
dev.off()

nn_RMSE_C_val<-rmse(nn_val_C$V1,nn_val_C$V2)
nn_rpd_val<-RPD(nn_val_C$V1,nn_val_C$V2)

nn_Carbon<-c(nn_RMSE_C_cal,r2cal,nn_rpd_cal,nn_RMSE_C_val,r2val,nn_rpd_val)
MIR_C_reg_NN<-qpcR:::cbind.na(nn_train_C,nn_val_C)
MIR_C_reg<-qpcR:::cbind.na(MIR_C_reg_pls,MIR_C_reg_NN)
colnames(MIR_C_reg)<-c("train_C", "pls_C_train","test_C","pls_C_test","NN_C_Train","train_C2","NN_C_Test", "test_C2")

## buidling the formula for the Iron Neural network --------------
name<-names(train.data[,-c(1:3,5:7),drop=F])
formule <- as.formula(paste("Fe ~", paste(name[!name %in% "Fe"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(MIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_Fe<-h + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_Fe.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Fe
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Fe_predicted"
Fe_predicted= (nn_result$Fe_predicted * (max(MIR1$Fe) - min(MIR1$Fe))) + min(MIR1$Fe)

nn_train_Fe<-as.data.frame(cbind(Fe_predicted,train.data1$Fe))
r_Fe_val<-cor(nn_result$Fe_predicted,train.data1$Fe)
r2cal<-r_Fe_val*r_Fe_val
nn_i<-ggplot(data=nn_train_Fe)+
  aes(x=V2)+
  geom_point(aes(y=Fe_predicted))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,Fe_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_val*r_Fe_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_i

tiff("MIR_nn_Fe_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_i
dev.off()

nn_RMSE_Fe_cal<-rmse(nn_result$Fe_predicted,train.data1$Fe)
nn_rpd_cal<-RPD(nn_result$Fe_predicted,train.data1$Fe)

## external validation of the neural network with testdata
nn_val_Fe<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_Fe= (nn_val_Fe * (max(MIR1$Fe) - min(MIR1$Fe))) + min(MIR1$Fe)

nn_val_Fe<-as.data.frame(cbind(nn_val_Fe,test.data1$Fe))
r_Fe_cal<-cor(nn_val_Fe$V1,nn_val_Fe$V2)
r2val<-r_Fe_cal*r_Fe_cal

nn_ii<-ggplot(data=nn_val_Fe)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Fe")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=40, y=10), label=paste("R^2 ==",round(r_Fe_cal*r_Fe_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,50)) + scale_y_continuous(limits=c(0,50))
nn_ii

tiff("MIR_nn_Fe_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_ii
dev.off()

nn_RMSE_Fe_val<-rmse(nn_val_Fe$V1,nn_val_Fe$V2)
nn_rpd_val<-RPD(nn_val_Fe$V1,nn_val_Fe$V2)

nn_iron<-c(nn_RMSE_Fe_cal,r2cal,nn_rpd_cal,nn_RMSE_Fe_val,r2val,nn_rpd_val)

## buidling the formula for the Manganese Neural network --------------
name<-names(train.data[,-c(1:4,6:7),drop=F])
formule <- as.formula(paste("Mn ~", paste(name[!name %in% "Mn"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(MIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_Mn<-j + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_Mn.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Mn
dev.off()

j3<-VIP_nn_Mn + theme(axis.title=element_blank())

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Mn_predicted"
Mn_predicted= (nn_result$Mn_predicted * (max(MIR1$Mn) - min(MIR1$Mn))) + min(MIR1$Mn)

nn_train_Mn<-as.data.frame(cbind(Mn_predicted,train.data1$Mn))
r_Mn_val<-cor(nn_result$Mn_predicted,train.data1$Mn)
r2cal<-r_Mn_val*r_Mn_val
nn_k<-ggplot(data=nn_train_Mn)+
  aes(x=V2)+
  geom_point(aes(y=Mn_predicted))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,Mn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_val*r_Mn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_k

tiff("MIR_nn_Mn_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_k
dev.off()

nn_RMSE_Mn_cal<-rmse(nn_result$Mn_predicted,train.data1$Mn)
nn_rpd_cal<-RPD(nn_result$Mn_predicted,train.data1$Mn)

## external validation of the neural network with testdata
nn_val_Mn<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_Mn= (nn_val_Mn * (max(MIR1$Mn) - min(MIR1$Mn))) + min(MIR1$Mn)

nn_val_Mn<-as.data.frame(cbind(nn_val_Mn,test.data1$Mn))
r_Mn_cal<-cor(nn_val_Mn$V1,nn_val_Mn$V2)
r2val<-r_Mn_cal*r_Mn_cal

nn_kk<-ggplot(data=nn_val_Mn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=20, y=5), label=paste("R^2 ==",round(r_Mn_cal*r_Mn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
nn_kk

tiff("MIR_nn_Mn_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_kk
dev.off()

nn_RMSE_Mn_val<-rmse(nn_val_Mn$V1,nn_val_Mn$V2)
nn_rpd_val<-RPD(nn_val_Mn$V1,nn_val_Mn$V2)

nn_manganese<-c(nn_RMSE_Mn_cal,r2cal,nn_rpd_cal,nn_RMSE_Mn_val,r2val,nn_rpd_val)
MIR_Mn_reg_NN<-qpcR:::cbind.na(nn_train_Mn,nn_val_Mn)
MIR_Mn_reg<-qpcR:::cbind.na(MIR_Mn_reg_pls,MIR_Mn_reg_NN)
colnames(MIR_Mn_reg)<-c("train_Mn","pls_Mn_train","test_Mn","pls_Mn_test","NN_Mn_Train","train_Mn2","NN_Mn_Test", "test_Mn2")

## buidling the formula for the cupper Neural network --------------
name<-names(train.data[,-c(1:5,7),drop=F])
formule <- as.formula(paste("Cu ~", paste(name[!name %in% "Cu"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=110)
MIR_nn$result.matrix[1,]
#plot(MIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_Cu<-l + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_Cu.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Cu
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Cu_predicted"
Cu_predicted= (nn_result$Cu_predicted * (max(MIR1$Cu) - min(MIR1$Cu))) + min(MIR1$Cu)

nn_train_Cu<-as.data.frame(cbind(Cu_predicted,train.data1$Cu))
r_Cu_val<-cor(nn_result$Cu_predicted,train.data1$Cu)
r2cal<-r_Cu_val*r_Cu_val
nn_m<-ggplot(data=nn_train_Cu)+
  aes(x=V2)+
  geom_point(aes(y=Cu_predicted))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,Cu_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_val*r_Cu_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_m

tiff("MIR_nn_Cu_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_m
dev.off()

nn_RMSE_Cu_cal<-rmse(nn_result$Cu_predicted,train.data1$Cu)
nn_rpd_cal<-RPD(nn_result$Cu_predicted,train.data1$Cu)

## external validation of the neural network with testdata
nn_val_Cu<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_Cu= (nn_val_Cu * (max(MIR1$Cu) - min(MIR1$Cu))) + min(MIR1$Cu)

nn_val_Cu<-as.data.frame(cbind(nn_val_Cu,test.data1$Cu))
r_Cu_cal<-cor(nn_val_Cu$V1,nn_val_Cu$V2)
r2val<-r_Cu_cal*r_Cu_cal

nn_mm<-ggplot(data=nn_val_Cu)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Cu")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=25, y=5), label=paste("R^2 ==",round(r_Cu_cal*r_Cu_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,40)) + scale_y_continuous(limits=c(0,40))
nn_mm

tiff("MIR_nn_Cu_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_mm
dev.off()

nn_RMSE_Cu_val<-rmse(nn_val_Cu$V1,nn_val_Cu$V2)
nn_rpd_val<-RPD(nn_val_Cu$V1,nn_val_Cu$V2)

nn_cupper<-c(nn_RMSE_Cu_cal,r2cal,nn_rpd_cal,nn_RMSE_Cu_val,r2val,nn_rpd_val)

## buidling the formula for the zinc Neural network --------------

name<-names(train.data[,-c(1:6),drop=F])
formule <- as.formula(paste("Zn ~", paste(name[!name %in% "Zn"], collapse = " + ")))
## Running the nnet
MIR_nn<-neuralnet(formule,data=train.data,hidden=c(5,3),linear.output = T,rep=80)
MIR_nn$result.matrix[1,]
#plot(MIR_nn)
olden(MIR_nn,rep="best")
nn_importance<-olden(MIR_nn,bar_plot=FALSE,rep="best")
nn_importance<-cbind(wavelength=rownames(nn_importance),nn_importance)
nn_importance <- nn_importance[with(nn_importance,order(-importance)),]
nn_importance<-nn_importance[c(1:5,1518:1522),,drop=F]
lines<-substr(rownames(nn_importance),2,16)
lines<-as.numeric(lines)
VIP_nn_Zn<-n + geom_vline(xintercept=lines,lty="dashed",size=1)

tiff("MIR_RI_nn_Zn.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
VIP_nn_Zn
dev.off()

## extracting the predicted results from the traindata
nn_result<-as.data.frame(MIR_nn$net.result[which.min(MIR_nn$result.matrix[1,])])
colnames(nn_result)<-"Zn_predicted"
Zn_predicted= (nn_result$Zn_predicted * (max(MIR1$Zn) - min(MIR1$Zn))) + min(MIR1$Zn)

nn_train_Zn<-as.data.frame(cbind(Zn_predicted,train.data1$Zn))
r_Zn_val<-cor(nn_result$Zn_predicted,train.data1$Zn)
r2cal<-r_Zn_val*r_Zn_val
nn_o<-ggplot(data=nn_train_Zn)+
  aes(x=V2)+
  geom_point(aes(y=Zn_predicted))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,Zn_predicted), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_val*r_Zn_val,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_o

tiff("MIR_nn_Zn_cal.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_o
dev.off()

nn_RMSE_Zn_cal<-rmse(nn_result$Zn_predicted,train.data1$Zn)
nn_rpd_cal<-RPD(nn_result$Zn_predicted,train.data1$Zn)

## external validation of the neural network with testdata
nn_val_Zn<-predict(MIR_nn,test.data, rep = which.min(MIR_nn$result.matrix[1,]))
nn_val_Zn= (nn_val_Zn * (max(MIR1$Zn) - min(MIR1$Zn))) + min(MIR1$Zn)

nn_val_Zn<-as.data.frame(cbind(nn_val_Zn,test.data1$Zn))
r_Zn_cal<-cor(nn_val_Zn$V1,nn_val_Zn$V2)
r2val<-r_Zn_cal*r_Zn_cal

nn_oo<-ggplot(data=nn_val_Zn)+
  aes(x=V2)+
  geom_point(aes(y=V1))+
  labs(y="nn Predicted", x= "Observed Zn")+
  geom_smooth(aes(V2,V1), method=lm, se=FALSE)+
  geom_text(aes(x=3, y=1), label=paste("R^2 ==",round(r_Zn_cal*r_Zn_cal,2)),size=6, parse =TRUE) +
  Temas+
  scale_x_continuous(limits=c(0,6)) + scale_y_continuous(limits=c(0,6))
nn_oo

tiff("MIR_nn_Zn_val.tiff", width = 6, height = 4.5, units = 'in', res = 300, compression = 'lzw')
nn_oo
dev.off()

nn_RMSE_Zn_val<-rmse(nn_val_Zn$V1,nn_val_Zn$V2)
nn_rpd_val<-RPD(nn_val_Zn$V1,nn_val_Zn$V2)

nn_zinc<-c(nn_RMSE_Zn_cal,r2cal,nn_rpd_cal,nn_RMSE_Zn_val,r2val,nn_rpd_val)

## NN final figure and tables only run once ----------------

NN_MIR<-as.data.frame(rbind(nn_phosphorus,nn_Nitrogen,nn_Carbon,nn_iron,nn_manganese,nn_cupper,nn_zinc))
colnames(NN_MIR)<-c("RMSE cal","R square cal","RPD cal","RMSE val","R square val","RPD val")
MIR_R2<-cbind(MIR_R,NN_MIR)
write.csv(MIR_R2,"Results/PLS and NN MIR Results.csv")

NN_MIR_plot<- ggarrange(VIP_nn_C,VIP_nn_N,VIP_nn_P,VIP_nn_Mn,VIP_nn_Fe,VIP_nn_Cu,VIP_nn_Zn, ncol=2, nrow=4, align="v",
                           labels = c("A)","B)","C)","D)","E)","F)","G)"),
                           font.label = list(size = 22, face = "bold"))
tiff("Results/NN_VIP_MIR.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_MIR_plot
dev.off()

NN_MIR_cal_plot<- ggarrange(nn_g,nn_e,nn_c,nn_k,nn_i,nn_m,nn_o, ncol=2, nrow=4, align="v",
                               labels = c("A)","B)","C)","D)","E)","F)","G)"),
                               font.label = list(size = 22, face = "bold"))
tiff("Results/NN_MIR_cal.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_MIR_cal_plot
dev.off()

NN_MIR_val_plot<- ggarrange(nn_gg,nn_ee,nn_cc,nn_kk,nn_ii,nn_mm,nn_oo, ncol=2, nrow=4, align="v",
                               labels = c("A)","B)","C)","D)","E)","F)","G)"),
                               font.label = list(size = 22, face = "bold"))
tiff("Results/NN_MIR_val.tiff", width = 12, height = 18, units = 'in', res = 300, compression = 'lzw')
NN_MIR_val_plot
dev.off()

## Final figures for the paper ----------
coef_PLS= cor(MIR_P_reg$train_P,MIR_P_reg$pls_P_train)
coef_nn=cor(MIR_P_reg$train_P,MIR_P_reg$NN_P_Train)

gg3I<-ggplot(data=MIR_P_reg)+
  aes(x=train_P)+
  geom_point(aes(y=pls_P_train), color = "red")+
  geom_point(aes(y=NN_P_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(train_P,pls_P_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_P,NN_P_Train), method=lm, se=FALSE, color = "blue")+
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse = TRUE, color= "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color= "blue") +
  geom_text(aes(x=1.5, y=4), label="Phosphorous",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3I

coef_PLS= cor(MIR_Mn_reg$train_P,MIR_Mn_reg$pls_Mn_train)
coef_nn=cor(MIR_Mn_reg$train_P,MIR_Mn_reg$NN_Mn_Train)

gg3L<-ggplot(data=MIR_Mn_reg)+
  aes(x=train_Mn)+
  geom_point(aes(y=pls_Mn_train), color = "red")+
  geom_point(aes(y=NN_Mn_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(train_Mn,pls_Mn_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_Mn,NN_Mn_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=20, y=5), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=20, y=3), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=7, y=20), label="Mangenese",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
gg3L

coef_PLS= cor(MIR_N_reg$train_N,MIR_N_reg$pls_N_train)
coef_nn=cor(MIR_N_reg$train_N,MIR_N_reg$NN_N_Train)

gg3F<-ggplot(data=MIR_N_reg)+
  aes(x=train_N)+
  geom_point(aes(y=pls_N_train), color = "red")+
  geom_point(aes(y=NN_N_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(train_N,pls_N_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_N,NN_N_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=0.4, y=0.1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=0.4, y=0.05), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=0.1, y=0.4), label="Nitrogen",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank()) +
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
gg3F


coef_PLS= cor(MIR_C_reg$train_C,MIR_C_reg$pls_C_train)
coef_nn=cor(MIR_C_reg$train_C,MIR_C_reg$NN_C_Train)

gg3c<-ggplot(data=MIR_C_reg)+
  aes(x=train_C)+
  geom_point(aes(y=pls_C_train), color = "red")+
  geom_point(aes(y=NN_C_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(train_C,pls_C_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_C,NN_C_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=1, y=4), label="Carbon",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3c

coef_PLS= cor(NIR_P_reg$train_P,NIR_P_reg$pls_P_train)
coef_nn=cor(NIR_P_reg$train_P,NIR_P_reg$NN_P_Train)

gg3H<-ggplot(data=NIR_P_reg)+
  aes(x=train_P)+
  geom_point(aes(y=pls_P_train), color = "red")+
  geom_point(aes(y=NN_P_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(train_P,pls_P_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_P,NN_P_Train), method=lm, se=FALSE, color = "blue")+
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse = TRUE, color= "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color= "blue") +
  geom_text(aes(x=1.5, y=4), label="Phosphorous",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3H

coef_PLS= cor(NIR_Mn_reg$train_Mn,NIR_Mn_reg$pls_Mn_train)
coef_nn=cor(NIR_Mn_reg$train_Mn,NIR_Mn_reg$NN_Mn_Train)

gg3K<-ggplot(data=NIR_Mn_reg)+
  aes(x=train_Mn)+
  geom_point(aes(y=pls_Mn_train), color = "red")+
  geom_point(aes(y=NN_Mn_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(train_Mn,pls_Mn_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_Mn,NN_Mn_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=20, y=5), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=20, y=3), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=7, y=20), label="Mangenese",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
gg3K

coef_PLS= cor(NIR_N_reg$train_N,NIR_N_reg$pls_N_train,use="complete.obs")
coef_nn=cor(NIR_N_reg$train_N,NIR_N_reg$NN_N_Train,use="complete.obs")

gg3E<-ggplot(data=NIR_N_reg)+
  aes(x=train_N)+
  geom_point(aes(y=pls_N_train), color = "red")+
  geom_point(aes(y=NN_N_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(train_N,pls_N_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_N,NN_N_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=0.4, y=0.1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=0.4, y=0.05), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=0.1, y=0.4), label="Nitrogen",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank())+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
gg3E

names(NIR_C_reg)
coef_PLS= cor(NIR_C_reg$train_C,NIR_C_reg$pls_C_train,use="complete.obs")
coef_nn=cor(NIR_C_reg$train_C,NIR_C_reg$NN_C_Train,use="complete.obs")

gg3B<-ggplot(data=NIR_C_reg)+
  aes(x=train_C)+
  geom_point(aes(y=pls_C_train), color = "red")+
  geom_point(aes(y=NN_C_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(train_C,pls_C_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_C,NN_C_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=1, y=4), label="Carbon",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank(), axis.text.y=element_blank()) +
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3B

coef_PLS= cor(VIS_P_reg$train_P,VIS_P_reg$pls_P_train)
coef_nn=cor(VIS_P_reg$train_P,VIS_P_reg$NN_P_Train)

gg3G<-ggplot(data=VIS_P_reg)+
  aes(x=train_P)+
  geom_point(aes(y=pls_P_train), color = "red")+
  geom_point(aes(y=NN_P_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed P")+
  geom_smooth(aes(train_P,pls_P_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_P,NN_P_Train), method=lm, se=FALSE, color = "blue")+
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse = TRUE, color= "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color= "blue") +
  geom_text(aes(x=1.5, y=4), label="Phosphorous",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank())+
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3G

coef_PLS= cor(VIS_Mn_reg$train_Mn,VIS_Mn_reg$pls_Mn_train)
coef_nn=cor(VIS_Mn_reg$train_Mn,VIS_Mn_reg$NN_Mn_Train)

gg3J<-ggplot(data=VIS_Mn_reg)+
  aes(x=train_Mn)+
  geom_point(aes(y=pls_Mn_train), color = "red")+
  geom_point(aes(y=NN_Mn_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed Mn")+
  geom_smooth(aes(train_Mn,pls_Mn_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_Mn,NN_Mn_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=20, y=5), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=20, y=3), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=7, y=20), label="Mangenese",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank())+
  scale_x_continuous(limits=c(0,30)) + scale_y_continuous(limits=c(0,30))
gg3J

coef_PLS= cor(VIS_N_reg$train_N,VIS_N_reg$pls_N_train,use="complete.obs")
coef_nn=cor(VIS_N_reg$train_N,VIS_N_reg$NN_N_Train,use="complete.obs")

gg3D<-ggplot(data=VIS_N_reg)+
  aes(x=train_N)+
  geom_point(aes(y=pls_N_train), color = "red")+
  geom_point(aes(y=NN_N_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed N")+
  geom_smooth(aes(train_N,pls_N_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_N,NN_N_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=0.4, y=0.1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=0.4, y=0.05), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=0.1, y=0.4), label="Nitrogen",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank())+
  scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
gg3D

coef_PLS= cor(VIS_C_reg$train_C,VIS_C_reg$pls_C_train,use="complete.obs")
coef_nn=cor(VIS_C_reg$train_C,VIS_C_reg$NN_C_Train,use="complete.obs")

gg3A<-ggplot(data=VIS_C_reg)+
  aes(x=train_C)+
  geom_point(aes(y=pls_C_train), color = "red")+
  geom_point(aes(y=NN_C_Train), color = "blue") +
  labs(y="nn Predicted", x= "Observed C")+
  geom_smooth(aes(train_C,pls_C_train), method=lm, se=FALSE, color = "red")+
  geom_smooth(aes(train_C,NN_C_Train), method=lm, se=FALSE, color = "blue")+ 
  geom_text(aes(x=4, y=1), label=paste("PLS~R^2 ==",round(coef_PLS*coef_PLS,2)),size=6, parse =TRUE, color = "red") +
  geom_text(aes(x=4, y=0.5), label=paste("NN~R^2 ==",round(coef_nn*coef_nn,2)),size=6, parse =TRUE, color = "blue") +
  geom_text(aes(x=1, y=4), label="Carbon",size=10, color = "black", fontface="bold") +
  Temas + theme(axis.title=element_blank()) +
  scale_x_continuous(limits=c(0,5)) + scale_y_continuous(limits=c(0,5))
gg3A

Figure3<- ggarrange(gg3A,gg3B,gg3c,gg3D,gg3E,gg3F,gg3G,gg3H,gg3I,gg3J,gg3K,gg3L, ncol=3, nrow=4, align="v",
                        labels = c("A)","B)","C)","D)","E)","F)","G)","H)","I)","J)","K)","L)"),
                        font.label = list(size = 22, face = "bold"),
                    hjust=c(-2,-2,-2,-2,-2,-2,-2,-2,-3,-2,-2,-2),vjust=2)
tiff("results/Figure3.tiff", width = 18, height = 18, units = 'in', res = 300, compression = 'lzw')
annotate_figure(Figure3,
                top = text_grob("VIS                                          NIR                                          MIR", face = "bold", size = 32),
                bottom = text_grob("Observed data", 
                                   face = "bold", size = 32),
                left = text_grob("Predicted data", rot = 90,
                                 face = "bold", size = 32)
)
dev.off()


Figure2<- ggarrange(f2,f1,f3,d2,d1,d3,b2,b1,b3,j2,j1,j3, ncol=3, nrow=4, align="v",
                    labels = c("A)","B)","C)","D)","E)","F)","G)","H)","I)","J)","K)","L)"),
                    font.label = list(size = 22, face = "bold"),
                    hjust=c(-17,-16.5,-16,-16.5,-17.5,-18.5, -15.5,-16.8,-28,-20,-16.5,-18), vjust=2)

tiff("results/Figure2.tiff", width = 18, height = 18, units = 'in', res = 300, compression = 'lzw')
annotate_figure(Figure2,
                top = text_grob("VIS                                          NIR                                          MIR", face = "bold", size = 32),
                bottom = text_grob("Wavelength (nn) \n Variable of Importance in Prediction", 
                                   face = "bold", size = 32),
                left = text_grob("     Mangenese                   Phosphorous                 Nitrogen                     Carbon      ", rot = 90,
                                 face = "bold", size = 30)
                )
dev.off()
