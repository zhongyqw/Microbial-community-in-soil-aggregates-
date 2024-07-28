library(ggpubr)
library(reshape2)
library(ggsci)
library(devtools)
library(vegan)
library(ggplot2)
library("RColorBrewer")
library(dplyr)
library(reshape)

otu1 <- read.csv('D:\\Program Files (x86)\\RStudio\\example\\otutab.csv',header=T,row.names = 1)
otu1<-data.frame(t(otu1))
observed_species <- rowSums(otu1 > 0)
Shannon <- diversity(otu1, index = 'shannon', base = exp(1))
alpha<-data.frame(Shannon,observed_species)
write.csv(alpha, 'alpha_diversity.csv', quote = FALSE)


data=read.csv('alpha_diversity.csv',header=T,row.names = 1)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
GROUP <- select(GROUP,c(2,3))
data <- data[rownames(GROUP), ]
data <- cbind( GROUP,data)
colnames(data)[4] <- c("Richness")

data_melt <- melt(data, id.vars = c("Stage","Group"), variable.name = "index",value.name = "y") 
colnames(data_melt)[c(3,4)] <- c("index","y")
data_melt$Group[data_melt$Group=="Silt-clay particles"]<-"Silt_clay particles"
data_melt$Group <- factor(data_melt$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt_clay particles"),labels= c("Total","LMA","SMA","MIA","SCP"))
data_melt_1<-data_melt[which(data_melt$index=="Shannon"),]

data=read.csv('alpha_diversity_ITS.csv',header=T,row.names = 1)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
GROUP <- select(GROUP,c(2,3))
data <- data[rownames(GROUP), ]
data <- cbind( GROUP,data)
data<-data[-c(39,44,63),]
colnames(data)[4] <- c("Richness")

data_melt <- melt(data, id.vars = c("Stage","Group"), variable.name = "index",value.name = "y")  
colnames(data_melt)[c(3,4)] <- c("index","y")
data_melt$Group[data_melt$Group=="Silt-clay particles"]<-"Silt_clay particles"
data_melt$Group <- factor(data_melt$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt_clay particles"),labels= c("Total","LMA","SMA","MIA","SCP"))
data_melt_2<-data_melt[which(data_melt$index=="Shannon"),]



data_melt_1$class<-"Bacteria"
data_melt_2$class<-"Fungi"
data_mel<-rbind(data_melt_1,data_melt_2)
data_mel$Stage <- factor(data_mel$Stage,levels = c("S1","S2","S3","S5"),labels= c("S1","S2","S3","S4"))

data_lsmp_1<-data_mel[which(data_mel$Group=="LMA"|data_mel$Group=="SMA"|data_mel$Group=="MIA"|data_mel$Group=="SCP"),]
data_lsmp<- data_lsmp_1
data_lsmp$Stage <- data_lsmp$Group
data_lsmp$Group <- "Aggregates"
data_mel<-rbind(data_lsmp,data_lsmp_1)
data_mel$Group <- factor(data_mel$Group,levels = c("Aggregates","LMA","SMA","MIA","SCP"))
blank_data <- data.frame( class= c("Bacteria","Bacteria","Fungi","Fungi"),
                          x=0,y = c(5, 8,2.1, 5))
p1<- ggboxplot(data_mel, x="Stage", y="y",  fill  =   "Stage",outlier.shape=NA,size = 0.1)+  
  facet_grid(class~Group, scales="free")+    
  theme_bw(base_size = 8)+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position = "none")+
  ylab("diversity_index_value")+
  labs( title = "ITS-GROUP") +
  scale_fill_manual(values =c(  '#709AE1FF','#71D0F5FF','#FF6F00FF','#EE4C97E5', "#f4b384",'#f78280','#77bbec',"#b6df9a"))+ 
  theme(strip.background=element_rect(colour="black",
                                      fill="white"))+
  geom_blank(data = blank_data, aes(x=x, y = y)) +
  
  theme(plot.title = element_text(hjust = 0.5),axis.text =element_text(size=10))
p1

options(rgl.useNULL=TRUE)
library("export")
graph2ppt(p1,file="z5.pptx",append=T,width=7,height=3.8)







library(vegan)
library(ggplot2)
library(ggsci)
library(cowplot)
library(egg)
library(ape)

data=read.csv('otutab.csv',header=T,row.names = 1)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
data<-t(data)

bray_dist<-vegdist(data,method = "bray")
data.pcoa<-pcoa(bray_dist,correction = "cailliez")
xlab <- paste0("PCOA1(",round(data.pcoa$values$Rel_corr_eig[1]*100,2),"%)")
ylab <- paste0("PCOA2(",round(data.pcoa$values$Rel_corr_eig[2]*100,2),"%)")

data_pcs <- data.frame(data.pcoa$vectors[,c(1,2)],Group=GROUP$Group,Stage=GROUP$Stage)
data_pcs$Axis.1<- as.numeric(data_pcs$Axis.1)
data_pcs$Axis.2<- as.numeric(data_pcs$Axis.2)
data_pcs$Group<- factor(data_pcs$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
data_pcs$Stage<- as.factor(data_pcs$Stage)
write.csv(data_pcs, 'data_pcs.csv', quote = FALSE)

data_pcs <- data_pcs[which(data_pcs$Group != "Total"),]
m1 = ggplot(data_pcs,aes(Axis.1,Axis.2,color=Group))+
  geom_point(size=1.1,alpha=0.8)+ 
  theme_bw()+
  labs(x=xlab,y=ylab,color="Group",title =expression("R"^"2"=="0.309  P=0.001"))+
  theme(legend.key.size = unit(10, "pt"),legend.background = element_rect(fill = "transparent"),legend.title=element_text(size=6),legend.text=element_text(size=4.8))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="right")+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  stat_ellipse(level = 0.95)+ 
  scale_color_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#EE4C97E5','#FF6F00FF',  '#71D0F5FF', '#709AE1FF')))

m1

m2 = ggplot(data_pcs,aes(Axis.1,Axis.2,color=Stage))+
  geom_point(size=1.1,alpha=1)+ 
  theme_bw(base_size=6)+
  labs(x=xlab,y=ylab,color="Stage",title =expression("R"^"2"=="0.362  P=0.001"))+
  theme(legend.key.size = unit(15, "pt"),legend.background = element_rect(fill = "transparent"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="right")+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  stat_ellipse(level = 0.95)+ 
  scale_color_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))
m2



data=read.csv('otutab_ITS.csv',header=T,row.names = 1)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
data<-t(data)
bray_dist<-vegdist(data,method = "bray")
data.pcoa<-pcoa(bray_dist,correction = "cailliez")
xlab <- paste0("PCOA1(",round(data.pcoa$values$Rel_corr_eig[1]*100,2),"%)")
ylab <- paste0("PCOA2(",round(data.pcoa$values$Rel_corr_eig[2]*100,2),"%)")
data_pcs <- data.frame(data.pcoa$vectors[,c(1,2)],Group=GROUP$Group,Stage=GROUP$Stage)
data_pcs$Axis.1<- as.numeric(data_pcs$Axis.1)
data_pcs$Axis.2<- as.numeric(data_pcs$Axis.2)
data_pcs$Group<- factor(data_pcs$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
data_pcs$Stage<- as.factor(data_pcs$Stage)

write.csv(data_pcs, 'data_pcs.csv', quote = FALSE)

data_pcs <- data_pcs[which(data_pcs$Group != "Total"),]
m3 = ggplot(data_pcs,aes(Axis.1,Axis.2,color=Group))+
  geom_point(size=1.1,alpha=0.8)+ 
  theme_bw(base_size=6)+
  labs(x=xlab,y=ylab,color="Group",title =expression("R"^"2"=="0.2  P=0.001"))+
  theme(legend.key.size = unit(15, "pt"),legend.background = element_rect(fill = "transparent"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="right")+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  stat_ellipse(level = 0.95)+
  scale_color_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#EE4C97E5','#FF6F00FF',  '#71D0F5FF', '#709AE1FF')))

m3

m4 = ggplot(data_pcs,aes(Axis.1,Axis.2,color=Stage))+
  geom_point(size=1.1,alpha=1)+ 
  theme_bw(base_size=6)+
  labs(x=xlab,y=ylab,color="Stage",title =expression("R"^"2"=="0.313  P=0.001"))+
  theme(legend.key.size = unit(15, "pt"),legend.background = element_rect(fill = "transparent"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),legend.position="right")+
  geom_hline(yintercept=0,linetype="dashed",size=0.3)+
  geom_vline(xintercept=0,linetype="dashed",size=0.3)+
  stat_ellipse(level = 0.95)+
  scale_color_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))
m4

m1 <- egg::set_panel_size(m1, width=unit(1.15, "in"), height=unit(1.15, "in"))
m2 <- egg::set_panel_size(m2, width=unit(1.15, "in"), height=unit(1.15, "in"))
m3 <- egg::set_panel_size(m3, width=unit(1.15, "in"), height=unit(1.15, "in"))
m4 <- egg::set_panel_size(m4, width=unit(1.15, "in"), height=unit(1.15, "in"))

p2 <- plot_grid(m1,m3,m2,m4,labels=c('16S_PCOA_Group_facet', '                16S_PCOA_Group','16S_PCOA_Stage_facet',"                16S_PCOA_Stage"),ncol = 2)  











library(ggplot2)
library(reshape)
library(ggalluvial)
library(ggsci)
library(ggpubr)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

otu1 <- read.csv('D:\\Program Files (x86)\\RStudio\\example\\otutab.csv',header=T,row.names = 1)
otu1 <- otu1[which(rowSums(otu1) >= 22), ]
otu1<-data.frame(t(otu1))
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
GROUP <- select(GROUP,c(3,2))
otu <- otu1[rownames(GROUP), ]
otu1 <- cbind(GROUP,otu)
otu3<-otu1
otu3<-otu3[which(otu3$Group!="Total"),]
otu3$Group<-as.factor(otu3$Group)
otu3$Stage<-as.factor(otu3$Stage)
write.csv(otu3, 'otu_l+s+m+s.csv', quote = FALSE)
otu3 <- read.csv('otu_l+s+m+s.csv',header=T,row.names = 1,stringsAsFactors = T)



otu1<-otu3[which(otu3$Group=="Large macroaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
library(randomForest)
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
p<-plot(otu_test$Stage, plant_predict, main = '测试集',
        xlab = 'Plant age (days)', ylab = 'Predict')
p+abline(1, 1)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 16)
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:16, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, '16S.Stage.otu.select_LMA.csv', quote = FALSE)


otu1<-otu3[which(otu3$Group=="Small macroaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
p<-plot(otu_test$Stage, plant_predict, main = '测试集',
        xlab = 'Plant age (days)', ylab = 'Predict')
p+abline(1, 1)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20 )
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:16, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, '16S.Stage.otu.select_SMA.csv', quote = FALSE)


otu1<-otu3[which(otu3$Group=="Microaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:37, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, '16S.Stage.otu.select_MIA.csv', quote = FALSE)



otu1<-otu3[which(otu3$Group=="Silt-clay particles"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:125, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, '16S.Stage.otu.select_SCP.csv', quote = FALSE)


data_frame_1 <-read.csv('16S.Stage.otu.select_LMA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\16_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]
data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
}

data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)

da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))
data<- aggregate( .~ Stage+Group, data = data, mean)
da1<-data[which(data$Group=="LMA"),]
da_melt_1 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y") 
colnames(da_melt_1)[c(3,4)] <- c("Phylum","y")



data_frame_1 <-read.csv('16S.Stage.otu.select_SMA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\16_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]
data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
} 
data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))

data<- aggregate( .~ Stage+Group, data = data, mean)

da1<-data[which(data$Group=="SMA"),]
da_melt_2 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y")
colnames(da_melt_2)[c(3,4)] <- c("Phylum","y")


data_frame_1 <-read.csv('16S.Stage.otu.select_MIA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\16_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]

data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
} 
data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))
data<- aggregate( .~ Stage+Group, data = data, mean)
da1<-data[which(data$Group=="MIA"),]
da_melt_3 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y") 
colnames(da_melt_3)[c(3,4)] <- c("Phylum","y")

data_frame_1 <-read.csv('16S.Stage.otu.select_SCP.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\16_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]

data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
} 
data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)

da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))
data<- aggregate( .~ Stage+Group, data = data, mean)
da1<-data[which(data$Group=="SCP"),]
da_melt_4 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y")  
colnames(da_melt_4)[c(3,4)] <- c("Phylum","y")


da_melt<-rbind(da_melt_1,da_melt_2,da_melt_3,da_melt_4)
da_melt$Phylum<-as.character(da_melt$Phylum)

d<- da_melt %>%
  group_by(Phylum) %>%
  summarize(sum_value = sum(y)) 
name<-d[order(d$sum_value,decreasing = F)[1:13],]$Phylum 
da_melt$Phylum[da_melt$Phylum %in% name]<-"others"            
da_melt<- aggregate( .~ Stage+Group+Phylum, data = da_melt, sum)
unique(da_melt$Phylum)
da_melt$Phylum<-factor(da_melt$Phylum,levels=c("Acidobacteria","Actinobacteria","Bacteroidetes","Chloroflexi","Gemmatimonadetes","Nitrospirae",
                                               "Planctomycetes","Proteobacteria","Latescibacteria","SM2F11","Verrucomicrobia","Unassigned","others"))

p3<-ggplot(da_melt,
           aes(x=Stage,
               y=y,
               fill=Phylum)) +   
  facet_grid(.~Group, scales = "free_x", space = "free_x")+
  scale_fill_manual(values =  rev(c('#8DD3C7', '#FFFFB3', '#BEBADA',  '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5','#FFFFB3',"#E41A1C","#1E90FF","#FF8C00","#4DAF4A",
                                    "#40E0D0","#FFC0CB","#00BFFF","#FFDEAD","#90EE90",
                                    "#00FFFF","#F0A3FF",  
                                    "#4C005C","#2BCE48","#FFCC99","#984EA3",
                                    "#808080","#8F7C00",
                                    "#C20088","#9DCC00",'#FCB2AF',"#003380",'#FB8072' ,"#0075DC", "#EE82EE",
                                    "#426600","#FFA405",'#80B1D3',"#FF0010",
                                    "#740AFF","#90EE90")))+
  geom_bar(stat='identity', width=0.8) +
   labs(x='Stage', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme_bw()+ 
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth=0.5, colour="black"),axis.text =element_text(size=6) ,legend.key.size = unit(10, "pt") )

p3



otu1 <- read.csv('D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.csv',header=T,row.names = 1)
otu1<-data.frame(t(otu1))
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
GROUP <- select(GROUP,c(3,2))
otu <- otu1[rownames(GROUP), ]
otu1 <- cbind(GROUP,otu)
otu3<-otu1
otu3 <- otu3[-c(39,44,63),]
otu3<-otu3[which(otu3$Group!="Total"),]
otu3$Group<-as.factor(otu3$Group)
otu3$Stage<-as.factor(otu3$Stage)
write.csv(otu3, 'otu_l+s+m+s.csv', quote = FALSE)
otu3 <- read.csv('otu_l+s+m+s.csv',header=T,row.names = 1,stringsAsFactors = T)



otu1<-otu3[which(otu3$Group=="Large macroaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:62, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, 'ITS.Stage.otu.select_LMA.csv', quote = FALSE)



otu1<-otu3[which(otu3$Group=="Small macroaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:42, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, 'ITS.Stage.otu.select_SMA.csv', quote = FALSE)



otu1<-otu3[which(otu3$Group=="Microaggregates"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:28, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, 'ITS.Stage.otu.select_MIA.csv', quote = FALSE)


otu1<-otu3[which(otu3$Group=="Silt-clay particles"),]
otu1<-otu1[,-2]
set.seed(123)
train <- sample(nrow(otu1), nrow(otu1)*0.7)
otu_train <- otu1[train, ]
otu_test <- otu1[-train, ]
set.seed(123)
otu_train.forest <-randomForest(Stage~.,data=otu_train,ntree=500,importance=TRUE,proximity=TRUE)    
otu_train.forest 
plant_predict <- predict(otu_train.forest, otu_train)
plant_predict <- predict(otu_train.forest, otu_test)
plant_predict1 <- predict(otu_train.forest, otu_train)
plant_predict2 <- predict(otu_train.forest, otu_test)
data1<-data.frame(otu_train$Stage,plant_predict1, check.names = FALSE)
data2<-data.frame(otu_test$Stage,plant_predict2, check.names = FALSE)
importance_otu <- data.frame(importance(otu_train.forest), check.names = FALSE)
head(importance_otu)
varImpPlot(otu_train.forest, n.var = min(30, nrow(otu_train.forest$importance)),
           main = 'Top 30 - variable importance')
set.seed(123)
otu_train.cv <- replicate(5, rfcv(otu1[-ncol(otu1)], otu1$Stage, cv.fold = 10, step = 1.5), simplify = FALSE)
otu_train.cv
otu_train.cv$error.cv
otu_train.cv <- data.frame(sapply(otu_train.cv, '[[', 'error.cv'))
otu_train.cv$otus <- rownames(otu_train.cv)
otu_train.cv <- reshape2::melt(otu_train.cv, id = 'otus')
otu_train.cv$otus <- as.numeric(as.character(otu_train.cv$otus))
otu_train.cv.mean <- aggregate(otu_train.cv$value, by = list(otu_train.cv$otus), FUN = mean)
head(otu_train.cv.mean, 20)
write.csv(otu_train.cv.mean, '16S.Stage.otu_train.cv.mean.csv', quote = FALSE)
importance_otu <- importance_otu[order(importance_otu$MeanDecreaseAccuracy, decreasing = TRUE), ]
importance_otu.select <- importance_otu[1:94, ]
importance_otu.select
otu_id.select <- rownames(importance_otu.select)
otu.select <- otu1[ ,c(otu_id.select, 'Stage')]
otu.select$Stage<-as.factor(otu.select$Stage)
write.csv(otu.select, 'ITS.Stage.otu.select_SCP.csv', quote = FALSE)

data_frame_1 <-read.csv('ITS.Stage.otu.select_LMA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\ITS_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names = 1)

data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]

data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
}

data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)

da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data <- data[-c(39,44,63),]
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))

data<- aggregate( .~ Stage+Group, data = data, mean)
da1<-data[which(data$Group=="LMA"),]
da_melt_1 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y")  
colnames(da_melt_1)[c(3,4)] <- c("Phylum","y")


data_frame_1 <-read.csv('ITS.Stage.otu.select_SMA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\ITS_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names = 1)

data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]
data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
}

data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data <- data[-c(39,44,63),]

data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))
data<- aggregate( .~ Stage+Group, data = data, mean)
da1<-data[which(data$Group=="SMA"),]
da_melt_2 <- melt(da1, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y")  
colnames(da_melt_2)[c(3,4)] <- c("Phylum","y")
data_frame_1 <-read.csv('ITS.Stage.otu.select_MIA.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\ITS_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]
data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
} 

data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data <- data[-c(39,44,63),]
data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))
data<- aggregate( .~ Stage+Group, data = data, mean)
da2<-data[which(data$Group=="MIA"),]
da_melt_3 <- melt(da2, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y") 
colnames(da_melt_3)[c(3,4)] <- c("Phylum","y")
data_frame_1 <-read.csv('ITS.Stage.otu.select_SCP.csv',header=T,row.names = 1)
data_frame_2<-read.csv('D:\\Program Files (x86)\\RStudio\\example\\column graph\\important species\\ITS_stage\\otutax_norm.csv',header=T,row.names = 1)
data_frame_3<-read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names = 1)
data_frame_2<-cbind(data_frame_2[,c(1:7)],data_frame_3)
data_frame<-data_frame_2[colnames(data_frame_1),]
data_frame<-data_frame[,-c(1,3:7)]
data_frame<- aggregate( .~ Phylum, data = data_frame, sum)
rownames(data_frame)<-data_frame[,1]
data_frame<-data_frame[,-1]
m=1
for (m in 1:80){
  data_frame[,m]<-as.data.frame((data_frame[,m]/colSums(data_frame)[m])*100)
} 

data_frame <-t(data_frame)
GROUP=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
da <- cbind(GROUP[,c(3,2)], data_frame)
data<-da
data <- data[-c(39,44,63),]

data<-data[which(data$Group!="Total"),]
data$Group <- factor(data$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"),labels = c("LMA","SMA","MIA","SCP"))
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S5"),labels =c("S1","S2","S3","S4"))


data<- aggregate( .~ Stage+Group, data = data, mean)

da3<-data[which(data$Group=="SCP"),]
da_melt_4 <- melt(da3, id.vars = c("Stage","Group"), variable.name = "Phylum",value.name = "y")  
colnames(da_melt_4)[c(3,4)] <- c("Phylum","y")

da_melt<-rbind(da_melt_1,da_melt_2,da_melt_3,da_melt_4)
da_melt$Phylum<-as.character(da_melt$Phylum)
d<- da_melt %>%
  group_by(Phylum) %>%
  summarize(sum_value = sum(y)) 
name<-d[order(d$sum_value,decreasing = F)[1:4],]$Phylum 
da_melt$Phylum[da_melt$Phylum %in% name]<-"others"            
da_melt<- aggregate( .~ Stage+Group+Phylum, data = da_melt, sum)
unique(da_melt$Phylum)
da_melt$Phylum<-factor(da_melt$Phylum,levels=c("Ascomycota","Basidiomycota","Mortierellomycota","Anthophyta","Chlorophyta","Mucoromycota",
                                               "Ciliophora","Ochrophyta", "Chytridiomycota","Rozellomycota" ,"Unassigned","others"))

A <- da_melt[which(da_melt$Group=="SCP" ),] %>% group_by(Stage,Phylum) %>%summarise(w=mean(y), sd = sd(y))
p4<-ggplot(da_melt,
           aes(x=Stage,
               y=y,
               fill=Phylum)) +   
  facet_grid(.~Group, scales = "free_x", space = "free_x")+
  scale_fill_manual(values = c("#A6CEE3" ,"#1F78B4", "#B2DF8A", "#33A02C" ,"#FB9A99", "#E31A1C" ,"#FDBF6F", "#FF7F00", "#CAB2D6" ,"#6A3D9A","#B15928","#DC7094","#B15928"))+
  geom_bar(stat='identity', width=0.8) +
  labs(x='Stage', y='Relative Abundance (%)')+
  scale_y_continuous(expand=c(0, 0))+
  theme_bw()+ 
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth=0.5, colour="black"),axis.text =element_text(size=6) ,legend.key.size = unit(10, "pt") )

p4

