library(ggplot2)
library(reshape2)
library(ggalluvial)
library(ggsci)
library(ggpubr)
library(tidyverse)
library("RColorBrewer")
library(readxl)

data_frame <-read_excel('distribution.xlsx')
data_frame$Group <- factor(data_frame$Group,levels = c("Large macroaggregates","Small macroaggregates","Microaggregates","Silt_clay particles"),labels=c("LMA","SMA","MIA","SCP"))
data_frame$Stage <- factor(data_frame$Stage,levels = c("S1","S2","S3","S4"))
data<- data_frame %>% group_by(Stage,Group) %>%summarise(w=mean(distribution), se=mean(se)) 
data=read.csv('mean.csv',header=T, row.names=1)
data$Stage <- factor(data$Stage,levels = c("S1","S2","S3","S4"))
data$Group <- factor(data$Group,levels = c("LMA","SMA","MIA","SCP"))


p1<-ggplot(data,
           aes(x=Group,
               y=w,
               fill=Stage)) +   
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+  
  geom_errorbar(aes(ymin=w-se,ymax=w+se),
                width=0.3,  
                linewidth=0.4,    
                stat="identity", 
                position=position_dodge(0.86))+
  geom_bar(stat = 'identity',
           position = 'dodge', 
           width = 0.86,      
           color='black',
           linewidth=0.2,
           alpha=1)+
  geom_text(aes(y=w+se,label = c("a","c","d","a","d","b","b","b","b","a","c","d","c","d","a","c")), 
            stat='identity',
            vjust=-0.4,
            size=3,
            position=position_dodge(0.86))+
  labs(x='Stage', y='Aggregate size distribution (%)')+
  scale_y_continuous(expand=c(0, 0),limits =c(0, 60))+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth=0.5, colour="black"),axis.text =element_text(size=10))
p1
ggsave("aggregate_distribution.pdf", width = 3.5, height = 2.5) 




library(ggplot2)
library(reshape2)
library(ggalluvial)
library(ggsci)
library(ggpubr)
library(tidyverse)
library("RColorBrewer")
library(multcompView)
library(readxl)
library(openxlsx)
env<- read.xlsx("SOC.xlsx", sheet = 1,rowNames = T,) 
env <- env[,c(2,3,4)]
env$Group[which(env$Group =='Large macroaggregates')] <- 'LMA'
env$Group[which(env$Group =='Small macroaggregates')] <- 'SMA'
env$Group[which(env$Group =='Microaggregates')] <- 'MIA'
env$Group[which(env$Group =='Silt-clay particles')] <- 'SCP'
env$Stage[which(env$Stage =='S5')] <- 'S4'

env$Stage <- factor(env$Stage,levels = c("S1","S2","S3","S4"))
env$Group <- factor(env$Group,levels = c("Total","LMA","SMA","MIA","SCP"))

amino_melt<- env %>% group_by(Stage,Group) %>%summarise(w=mean(SOC), sd = sd(SOC))
amino_melt$Stage <- factor(amino_melt$Stage,levels = c("S1","S2","S3","S4"))
amino_melt$Group <- factor(amino_melt$Group,levels = c("Total","LMA","SMA","MIA","SCP"))



p2<-ggplot(amino_melt,
           aes(x=Group,
               y=w,
               fill=Stage)) + 
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+ 
   geom_errorbar(aes(ymin=w-sd,ymax=w+sd),
                width=0.3,  
                linewidth=0.4,   
                stat="identity", 
                position=position_dodge(0.86))+
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.86,    
           linewidth=0.2,
           alpha=1)+
  geom_text(aes(y=w+sd,label = c("c","a","b","d","b","c","c","b","c","c","b","a","a","b","d","a","b","a","a","a")),
            stat='identity',
            vjust=-0.4,
            size=3,
            position=position_dodge(0.86))+
  labs(x='Stage', y=expression("SOC content"* " (g kg"^"-1" *")"))+
  
  scale_y_continuous(expand=c(0, 0),limits =c(0, 25))+
  theme_bw()+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth=0.5, colour="black"),axis.text =element_text(size=10)) 
p2
ggsave("aggregate_SOC.pdf", width = 3.58, height = 2.5) 




env<- read.xlsx("SOC.xlsx", sheet = 1,rowNames = T,) 
env <- env[,c(2,3,4)]
env$Group[which(env$Group =='Large macroaggregates')] <- 'LMA'
env$Group[which(env$Group =='Small macroaggregates')] <- 'SMA'
env$Group[which(env$Group =='Microaggregates')] <- 'MIA'
env$Group[which(env$Group =='Silt-clay particles')] <- 'SCP'
env$Stage[which(env$Stage =='S5')] <- 'S4'

env$Stage <- factor(env$Stage,levels = c("S1","S2","S3","S4"))
env$Group <- factor(env$Group,levels = c("Total","LMA","SMA","MIA","SCP"))
env<- env[which(env$Group !='Total'),]

env$r <- rep(c(1,2,3),times=16)
env_1<- env[which(env$Stage !='S1'),]
env_2<- env[which(env$Stage =='S1'),]

env_3<-merge(env_1,env_2,by=c("Group","r")) 
env_3$SOC <- (env_3$SOC.x-env_3$SOC.y)*100/env_3$SOC.y 
names(env_3)[3] <- "Stage"
ENV_4<- merge(env,env_3[,c(1,2,3,7)],by=c("Stage","Group","r"),all = F) %>% select(.,c(1,2,3,5))
names(ENV_4)[4] <- "y"
amino_melt<- ENV_4 %>% group_by(Stage,Group) %>%summarise(w=mean(y), sd = sd(y)) 
amino_S1<- env %>% group_by(Stage,Group) %>%summarise(w=mean(SOC), sd = sd(SOC)) %>% .[which(.$Stage=="S1"),]
amino_melt<- rbind(amino_S1,amino_melt) 
amino_melt[which(amino_melt$Stage=="S1"),]$w <- 0

p3<-ggplot(amino_melt,
            aes(x=Stage,
                y=w,
                group = Group,
                fill=Group,
                colour = Group)) + 
  scale_color_manual(values =rev(c( '#EE4C97E5','#FF6F00FF',  '#71D0F5FF', '#709AE1FF')))+ 
  geom_errorbar(aes(ymin=w-sd,ymax=w+sd),
                width=0.1, 
                linewidth=0.2,    
                stat="identity" )+
  geom_line(mapping=aes(x=Stage,
                        y=w,
                        colour = Group), 
           stat = 'identity',
            linetype = 1)+
  geom_point(aes(x=Stage,
                 y=w,
                 colour = Group))+
  geom_text(aes(y=w+sd,label = c(" "," "," "," ","D","B","A","C","C","B","A","D","C","B","A","B")),
            stat='identity',
            vjust=-0.4,
            size=2.6)+
  labs(x='Stage', y=expression("Changes in SOC (%)"))+
  geom_hline(yintercept=0,linetype="dashed",linewidth=0.3)+
  scale_y_continuous(limits =c(-50, 150))+ 
  theme_bw()+ 
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth=0.5, colour="black"),axis.text =element_text(size=10))
p3

ggsave("SOC_CHANGE.pdf", width = 3.5, height = 2.5) 

