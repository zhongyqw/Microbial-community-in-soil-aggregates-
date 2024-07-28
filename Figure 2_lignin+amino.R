library(linkET)
library(tidyverse)
library(ggplot2)
library(vegan)
library(openxlsx)
library(ggsci)
library(reshape2)
library(dplyr)
library(ggprism)
library(ggbeeswarm)
library(rstatix)

env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_Z.xlsx", sheet = 6,rowNames = T,)
env<-env[,-c(1,10,30)]
names(env)[c(10,27)]<-c("Amino_sugar","Lignin_phenol")
names(env)[c(19:26)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")
env$Group[which(env$Group =='Large macroaggregates')] <- 'LMA'
env$Group[which(env$Group =='Small macroaggregates')] <- 'SMA'
env$Group[which(env$Group =='Microaggregates')] <- 'MIA'
env$Group[which(env$Group =='Silt-clay particles')] <- 'SCP'
env$Stage[which(env$Stage =='S5')] <- 'S4'
env$Stage <- factor(env$Stage,levels = c("S1","S2","S3","S4"))
env$Group <- factor(env$Group,levels = c("Total","LMA","SMA","MIA","SCP"))

amino<- env %>%select(Group,Stage,SOC,Vanillin,Acetovanillone,Acetosyringone,Vanillic_acid,Syringaldehyde,Syringic_acid,Hydroxycinnamic_acid,Ferulic_acid,Lignin_phenol) %>%
  mutate(V = Vanillin+Acetovanillone+Vanillic_acid,S=Acetosyringone+Syringaldehyde+Syringic_acid,C=Hydroxycinnamic_acid+Ferulic_acid) %>%
  mutate(SV = S/V ,CV = C/V ,AdAlv=Vanillic_acid/Vanillin,AdAls=Syringic_acid/Syringaldehyde)%>%
  mutate(Lig_C = (Lignin_phenol/(1000*SOC))*100 ,V_C = (V/(1000*SOC))*100  ,S_C = (S/(1000*SOC))*100,C_C = (C/(1000*SOC))*100)%>%
  select(-SOC)

amino_me <- melt(amino, id.vars = c("Stage","Group"), variable.name = "A",value.name = "y")
amino_melt<- amino_me %>% group_by(Stage,Group,A) %>%summarise(w=mean(y), sd = sd(y)) 
amino_m <- amino_melt
library(multcompView)

aa<- c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid","Lignin_phenol","V","S","C","SV","CV","AdAlv","AdAls","Lig_C","V_C","S_C","C_C")
gg<-  c("Total","LMA","SMA","MIA","SCP")

for(i in aa)
{
  for(j in gg)
  { #browser()
    AG<- amino_me[which(amino_me$A ==i & amino_me$Group ==j),]
    anova <- aov(y~Stage,data=AG)
    tukey <- TukeyHSD(anova)
    cld <- multcompLetters4(anova,tukey)
    cld
    letter <- as.data.frame.list(cld$Stage) %>% select(1) %>% rownames_to_column("Stage") %>% mutate(A=i,Group=j)
    amino_m <- amino_m %>%left_join(.,letter,by=c("Stage","Group","A"))
  }
}
amino_m <- amino_m %>%mutate(letters=0)

for(k in 1:nrow(amino_m))
{
  for(l in 6:ncol(amino_m)-1)
  {
    if(!is.na(amino_m[k,l])){
      amino_m$letters[k]<- amino_m[k,l]
    }
    
  }
}


a <- data.frame()
for(q in aa)
{
  
  AG<- amino_me[which(amino_me$A ==q),]
  anova <- aov(y~Group,data=AG)
  tukey <- TukeyHSD(anova)
  cld <- multcompLetters4(anova,tukey)
  cld
  cld <- as.data.frame.list(cld$Group)
  cld$Letters <- toupper(cld$Letters)
  letter <- cld %>% select(1) %>% rownames_to_column("Group") %>% mutate(A=q)
  a <- rbind(a,letter)
  
}


amino <- amino_m [,c("Stage","Group","A","w","sd","letters")]
amino$Stage <- factor(amino$Stage,levels = c("S1","S2","S3","S4"))
amino$Group <- factor(amino$Group,levels = c("Total","LMA","SMA","MIA","SCP"))
amino$A <- factor(amino$A,levels = c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid","Lignin_phenol","V","S","C","SV","CV","AdAlv","AdAls","Lig_C","V_C","S_C","C_C"))


amino_lignin <- amino[which( amino$A=="Lig_C"),]
P1<-ggplot(amino_lignin,aes(x=Group,  y=w, fill=Stage))+
  facet_wrap(~A, scales="free",ncol = 4, nrow = 6)+
  geom_errorbar(aes(x=Group,ymin=w-sd,ymax=w+sd),
                width=0.2,
                stat="identity", position=position_dodge(0.86))+
  geom_bar(stat = 'identity',
           position = 'dodge',
           width = 0.86,     
           color='black',
           linewidth=0.1,
           alpha=1)+
  geom_text(aes(x=Group,y=w+sd,label=letters),
            vjust=-0.4,
            size=2.6,
            stat="identity", position=position_dodge(0.86))+  
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_y_continuous(expand=c(0, 0),limits = c(0, 5))+
  labs(x='Group', y=expression("Lignin phenols / SOC (%)"))+
  scale_x_discrete(labels  = c("Total","LMA","SMA","MIA","SCP"))+
  theme_bw(base_size = 12)+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth = 0.5, colour="black"),axis.text =element_text(size=10))

P1


amino_AdAl <- amino[which(amino$A=="AdAlv" | amino$A=="AdAls"),]
P2<-ggplot(amino_AdAl,aes(x=Group,  y=w, fill=Stage))+
  facet_wrap(~A, scales="free",ncol = 4, nrow = 6)+
  geom_errorbar(aes(x=Group,ymin=w-sd,ymax=w+sd),
                width=0.2,
                stat="identity", position=position_dodge(0.86))+
  geom_bar(stat = 'identity',
           position = 'dodge', 
           width = 0.86,     
           color='black',
           linewidth=0.1,
           alpha=1)+
  geom_text(aes(x=Group,y=w+sd,label=letters),
            vjust=-0.4,
            size=2.6,
            stat="identity", position=position_dodge(0.86))+  
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+  
  theme(axis.text = element_text(colour = 'black'))+
  scale_y_continuous(expand=c(0, 0),limits = c(0, 0.5))+
  labs(x='Group', y=expression("Degree of lignin degradation"))+
  scale_x_discrete(labels  = c("Total","LMA","SMA","MIA","SCP"))+
  theme_classic(base_size = 12)+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth = 0.5, colour="black"),axis.text =element_text(size=10))

P2


library(linkET)
library(tidyverse)
library(ggplot2)
library(vegan)
library(openxlsx)
library(ggsci)
library(reshape2)
library(dplyr)
library(ggprism)
library(ggbeeswarm)
library(rstatix)

env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_Z.xlsx", sheet = 5,rowNames = T,) 
env<-env[,-c(1,10,30)]
names(env)[c(10,27)]<-c("Amino_sugar","Lignin_phenol")
names(env)[c(19:26)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")
env$Group[which(env$Group =='Large macroaggregates')] <- 'LMA'
env$Group[which(env$Group =='Small macroaggregates')] <- 'SMA'
env$Group[which(env$Group =='Microaggregates')] <- 'MIA'
env$Group[which(env$Group =='Silt-clay particles')] <- 'SCP'
env$Stage[which(env$Stage =='S5')] <- 'S4'
env$Stage <- factor(env$Stage,levels = c("S1","S2","S3","S4"))
env$Group <- factor(env$Group,levels = c("Total","LMA","SMA","MIA","SCP"))


NC <- env[,c(1,2,12,13,14)]
library(reshape2)

amino_me <- melt(NC, id.vars = c("Stage","Group"), variable.name = "A",value.name = "y")  
amino_melt<- amino_me %>% group_by(Stage,Group,A) %>%summarise(w=mean(y), sd = sd(y)) 
amino_m <- amino_melt
library(multcompView)
aa<- c("BNC","FNC","TotalNC")
gg<-  c("Total","LMA","SMA","MIA","SCP")

for(i in aa)
{
  for(j in gg)
  { 
    AG<- amino_me[which(amino_me$A ==i & amino_me$Group ==j),]
    anova <- aov(y~Stage,data=AG)
    tukey <- TukeyHSD(anova)
    cld <- multcompLetters4(anova,tukey)
    cld
    letter <- as.data.frame.list(cld$Stage) %>% select(1) %>% rownames_to_column("Stage") %>% mutate(A=i,Group=j)
    amino_m <- amino_m %>%left_join(.,letter,by=c("Stage","Group","A"))
  }
}
amino_m <- amino_m %>%mutate(letters=0)
for(k in 1:nrow(amino_m))
{
  for(l in 6:ncol(amino_m)-1)
  {
    if(!is.na(amino_m[k,l])){
      amino_m$letters[k]<- amino_m[k,l]
    }
    
  }
}

a <- data.frame()
for(q in aa)
{
  
  AG<- amino_me[which(amino_me$A ==q),]
  anova <- aov(y~Group,data=AG)
  tukey <- TukeyHSD(anova)
  cld <- multcompLetters4(anova,tukey)
  cld
  cld <- as.data.frame.list(cld$Group)
  cld$Letters <- toupper(cld$Letters)
  
  letter <- cld %>% select(1) %>% rownames_to_column("Group") %>% mutate(A=q)
  a <- rbind(a,letter)
  
}

amino <- amino_m [,c("Stage","Group","A","w","sd","letters")]
amino$Stage <- factor(amino$Stage,levels = c("S1","S2","S3","S4"))
amino$Group <- factor(amino$Group,levels = c("Total","LMA","SMA","MIA","SCP"))
amino$A <- factor(amino$A,levels = c("BNC","FNC","TotalNC"))


P3<-ggplot(amino,aes(x=Group,  y=w, fill=Stage))+
  facet_wrap(~A, scales="free",ncol = 4, nrow = 1)+
  geom_errorbar(aes(x=Group,ymin=w-sd,ymax=w+sd),
                width=0.2,
                stat="identity", position=position_dodge(0.8))+
  geom_bar(stat = 'identity',aes(color=Stage),
           position = 'dodge',
           width = 0.8,      
           color='black',
           linewidth=0.1)+
  geom_text(aes(x=Group,y=w+sd,label=letters),
            vjust=-0.4,
            size=3,
            stat="identity", position=position_dodge(0.8))+   
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+  
  theme(axis.text = element_text(colour = 'black'))+
  labs(x='Group', y='Total Content')+
  scale_x_discrete(labels  = c("Total","LMA","SMA","MIA","SCP"))+
  theme_bw(base_size = 12)+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth = 0.5, colour="black"))

P3


NC <- env[,c(1,2,15,16,17)]
library(reshape2)
NC<- NC[which(NC$`TotalNC/SOC` <1 ),]
amino_me <- melt(NC, id.vars = c("Stage","Group"), variable.name = "A",value.name = "y") 
amino_me$y <- amino_me$y*100

amino_melt<- amino_me %>% group_by(Stage,Group,A) %>%summarise(w=mean(y), sd = sd(y)) 
amino_m <- amino_melt
library(multcompView)

aa<- c("BNC/SOC","FNC/SOC","TotalNC/SOC")
gg<-  c("Total","LMA","SMA","MIA","SCP")

for(i in aa)
{
  for(j in gg)
  { 
    AG<- amino_me[which(amino_me$A ==i & amino_me$Group ==j),]
    anova <- aov(y~Stage,data=AG)
    tukey <- TukeyHSD(anova)
    cld <- multcompLetters4(anova,tukey)
    cld
    letter <- as.data.frame.list(cld$Stage) %>% select(1) %>% rownames_to_column("Stage") %>% mutate(A=i,Group=j)
    amino_m <- amino_m %>%left_join(.,letter,by=c("Stage","Group","A"))
  }
}
amino_m <- amino_m %>%mutate(letters=0)
for(k in 1:nrow(amino_m))
{
  for(l in 6:ncol(amino_m)-1)
  {
    if(!is.na(amino_m[k,l])){
      amino_m$letters[k]<- amino_m[k,l]
    }
    
  }
}

a <- data.frame()
for(q in aa)
{
  
  AG<- amino_me[which(amino_me$A ==q),]
  anova <- aov(y~Group,data=AG)
  tukey <- TukeyHSD(anova)
  cld <- multcompLetters4(anova,tukey)
  cld
  cld <- as.data.frame.list(cld$Group)
  cld$Letters <- toupper(cld$Letters)
  
  letter <- cld %>% select(1) %>% rownames_to_column("Group") %>% mutate(A=q)
  a <- rbind(a,letter)
  
}

amino <- amino_m [,c("Stage","Group","A","w","sd","letters")]
amino <- merge(amino,a,by=c("Group","A"))
amino <- amino %>% group_by(Group,A) %>%summarise(LW=mean(w), LSD = max(sd)) %>% merge(amino,by=c("Group","A"))
amino$Stage <- factor(amino$Stage,levels = c("S1","S2","S3","S4"))
amino$Group <- factor(amino$Group,levels = c("Total","LMA","SMA","MIA","SCP"))
amino$A <- factor(amino$A,levels = c("BNC/SOC","FNC/SOC","TotalNC/SOC"))


P4<-ggplot(amino,aes(x=Group,  y=w, fill=Stage))+
  facet_wrap(~A, scales="free",ncol = 4, nrow = 1)+
  geom_errorbar(aes(x=Group,ymin=w-sd,ymax=w+sd),
                width=0.2,
                stat="identity", position=position_dodge(0.8))+
  geom_bar(stat = 'identity',aes(color=Stage),
           position = 'dodge',
           width = 0.8,  
           color='black',
           linewidth=0.1)+
  geom_text(aes(x=Group,y=w+sd,label=letters),
            vjust=-0.4,
            size=3,
            stat="identity", position=position_dodge(0.8))+  
  scale_fill_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))+ 
  theme(axis.text = element_text(colour = 'black'))+
  labs(x='Group', y='Total Content')+
  scale_x_discrete(labels  = c("Total","LMA","SMA","MIA","SCP"))+
  theme_classic(base_size = 12)+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth = 0.5, colour="black"))
P4


p1 <- egg::set_panel_size(P1, width=unit(2, "in"), height=unit(2, "in"))
p2 <- egg::set_panel_size(P2, width=unit(2, "in"), height=unit(2, "in"))
p3 <- egg::set_panel_size(P3, width=unit(2, "in"), height=unit(2, "in"))
p4 <- egg::set_panel_size(P4, width=unit(2, "in"), height=unit(2, "in"))
library(cowplot)
p <- plot_grid(p1,p2,p3,p4,p5,ncol = 1,nrow = 5)  
ggsave("ass_4.pdf", egg::set_panel_size(p, width=unit(10, "in"), height=unit(15, "in")), 
       width=12, height = 16, units = 'in')
