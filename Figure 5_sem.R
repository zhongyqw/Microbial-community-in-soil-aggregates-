library(reshape2)
library(openxlsx)
library(magrittr)
library(tidyverse)
library(vegan)
library(dplyr)
ENV<- read.xlsx("ENV_Z_SEM.xlsx", sheet =4,rowNames = T,) 
names(ENV)[c(12:20)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid","Lignin_phenol")
ENV<- ENV %>%mutate(V = Vanillin+Acetovanillone+Vanillic_acid,
                    S=Acetosyringone+Syringaldehyde+Syringic_acid,
                    C=Hydroxycinnamic_acid+Ferulic_acid, 
                    AdAlv=Vanillic_acid/Vanillin) 
ENV <- ENV[,-c(12:19)]
names(ENV)[12] <- "lig"
otutab <- read.csv("data_shanoon_16S.csv",header = T,row.names = 1) 
data<- merge(ENV,otutab[,c(1,2)],by=0) %>%column_to_rownames(., var = "Row.names") 
data<- dplyr::select(data,c(1:18))
otutab <- read.csv("data_shanoon_ITS.csv",header = T,row.names = 1)
data<- merge(data,otutab[,c(1,2)],by=0) %>%column_to_rownames(., var = "Row.names")

data <- data %>% filter(MIA ==100|SCP==100 )
dt1_psem <- psem(
  lm(lig ~  Stage, data = data),
  lm(Shannon.b ~Stage+lig, data = data),
  lm(Shannon.f ~ Stage+lig, data = data),
  lm(BNC ~Shannon.b   +Shannon.f   , data = data),
  lm(FNC ~Shannon.b   +Shannon.f , data = data),
  lm(SOC ~  Stage +BNC +FNC , data = data),
  Shannon.b %~~% Shannon.f,
  BNC %~~% FNC
)

summary(dt1_psem,  standardized = TRUE)
plot(dt1_psem,return = FALSE,alpha = 0.05,show = "std")


ENV<- read.xlsx("ENV_Z_SEM.xlsx", sheet =4,rowNames = T,) 
names(ENV)[c(12:20)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid","Lignin_phenol")
ENV<- ENV %>%mutate(V = Vanillin+Acetovanillone+Vanillic_acid,
                    S=Acetosyringone+Syringaldehyde+Syringic_acid,
                    C=Hydroxycinnamic_acid+Ferulic_acid, 
                    AdAlv=Vanillic_acid/Vanillin) 
ENV <- ENV[,-c(12:19)]
names(ENV)[12] <- "lig"
otutab <- read.csv("data_shanoon_16S.csv",header = T,row.names = 1) 
data<- merge(ENV,otutab[,c(1,2)],by=0) %>%column_to_rownames(., var = "Row.names") 
data<- dplyr::select(data,c(1:18))
otutab <- read.csv("data_shanoon_ITS.csv",header = T,row.names = 1)
data<- merge(data,otutab[,c(1,2)],by=0) %>%column_to_rownames(., var = "Row.names")
data <- data %>% filter(LMA ==100|SMA==100 )
dt1_psem <- psem(
  lm(lig ~  Stage, data = data),
  lm(Shannon.b ~Stage+lig, data = data),
  lm(Shannon.f ~ Stage+lig, data = data),
  
  lm(BNC ~Shannon.b    , data = data),
  lm(FNC ~Stage +Shannon.b   +Shannon.f , data = data),
  lm(SOC ~  Stage +BNC +FNC , data = data),
  Shannon.b %~~% Shannon.f,
  BNC %~~% FNC
)

summary(dt1_psem,  standardized = TRUE)
plot(dt1_psem,return = FALSE,alpha = 0.05,show = "std")


library(reshape2)
library(openxlsx)
library(magrittr)
library(tidyverse)
library(vegan)
library(dplyr)
ENV<- read.xlsx("zhijie.xlsx", sheet =1) 
ENV <- ENV[which(ENV$G!="BNC"),]
ENV$factor <- factor(ENV$factor,levels = c("succession","lignin","bacterial cummunity","fungal cummunity","BNC","FNC"),labels = c("succession","lignin phenols","bacterial","fungal","BNC","FNC"))
P1<-ggplot(ENV,aes(x=factor,  y=y, fill=type))+
  facet_wrap(~G, scales="free",ncol = 4, nrow = 1)+
  geom_bar(stat = 'identity',
           position = 'stack', 
           width = 0.8,    
           color='black',
           linewidth=0.1,
           alpha=1)+
scale_fill_manual(values =c(  "#FD9062",'#90A0D8'))+ 
  theme(axis.text = element_text(colour = 'black'))+
  theme_classic(base_size = 12)+
  theme(panel.grid=element_blank(),axis.line=element_line(linewidth = 0.5, colour="black"))+
  coord_flip()
P1
