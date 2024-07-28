library(devtools)
library(ggvegan)
library(vegan)
library(openxlsx)
library(ggrepel) 
library(ggplot2)


otu <- read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names=1)
otu <- data.frame(t(otu))
otu.hell <- decostand(otu, "hellinger")
env.data<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\RDA\\environment.xlsx", sheet = 1,rowNames = T,)
env.data<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_WU.xlsx", sheet = 1,rowNames = T,)
env.data<- env.data[,-c(7,10,11)]
env.data<-env.data[,-1]
names(env.data)[c(7,16)]<-c("Amino_sugar","Lignin_phenol")

env.data.log <- log1p(env.data)
env <- na.omit(env.data.log)

otu.hell <- decostand(otu, "hellinger")
idx = rownames(otu.hell) %in% rownames(env)
otu.hell = otu.hell[idx,]
env = env[rownames(otu.hell),]
sel <- decorana(otu.hell)
sel

db_rda <- capscale(otu.hell~., env, distance = 'bray', add = TRUE)
db_rda

db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.scaling1

r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared 
db_rda_noadj
db_rda_adj <- r2$adj.r.squared     
db_rda_adj

db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_exp_adj 
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi
db_rda_eig_adj

db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)
db_rda_test_axis
db_rda_test_term <- anova.cca(db_rda, by = 'term', permutations = 999)
db_rda_test_term
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')
db_rda_test_axis$`Pr(>F)`


vif.cca(db_rda)
db_rda_forward_pr <- ordiR2step(capscale(otu.hell~1, env, distance = 'bray', add = TRUE), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
db_rda_forward_pr <- capscale( otu.hell ~ Lignin_phenol  + Amino_sugar +TN +SOC+pH  , env, distance = 'bray', add = TRUE)
db_rda_forward_pr
vif.cca(db_rda_forward_pr)


summary(db_rda_forward_pr, scaling = 1)
RsquareAdj(db_rda)$adj.r.squared
RsquareAdj(db_rda_forward_pr)$adj.r.squared

db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)
db_rda_forward_pr_test

db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)
db_rda_forward_pr_test_axis

db_rda_forward_pr_test_term <- anova.cca(db_rda_forward_pr, by = 'term', permutations = 999)
db_rda_forward_pr_test_term
vif.cca(db_rda_forward_pr)


db_rda_forward_pr_site.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'wa')
db_rda_forward_pr_sp.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'sp')
db_rda_forward_pr_env.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'bp')

smry <- summary(db_rda_forward_pr, scaling = 1)
df1  <- data.frame(smry$sites[,1:2])       
df2  <- data.frame(smry$species[,1:2])     
df3  <- data.frame(smry$biplot[,1:2])

group=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
group$Group<- factor(group$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
df11<-merge(df1,group,by = "row.names", all = TRUE)
df11<-na.omit(df11)
exp_adj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
rda1_exp <- paste('RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('RDA2:', round(exp_adj[2]*100, 2), '%')

rda.plot <- ggplot(df11, aes(x=CAP1, y=CAP2,colour=Stage)) +
  geom_point(aes(shape=Group),size=3,stroke=0.5)+
  geom_hline(yintercept=0, linetype="dotted",lwd=0.3) +
  geom_vline(xintercept=0, linetype="dotted",lwd=0.3)+
  scale_shape_manual(values = c(12,16,17,3,18))

rda.plot
rda.biplot <- rda.plot +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2),color="BLACK",linewidth = 0.4, arrow=arrow(length=unit(0.03,"npc")))+
  geom_text_repel(data=df3, aes(x=CAP1,y=CAP2,label=rownames(df3)),color="BLACK", size=2)+
  labs(x = rda1_exp, y = rda2_exp) +
  theme_bw(base_size=10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_color_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))

P1<-rda.biplot
P1


library(tidyverse)
library("randomForest")
library("rfUtilities")
library("rfPermute")
library(ggsci)
library(patchwork)
library(cowplot)
library(openxlsx)

otu <- read.csv('data_pcs.csv', row.names = 1,header = TRUE, fill = TRUE,)
env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\RDA\\environment.xlsx", sheet = 1,rowNames = T,)
env<-env[,c(2,3,4,8,17)]
names(env)[c(4,5)]<-c("Amino_sugar","Lignin_phenol")
group <-read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
group <-group[,c(3,2)]
group$Group<- factor(group$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
group$Stage<- factor(group$Stage,levels = c("S1","S2","S3","S5"),labels = c("S1","S2","S3","S4"))
otu <- otu[rownames(env), ]
group <- group[rownames(env), ]
otu_P <- cbind(otu[,1:2], env)
RFdata2 <- otu_P[,-2]
set.seed(123)
richness_rf <- randomForest(Axis.1 ~ ., data= RFdata2,
                            importance=TRUE,proximity=TRUE)
richness_rf
set.seed(123)
richness_perm <- rf.significance(richness_rf, RFdata2[,-1], nperm=99, ntree=501)
richness_perm
set.seed(123)
richness_rfP<- rfPermute(Axis.1  ~ ., data = RFdata2, ntree = 500, 
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)
p2 <- richness_dat %>% 
  as_tibble(rownames = "names") %>% 
  data.frame() %>% 
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**", 
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>% 
  arrange(X.IncMSE) %>% 
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>% 
  ggplot(aes(x = names, y = X.IncMSE))+
  theme_bw(base_size = 10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.position = "none")+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(title =expression("R"^"2"=="0.346  P=0.01 "))+
  theme(plot.title = element_text(size = 8))+
  labs(x = NULL, y = "PC1")+
  scale_fill_manual(values =  rev(c( "#B0A875",'#AF7EC0')))+
  coord_flip()
p2

RFdata2 <- otu_P[,-1]
set.seed(123)
richness_rf <- randomForest(Axis.2 ~ ., data= RFdata2,
                            importance=TRUE,proximity=TRUE)
richness_rf
set.seed(123)
richness_perm <- rf.significance(richness_rf, RFdata2[,-1], nperm=99, ntree=501)
richness_perm
set.seed(123)
richness_rfP<- rfPermute(Axis.2  ~ ., data = RFdata2, ntree = 500, 
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)
p3 <- richness_dat %>% 
  as_tibble(rownames = "names") %>% 
  data.frame() %>% 
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**", 
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>% 
  arrange(X.IncMSE) %>% 
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>% 
  ggplot(aes(x = names, y = X.IncMSE))+
  theme_bw(base_size = 10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.position = "none")+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(title =expression("R"^"2"=="0.725  P=0.01 "))+
  theme(plot.title = element_text(size = 8))+
  labs(x = NULL, y = "PC2")+
  scale_fill_manual(values =  c( "#B0A875",'#AF7EC0'))+
  coord_flip()
p3


library(linkET)
library(ggtext)
library(dplyr)
library(ggplot2)
library(openxlsx)
env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_WU.xlsx", sheet =1,rowNames = T,) 
env<-env[,-c(1,5,7)]
names(env)[c(7,18)]<-c("Amino_sugar","Lignin_phenol")
names(env)[c(10:17)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")

amino<- env %>%
  mutate(V = Vanillin+Acetovanillone+Vanillic_acid,S=Acetosyringone+Syringaldehyde+Syringic_acid,C=Hydroxycinnamic_acid+Ferulic_acid) %>%
  mutate(SV = S/V ,CV = C/V ,AdAlv=Vanillic_acid/Vanillin,AdAls=Syringic_acid/Syringaldehyde)%>%
  mutate(Lig_C = (Lignin_phenol/(1000*SOC))*100 ,V_C = (V/(1000*SOC))*100  ,S_C = (S/(1000*SOC))*100,C_C = (C/(1000*SOC))*100)%>%
  select(c("pH","SOC","TN","Amino_sugar","BNC","FNC","Lignin_phenol","V","S","C","AdAlv","AdAls"))
env <- amino
names(env)[c(4:12)]<-c("Amino_sugar","BNC","FNC","Lignin_phenols","Vanillyl_phenols","Syringyl_phenols","Cinnamyl_phenols","AdAl_V","AdAl_S")

mantel_funciton1 <- function(class){
  library(linkET)
  library(tidyverse)
  library(ggplot2)
  library(vegan)
  library(openxlsx)
  otu <- read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab.CSV",header=T,row.names=1)
  env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_WU.xlsx", sheet =1,rowNames = T,) 
  env<-env[,-c(1,5,7)]
  names(env)[c(7,18)]<-c("Amino_sugar","Lignin_phenol")
  names(env)[c(10:17)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")
  amino<- env %>%
    mutate(V = Vanillin+Acetovanillone+Vanillic_acid,S=Acetosyringone+Syringaldehyde+Syringic_acid,C=Hydroxycinnamic_acid+Ferulic_acid) %>%
    mutate(SV = S/V ,CV = C/V ,AdAlv=Vanillic_acid/Vanillin,AdAls=Syringic_acid/Syringaldehyde)%>%
    mutate(Lig_C = (Lignin_phenol/(1000*SOC))*100 ,V_C = (V/(1000*SOC))*100  ,S_C = (S/(1000*SOC))*100,C_C = (C/(1000*SOC))*100)%>%
    select(c("pH","SOC","TN","Amino_sugar","BNC","FNC","Lignin_phenol","V","S","C","AdAlv","AdAls"))
  env <- amino
  names(env)[c(4:12)]<-c("Amino_sugar","BNC","FNC","Lignin_phenols","Vanillyl_phenols","Syringyl_phenols","Cinnamyl_phenols","AdAl_V","AdAl_S")
  group=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
  group<-group[,-1]
  group<-group[which(group$Group!="Total"),]
  group$Group[which(group$Group=="Large macroaggregates")]<-"LMA+SMA"
  group$Group[which(group$Group=="Small macroaggregates")]<-"LMA+SMA"
  group$Group[which(group$Group=="Microaggregates")]<-"MIA+SCP"
  group$Group[which(group$Group=="Silt-clay particles")]<-"MIA+SCP"
  group$Group<- factor(group$Group,levels = c("LMA+SMA","MIA+SCP"),labels = c("LMA+SMA","MIA+SCP"))
  group<-filter(group,Group == class)
  otu <- as.data.frame(t(otu))
  otu <- otu[rownames(group), ]
  env <- env[rownames(group), ]
  dist.abund <- vegdist(otu, method = 'bray')
  dist.env1 <- map(env,dist, method = 'euclidean')
  mantel1 <- function(dist){
    result <- mantel(dist.abund, dist, method = 'spearman', permutations = 9999, na.rm = TRUE)
    r<- c(result[["statistic"]],result[["signif"]])
    return(r) }
  mantelr <- map(dist.env1,mantel1)
  mantel_result<- t(as.data.frame(mantelr))
  scale.env <- scale(env, center = TRUE, scale = TRUE)
  dist.env <- dist(scale.env, method = 'euclidean')
  mantel <- cbind(class,row.names(mantel_result),mantel_result) %>% as.data.frame
  names(mantel) <- c('spec',"env","r",'P')
  return(mantel)
}

library(parallel)
clnum <- detectCores()
cl<- makeCluster(2)
cname1 <-  c("LMA+SMA","MIA+SCP")
result1 <- parLapply(cl,cname1, mantel_funciton1)
stopCluster(cl)
mantel_s <- rbind(result1[[1]],result1[[2]])
mantel_sn <-  as.data.frame(lapply(mantel_s[,3:4],as.numeric))
mantel_s <- cbind(mantel_s[1:2],mantel_sn)
mantel2 <- mantel_s%>% 
  filter(env!="env") %>% 
  filter(P < 0.05) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, Inf),
                  labels = c("< 0.4", "0.4 - 0.6", ">= 0.6")),
         pd = cut(P, breaks = c(-Inf, 0.01, 0.05),
                  labels = c("< 0.01", "0.01 - 0.05")))
p4<-qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel2, 
              curvature = nice_curvature()) +
  scale_size_manual(values = c(0.3, 0.6, 1.2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))+
  theme(axis.text=element_text(colour='black',size=9),legend.text=element_text(colour='black',size=7),legend.title = element_text(colour='black',size=8))+
  theme(panel.background = element_blank())+
  scale_fill_continuous(low = "#7A7EE3", high = "#DA6349")

p4


otu <- read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names=1)
otu <- data.frame(t(otu))
env.data<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\RDA\\environment.xlsx", sheet = 1,rowNames = T,)
names(env.data)[c(7,16)]<-c("Amino_sugar","Lignin_phenol")
env.data.log <- log1p(env.data)
env <- na.omit(env.data.log)

otu.hell <- decostand(otu, "hellinger")
idx = rownames(otu.hell) %in% rownames(env)
otu.hell = otu.hell[idx,]
env = env[rownames(otu.hell),]
sel <- decorana(otu.hell)
sel

db_rda <- capscale(otu.hell~., env, distance = 'bray', add = TRUE)
db_rda

db_rda.scaling1 <- summary(db_rda, scaling = 1)
db_rda.scaling1

r2 <- RsquareAdj(db_rda)
db_rda_noadj <- r2$r.squared 
db_rda_noadj
db_rda_adj <- r2$adj.r.squared     
db_rda_adj

db_rda_exp_adj <- db_rda_adj * db_rda$CCA$eig/sum(db_rda$CCA$eig)
db_rda_exp_adj 
db_rda_eig_adj <- db_rda_exp_adj * db_rda$tot.chi
db_rda_eig_adj

db_rda_test <- anova.cca(db_rda, permutations = 999)
db_rda_test
db_rda_test_axis <- anova.cca(db_rda, by = 'axis', permutations = 999)
db_rda_test_axis
db_rda_test_term <- anova.cca(db_rda, by = 'term', permutations = 999)
db_rda_test_term
db_rda_test_axis$`Pr(>F)` <- p.adjust(db_rda_test_axis$`Pr(>F)`, method = 'bonferroni')
db_rda_test_axis$`Pr(>F)`


vif.cca(db_rda)
db_rda_forward_pr <- ordiR2step(capscale(otu.hell~1, env, distance = 'bray', add = TRUE), scope = formula(db_rda), R2scope = TRUE, direction = 'forward', permutations = 999)
db_rda_forward_pr <- capscale( otu.hell ~ Lignin_phenol  + Amino_sugar +TN +SOC+pH  , env, distance = 'bray', add = TRUE)
db_rda_forward_pr
vif.cca(db_rda_forward_pr)


summary(db_rda_forward_pr, scaling = 1)
RsquareAdj(db_rda)$adj.r.squared
RsquareAdj(db_rda_forward_pr)$adj.r.squared

db_rda_forward_pr_test <- anova.cca(db_rda_forward_pr, permutations = 999)
db_rda_forward_pr_test

db_rda_forward_pr_test_axis <- anova.cca(db_rda_forward_pr, by = 'axis', permutations = 999)
db_rda_forward_pr_test_axis

db_rda_forward_pr_test_term <- anova.cca(db_rda_forward_pr, by = 'term', permutations = 999)
db_rda_forward_pr_test_term
vif.cca(db_rda_forward_pr)


db_rda_forward_pr_site.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'wa')
db_rda_forward_pr_sp.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'sp')
db_rda_forward_pr_env.scaling1 <- scores(db_rda_forward_pr, choices = 1:2, scaling = 1, display = 'bp')

smry <- summary(db_rda_forward_pr, scaling = 1)
df1  <- data.frame(smry$sites[,1:2])       
df2  <- data.frame(smry$species[,1:2])     
df3  <- data.frame(smry$biplot[,1:2])

group=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
group$Group<- factor(group$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
df11<-merge(df1,group,by = "row.names", all = TRUE)
df11<-na.omit(df11)
exp_adj <- RsquareAdj(db_rda_forward_pr)$adj.r.squared * db_rda_forward_pr$CCA$eig/sum(db_rda_forward_pr$CCA$eig)
rda1_exp <- paste('RDA1:', round(exp_adj[1]*100, 2), '%')
rda2_exp <- paste('RDA2:', round(exp_adj[2]*100, 2), '%')

rda.plot <- ggplot(df11, aes(x=CAP1, y=CAP2,colour=Stage)) +
  geom_point(aes(shape=Group),size=3,stroke=0.5)+
  geom_hline(yintercept=0, linetype="dotted",lwd=0.3) +
  geom_vline(xintercept=0, linetype="dotted",lwd=0.3)+
  scale_shape_manual(values = c(12,16,17,3,18))

rda.plot
rda.biplot <- rda.plot +
  geom_segment(data=df3, aes(x=0, xend=CAP1, y=0, yend=CAP2),color="BLACK",linewidth = 0.4, arrow=arrow(length=unit(0.03,"npc")))+
  geom_text_repel(data=df3, aes(x=CAP1,y=CAP2,label=rownames(df3)),color="BLACK", size=2)+
  labs(x = rda1_exp, y = rda2_exp) +
  theme_bw(base_size=10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+ 
  scale_color_manual(values =c(  "#f4b384",'#f78280','#77bbec',"#b6df9a"))

P5<-rda.biplot
P5




otu <- read.csv('data_pcs.csv', row.names = 1,header = TRUE, fill = TRUE,)
env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\RDA\\environment.xlsx", sheet = 1,rowNames = T,)
env<-env[,c(2,3,4,8,17)]
names(env)[c(4,5)]<-c("Amino_sugar","Lignin_phenol")
group <-read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
group <-group[,c(3,2)]
group$Group<- factor(group$Group,levels = c("Total","Large macroaggregates","Small macroaggregates","Microaggregates","Silt-clay particles"))
group$Stage<- factor(group$Stage,levels = c("S1","S2","S3","S5"),labels = c("S1","S2","S3","S4"))
otu <- otu[rownames(env), ]
group <- group[rownames(env), ]
otu_P <- cbind(otu[,1:2], env)
RFdata2 <- otu_P[,-2]
set.seed(123)
richness_rf <- randomForest(Axis.1 ~ ., data= RFdata2,
                            importance=TRUE,proximity=TRUE)
richness_rf
set.seed(123)
richness_perm <- rf.significance(richness_rf, RFdata2[,-1], nperm=99, ntree=501)
richness_perm
set.seed(123)
richness_rfP<- rfPermute(Axis.1  ~ ., data = RFdata2, ntree = 500, 
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)
p6 <- richness_dat %>% 
  as_tibble(rownames = "names") %>% 
  data.frame() %>% 
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**", 
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>% 
  arrange(X.IncMSE) %>% 
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>% 
  ggplot(aes(x = names, y = X.IncMSE))+
  theme_bw(base_size = 10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.position = "none")+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(title =expression("R"^"2"=="0.792  P=0.01 "))+
  theme(plot.title = element_text(size = 8))+
  labs(x = NULL, y = "PC1")+
  scale_fill_manual(values =  rev(c( "#B0A875",'#AF7EC0')))+
  coord_flip()
p6

RFdata2 <- otu_P[,-1]
set.seed(123)
richness_rf <- randomForest(Axis.2 ~ ., data= RFdata2,
                            importance=TRUE,proximity=TRUE)
richness_rf
set.seed(123)
richness_perm <- rf.significance(richness_rf, RFdata2[,-1], nperm=99, ntree=501)
richness_perm
set.seed(123)
richness_rfP<- rfPermute(Axis.2  ~ ., data = RFdata2, ntree = 500, 
                         na.action = na.omit, nrep = 100,num.cores = 1)
richness_dat <- importance(richness_rfP, sort.by = NULL, decreasing = TRUE)
p7 <- richness_dat %>% 
  as_tibble(rownames = "names") %>% 
  data.frame() %>% 
  mutate(label = if_else(X.IncMSE.pval < 0.001,"***",
                         if_else(X.IncMSE.pval <0.01,"**", 
                                 if_else(X.IncMSE.pval<0.05,"*","ns"))),
         X.IncMSE = as.numeric(X.IncMSE)) %>% 
  arrange(X.IncMSE) %>% 
  mutate(group = if_else(label=="ns","In_sig","Sig"),
         names = forcats::fct_inorder(names)) %>% 
  ggplot(aes(x = names, y = X.IncMSE))+
  theme_bw(base_size = 10)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"),
        legend.position = "none")+
  geom_bar(aes(fill = group),stat = "identity")+
  geom_text(aes(y = X.IncMSE + 1,label = label))+
  labs(title =expression("R"^"2"=="0.192  P=0.01 "))+
  theme(plot.title = element_text(size = 8))+
  labs(x = NULL, y = "PC2")+
  scale_fill_manual(values =  c( "#B0A875",'#AF7EC0'))+
  coord_flip()
p7

env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_WU.xlsx", sheet =1,rowNames = T,) 
env<-env[,-c(1,5,7)]
names(env)[c(7,18)]<-c("Amino_sugar","Lignin_phenol")
names(env)[c(10:17)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")
amino<- env %>%
  mutate(V = Vanillin+Acetovanillone+Vanillic_acid,S=Acetosyringone+Syringaldehyde+Syringic_acid,C=Hydroxycinnamic_acid+Ferulic_acid) %>%
  mutate(SV = S/V ,CV = C/V ,AdAlv=Vanillic_acid/Vanillin,AdAls=Syringic_acid/Syringaldehyde)%>%
  mutate(Lig_C = (Lignin_phenol/(1000*SOC))*100 ,V_C = (V/(1000*SOC))*100  ,S_C = (S/(1000*SOC))*100,C_C = (C/(1000*SOC))*100)%>%
  select(c("pH","SOC","TN","Amino_sugar","BNC","FNC","Lignin_phenol","V","S","C","AdAlv","AdAls"))
env <- amino
names(env)[c(4:12)]<-c("Amino_sugar","BNC","FNC","Lignin_phenols","Vanillyl_phenols","Syringyl_phenols","Cinnamyl_phenols","AdAl_V","AdAl_S")


mantel_funciton1 <- function(class){
  library(linkET)
  library(tidyverse)
  library(ggplot2)
  library(vegan)
  library(openxlsx)
  otu <- read.csv("D:\\Program Files (x86)\\RStudio\\example\\otutab_ITS.CSV",header=T,row.names=1)
  env<- read.xlsx("D:\\Program Files (x86)\\RStudio\\example\\FINAL\\ENV_WU.xlsx", sheet =1,rowNames = T,) 
  env<-env[,-c(1,5,7)]
  names(env)[c(7,18)]<-c("Amino_sugar","Lignin_phenol")
  names(env)[c(10:17)]<-c("Vanillin","Acetovanillone","Acetosyringone","Vanillic_acid","Syringaldehyde","Syringic_acid","Hydroxycinnamic_acid","Ferulic_acid")
  amino<- env %>%
    mutate(V = Vanillin+Acetovanillone+Vanillic_acid,S=Acetosyringone+Syringaldehyde+Syringic_acid,C=Hydroxycinnamic_acid+Ferulic_acid) %>%
    mutate(SV = S/V ,CV = C/V ,AdAlv=Vanillic_acid/Vanillin,AdAls=Syringic_acid/Syringaldehyde)%>%
    mutate(Lig_C = (Lignin_phenol/(1000*SOC))*100 ,V_C = (V/(1000*SOC))*100  ,S_C = (S/(1000*SOC))*100,C_C = (C/(1000*SOC))*100)%>%
    select(c("pH","SOC","TN","Amino_sugar","BNC","FNC","Lignin_phenol","V","S","C","AdAlv","AdAls"))
  env <- amino
  names(env)[c(4:12)]<-c("Amino_sugar","BNC","FNC","Lignin_phenols","Vanillyl_phenols","Syringyl_phenols","Cinnamyl_phenols","AdAl_V","AdAl_S")
  group=read.csv('D:\\Program Files (x86)\\RStudio\\example\\Group.csv',header=T,row.names = 1)
  group<-group[,-1]
  group<-group[which(group$Group!="Total"),]
  group$Group[which(group$Group=="Large macroaggregates")]<-"LMA+SMA"
  group$Group[which(group$Group=="Small macroaggregates")]<-"LMA+SMA"
  group$Group[which(group$Group=="Microaggregates")]<-"MIA+SCP"
  group$Group[which(group$Group=="Silt-clay particles")]<-"MIA+SCP"
  group$Group<- factor(group$Group,levels = c("LMA+SMA","MIA+SCP"),labels = c("LMA+SMA","MIA+SCP"))
  group<-filter(group,Group == class)
  otu <- as.data.frame(t(otu))
  otu <- otu[rownames(group), ]
  env <- env[rownames(group), ]
  dist.abund <- vegdist(otu, method = 'bray')
  dist.env1 <- map(env,dist, method = 'euclidean')
  mantel1 <- function(dist){
    result <- mantel(dist.abund, dist, method = 'spearman', permutations = 9999, na.rm = TRUE)
    r<- c(result[["statistic"]],result[["signif"]])
    return(r) }
  mantelr <- map(dist.env1,mantel1)
  mantel_result<- t(as.data.frame(mantelr))
  scale.env <- scale(env, center = TRUE, scale = TRUE)
  dist.env <- dist(scale.env, method = 'euclidean')
  mantel <- cbind(class,row.names(mantel_result),mantel_result) %>% as.data.frame
  names(mantel) <- c('spec',"env","r",'P')
  return(mantel)
}

library(parallel)
clnum <- detectCores()
cl<- makeCluster(2)
cname1 <-  c("LMA+SMA","MIA+SCP")
result1 <- parLapply(cl,cname1, mantel_funciton1)
stopCluster(cl)
mantel_s <- rbind(result1[[1]],result1[[2]])
mantel_sn <-  as.data.frame(lapply(mantel_s[,3:4],as.numeric))
mantel_s <- cbind(mantel_s[1:2],mantel_sn)
mantel2 <- mantel_s%>% 
  filter(env!="env") %>%
  filter(P < 0.05) %>% 
  mutate(rd = cut(r, breaks = c(-Inf, 0.4, 0.6, Inf),
                  labels = c("< 0.4", "0.4 - 0.6", ">= 0.6")),
         pd = cut(P, breaks = c(-Inf, 0.01, 0.05, Inf),
                  labels = c("< 0.01", "0.01 - 0.05", ">= 0.05")))
p8<-qcorrplot(correlate(env), type = "lower", diag = FALSE) +
  geom_square() +
  geom_couple(aes(colour = pd, size = rd), 
              data = mantel2, 
              curvature = nice_curvature()) +
  scale_size_manual(values = c(0.3, 0.6, 1.2)) +
  scale_colour_manual(values = color_pal(3)) +
  guides(size = guide_legend(title = "Mantel's r",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Mantel's p", 
                               override.aes = list(size = 3), 
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))+
  theme(axis.text=element_text(colour='black',size=9),legend.text=element_text(colour='black',size=7),legend.title = element_text(colour='black',size=8))+
  theme(panel.background = element_blank())+
  scale_fill_continuous(low = "#7A7EE3", high = "#DA6349")

p8








