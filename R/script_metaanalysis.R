#####Charge packages#####
library(here)
library (metafor)
library (robumeta)
library (vegan)
library(tidyverse)
library(patchwork)
library(lme4)


####Meta-analysis and Meta-regression 
####MACROINVERTEBRATES####
data_macroinv <- read.table(here ("data","processed","macroinvertebrates.txt"), h=T) ##Macroinvertebrate data
colnames(data_macroinv)
str(data_macroinv)

##Modelling de Hedge's H (Fischer Z)

dat.macroinv<- escalc(measure="ZCOR", ri=Effect_size, ni=N_sample, vtype = "LS",
                     data=data_macroinv, append=TRUE) #método para valores de correlação, porém, extraindo o valor de z de fischer, já que valor do r (vr = ((1-r^2)^2)/(n-1))

#Excluding relationships with small number 
dim(dat.macroinv) #initial dimension 
length(unique(dat.macroinv$ID)) #Number of studies
dat.macroinv_filt <- dat.macroinv[dat.macroinv$Transition!="Openagr",]
dat.macroinv_filt <- dat.macroinv_filt[dat.macroinv_filt$Transition!="Openurb",]
dim(dat.macroinv_filt)
length(unique(dat.macroinv_filt$ID)) #We lost 5 studies

#Testing colinearity of climate variables
colnames(data_macroinv)
correl_climate_mac <- cor(dat.macroinv_filt[,c("Temp_media","Precipitacao","Temp_sazonal","Prec_sazonal")]) 
as.dist(correl_climate_mac)

##Now, we will make 4 metaregressions, one for each spatial scale
#25km buffer

#Generating PCA axes
colnames(data_macroinv)
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c("Temp_media_25","Precipitacao_25","Temp_sazonal_25","Prec_sazonal_25")]) #climatic data to a data frame
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #PCA with climatic data
pca_result_mac
plot(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"])
text(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"], dat.macroinv_filt[,"ID"])

PCA1_mac_25<- as.data.frame(pca_result_mac$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_mac_25 <- PCA1_mac_25 #Climatic matrix
dat.macroinv_filt_25 <- cbind(dat.macroinv_filt, Clim_mac_25) #Climatic matrix added to macroinv matrix

#Metaregression
colnames(dat.macroinv_filt_25)
round(as.dist(cor(dat.macroinv_filt1[,c("Range", "Position", "LUI_25", "Topography", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_macroinv_25<- robu(yi ~ Transition+LUI_25+PC1+PC2+Range,
                                data = dat.macroinv_filt1, 
                                modelweights = "CORR", studynum = ID,  
                                var.eff.size = vi , small = T)

model_macroinv_25

#50k buffer

#Generating PCA axes
colnames(data_macroinv)
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c("Temp_media_50","Precipitacao_50","Temp_sazonal_50","Prec_sazonal_50")]) #climatic data to a data frame
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #PCA with climatic data
pca_result_mac
plot(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"])
text(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"], dat.macroinv_filt[,"ID"])

PCA1_mac_50<- as.data.frame(pca_result_mac$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_mac_50 <- PCA1_mac_50 #Climatic matrix
dat.macroinv_filt_50 <- cbind(dat.macroinv_filt, Clim_mac_50) #Climatic matrix added to macroinv matrix

#Metaregression
colnames(dat.macroinv_filt_50)
round(as.dist(cor(dat.macroinv_filt1[,c("Range", "Position", "LUI_50", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_macroinv_50<- robu(yi ~ Transition+LUI_50+PC1+PC2+Range,
                         data = dat.macroinv_filt1, 
                         modelweights = "CORR", studynum = ID,  
                         var.eff.size = vi , small = T)

model_macroinv_50


#75k buffer

#Generating PCA axes
colnames(data_macroinv)
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c("Temp_media_75","Precipitacao_75","Temp_sazonal_75","Prec_sazonal_75")]) #climatic data to a data frame
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #PCA with climatic data
pca_result_mac
plot(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"])
text(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"], dat.macroinv_filt[,"ID"])

PCA1_mac_75<- as.data.frame(pca_result_mac$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_mac_75 <- PCA1_mac_75 #Climatic matrix
dat.macroinv_filt_75 <- cbind(dat.macroinv_filt, Clim_mac_75) #Climatic matrix added to macroinv matrix

#Metaregression
colnames(dat.macroinv_filt_75)
round(as.dist(cor(dat.macroinv_filt1[,c("Range", "Position", "LUI_75", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_macroinv_75<- robu(yi ~ Transition+LUI_75+PC1+PC2+Range,
                         data = dat.macroinv_filt1, 
                         modelweights = "CORR", studynum = ID,  
                         var.eff.size = vi , small = T)

model_macroinv_75

#100k buffer

#Generating PCA axes
colnames(data_macroinv)
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c("Temp_media_100","Precipitacao_100","Temp_sazonal_100","Prec_sazonal_100")]) #climatic data to a data frame
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #PCA with climatic data
pca_result_mac
plot(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"])
text(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"], dat.macroinv_filt[,"ID"])

PCA1_mac_100<- as.data.frame(pca_result_mac$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_mac_100 <- PCA1_mac_100 #Climatic matrix
dat.macroinv_filt_100 <- cbind(dat.macroinv_filt, Clim_mac_100) #Climatic matrix added to macroinv matrix

#Metaregression
colnames(dat.macroinv_filt_100)
round(as.dist(cor(dat.macroinv_filt1[,c("Range", "Position", "LUI_100", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_macroinv_100<- robu(yi ~ Transition+LUI_100+PC1+PC2+Range,
                         data = dat.macroinv_filt1, 
                         modelweights = "CORR", studynum = ID,  
                         var.eff.size = vi , small = T)

model_macroinv_100





boxplot(dat.macroinv_filt1$yi~dat.macroinv_filt1$Transition)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$Topography)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$LUI)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$ADR)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$PC1)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$PC2)
plot (yi ~ Range, data = dat.macroinv_filt1)

#Safe-number
mean_yi_macro <- as.numeric(unlist(by(dat.macroinv_filt1$yi,INDICES = dat.macroinv_filt1$ID,FUN = mean)))
mean_vi_macro <- as.numeric(unlist(by(dat.macroinv_filt1$vi,INDICES = dat.macroinv_filt1$ID,FUN = mean)))

fsn (yi=mean_yi_macro, vi=mean_vi_macro,type="Orwin")

#Others
write.table(dat.macroinv_filt1, file= (here("output", "macroinv_filt.csv")),sep = ",", quote = TRUE,  row.names=F)

####FISHES####

data_fish <- read.table(here ("data","processed","fishes.txt"), h=T) ##Macroinvertebrate data
colnames(data_fish)
str(data_fish)

##Modelling de Hedge's H (Fisher's z)

dat.fish<- escalc(measure="ZCOR", ri=Effect_size, ni=N_sample, vtype = "LS",
                      data=data_fish, append=TRUE) #método para valores de correlação, porém, extraindo o valor de z de fischer, já que valor do r (vr = ((1-r^2)^2)/(n-1))

#Excluding relationships with small number
dim(dat.fish) #initial dimension 
length(unique(dat.fish$ID)) #Number of studies
dat.fish_filt <- dat.fish[dat.fish$Transition!="Openagr",]
dat.fish_filt <- dat.fish_filt[dat.fish_filt$Transition!="Openurb",]
dim(dat.fish_filt)
length(unique(dat.fish_filt$ID)) #We lost 3 studies

#Testing colinearity of climate variables
colnames(data_fish)
correl_climate_fish <- cor(dat.fish_filt[,c("Temp_media","Precipitacao","Temp_sazonal","Prec_sazonal")]) 
as.dist(correl_climate_fish) # Precipitacao and Temp_sazonal will be excluded from the analisys
######dat.fish_filt1 <- subset(dat.fish_filt, select = -c(Precipitacao,Temp_sazonal))
#####colnames(dat.fish_filt1)

##Now, we will make 4 metaregressions, one for each spatial scale
#25km buffer

#Generating PCA axes
colnames(data_fish)
dat.fish_frame <- as.data.frame(dat.fish_filt[,c("Temp_media_25","Precipitacao_25","Temp_sazonal_25","Prec_sazonal_25")]) #climatic data to a data frame
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #PCA with climatic data
pca_result_fish
plot(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"])
text(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"], dat.fish_filt[,"ID"])

PCA1_fish_25<- as.data.frame(pca_result_fish$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_fish_25 <- PCA1_fish_25 #Climatic matrix
dat.fish_filt_25 <- cbind(dat.fish_filt, Clim_fish_25) #Climatic matrix added to fish matrix

#Metaregression
colnames(dat.fish_filt_25)
head(dat.fish_filt_25)
round(as.dist(cor(dat.fish_filt_25[,c("Range", "Position", "LUI_25", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_fish_25<- robu(yi ~ Transition+LUI_25+PC1+PC2+Range,
                         data = dat.fish_filt_25, 
                         modelweights = "CORR", studynum = ID,  
                         var.eff.size = vi , small = T)

model_fish_25

#50km buffer

#Generating PCA axes
colnames(data_fish)
dat.fish_frame <- as.data.frame(dat.fish_filt[,c("Temp_media_50","Precipitacao_50","Temp_sazonal_50","Prec_sazonal_50")]) #climatic data to a data frame
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #PCA with climatic data
pca_result_fish
plot(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"])
text(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"], dat.fish_filt[,"ID"])

PCA1_fish_50<- as.data.frame(pca_result_fish$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_fish_50 <- PCA1_fish_50 #Climatic matrix
dat.fish_filt_50 <- cbind(dat.fish_filt, Clim_fish_50) #Climatic matrix added to fish matrix

#Metaregression
colnames(dat.fish_filt_50)
head(dat.fish_filt_50)
round(as.dist(cor(dat.fish_filt_50[,c("Range", "Position", "LUI_50", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_fish_50<- robu(yi ~ Transition+LUI_50+PC1+PC2+Range,
                     data = dat.fish_filt_50, 
                     modelweights = "CORR", studynum = ID,  
                     var.eff.size = vi , small = T)

model_fish_50

#75km buffer

#Generating PCA axes
colnames(data_fish)
dat.fish_frame <- as.data.frame(dat.fish_filt[,c("Temp_media_75","Precipitacao_75","Temp_sazonal_75","Prec_sazonal_75")]) #climatic data to a data frame
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #PCA with climatic data
pca_result_fish
plot(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"])
text(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"], dat.fish_filt[,"ID"])

PCA1_fish_75<- as.data.frame(pca_result_fish$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_fish_75 <- PCA1_fish_75 #Climatic matrix
dat.fish_filt_75 <- cbind(dat.fish_filt, Clim_fish_75) #Climatic matrix added to fish matrix

#Metaregression
colnames(dat.fish_filt_75)
head(dat.fish_filt_75)
round(as.dist(cor(dat.fish_filt_75[,c("Range", "Position", "LUI_75", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_fish_75<- robu(yi ~ Transition+LUI_75+PC1+PC2+Range,
                     data = dat.fish_filt_75, 
                     modelweights = "CORR", studynum = ID,  
                     var.eff.size = vi , small = T)

model_fish_75

#100km buffer

#Generating PCA axes
colnames(data_fish)
dat.fish_frame <- as.data.frame(dat.fish_filt[,c("Temp_media_100","Precipitacao_100","Temp_sazonal_100","Prec_sazonal_100")]) #climatic data to a data frame
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #PCA with climatic data
pca_result_fish
plot(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"])
text(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"], dat.fish_filt[,"ID"])

PCA1_fish_100<- as.data.frame(pca_result_fish$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_fish_100 <- PCA1_fish_100 #Climatic matrix
dat.fish_filt_100 <- cbind(dat.fish_filt, Clim_fish_100) #Climatic matrix added to fish matrix

#Metaregression
colnames(dat.fish_filt_100)
head(dat.fish_filt_100)
round(as.dist(cor(dat.fish_filt_100[,c("Range", "Position", "LUI_100", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

model_fish_100<- robu(yi ~ Transition+LUI_100+PC1+PC2+Range,
                     data = dat.fish_filt_100, 
                     modelweights = "CORR", studynum = ID,  
                     var.eff.size = vi , small = T)

model_fish_100



boxplot(dat.fish_filt1$yi~dat.fish_filt1$Transition)
plot(dat.fish_filt1$yi~dat.fish_filt1$Topography)
plot(dat.fish_filt1$yi~dat.fish_filt1$LUI)
plot(dat.fish_filt1$yi~dat.fish_filt1$Temp_media)
plot(dat.fish_filt1$yi~dat.fish_filt1$Prec_sazonal)
plot(dat.fish_filt1$yi~dat.fish_filt1$Range)

#Safe-number
mean_yi_fish <- as.numeric(unlist(by(dat.fish_filt$yi,INDICES = dat.fish_filt$ID,FUN = mean)))
mean_vi_fish <- as.numeric(unlist(by(dat.fish_filt$vi,INDICES = dat.fish_filt$ID,FUN = mean)))

fsn (yi=mean_yi_fish, vi=mean_vi_fish,type="Orwin")

#Others
write.table(dat.fish_filt, file= (here("output", "dados_plot_fish.csv")),sep = ",", quote = TRUE)

#####PLOTS####
###Effect Size Transition
data_ <- as_tibble(dat.macroinv_filt1)

theme_set(theme_minimal(base_size = 15, base_family = "Arial"))

theme_update(
  panel.grid.major = element_line(color = "grey92", size = .4),
  panel.grid.minor = element_blank(),
  axis.title.x = element_text(color = "black", margin = margin(t = 7)),
  axis.title.y = element_text(color = "black", margin = margin(r = 7)),
  axis.text = element_text(color = "black"),
  axis.ticks =  element_line(color = "grey92", size = .4),
  axis.ticks.length = unit(.6, "lines"),
  legend.position = "top",
  plot.subtitle = element_text(hjust = 0.5, face = "bold", color = "black",
                               family = "Arial", 
                               size = 13, margin = margin(0, 0, 10, 0)),
  plot.margin = margin(rep(20, 4)),
  axis.line = element_line(colour="black")
)

macro_plot <-ggplot(dat.macroinv_filt1,aes(x=Transition, y=yi, size = vi))+
  geom_point(shape =21, alpha =.7,
             color = "black", fill = "grey40",
             position = position_jitter(w=.25, seed =1))+
  ylim(-1.6,1.3)+
  scale_size(range = c(7.3,2.5), guide = F)+
  labs(y = "Effect Size (Fisher's Z)", x = "\nTransition categories",
       subtitle = "Macroinvertebrates")+
  scale_x_discrete(labels =c("Forest to Agriculture", "Forest to Urban"))+
  geom_hline(yintercept = 0,alpha = 0.8, linetype = "dashed", colour = "grey50")+
  annotate("rect",xmin = 0.71, xmax = 1.30, ymin = -0.623, ymax= 0.210, 
           fill = "grey45", alpha = 0.2)+
  annotate("rect",xmin = 1.74, xmax = 2.27, ymin = -0.660, ymax = 0.177,
           fill = "grey45", alpha = 0.2)+
  annotate("segment",x = 0.71, y = -0.20631, xend = 1.30,  yend = -0.20631, size = 1)+
  annotate("segment",x = 1.74, y = -0.24180, xend = 2.27, yend = -0.24180, size = 1)


fish_plot <- ggplot(dat.fish_filt1,aes(x=Transition, y=yi, size = vi))+
  geom_point(shape =21, alpha =.7,
             color = "black", fill = "grey40",
             position = position_jitter(w=.25, seed =2))+
  ylim(-1.6,1.3)+
  scale_size(range = c(7.3,2.5), guide = F)+
  labs(x = "\nTransition categories",subtitle = "Fishes")+
  scale_x_discrete(labels =c("Forest to Agriculture", "Forest to Urban"))+
  geom_hline(yintercept = 0,alpha = 0.8, linetype = "dashed", colour = "grey50")+
  annotate("rect",xmin = 0.71, xmax = 1.30, ymin = -1.483, ymax= 0.351, 
           fill = "grey45", alpha = 0.2)+
  annotate("rect",xmin = 1.70, xmax = 2.27, ymin = -1.5, ymax = 0.746,
           fill = "grey45", alpha = 0.2)+
  annotate("segment",x = 0.71, y = -0.566, xend = 1.30,  yend = -0.566, size = 1)+
  annotate("segment",x = 1.70, y = -0.492, xend = 2.27, yend = -0.492, size = 1)+
  theme(axis.line.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y= element_blank())


p <- macro_plot + fish_plot


gridExtra::grid.arrange(egg::set_panel_size(p=p,
                                            width=unit(10, "in"), height=unit(7, "cm")))

ggsave (here("output", "Fig4.png"), width = 10, height = 7, dpi =600)

#Range plot
dat.macroinv_filt1

relacao_plot <- ggplot(dat.macroinv_filt1,aes(x=Range, y=yi, size = vi))+
  geom_smooth(method = "lm", color = "black")+
  geom_point(shape =21, alpha =.7,
             color = "black", fill = "grey40",
             position = position_jitter(w=.25, seed =3))+
  scale_size(range = c(7.3,2.5), guide = F)+
  labs(y = "Effect Size (Fisher's Z)", x = "\nRange (%)")


ggsave (here("output", "fig5.png"), width = 8, height = 7, dpi = 600)

###Funnel plot
fit.macroinv <- rma (dat.macroinv, yi, vi)
fit.fish <-  rma (dat.fish, yi, vi)

png(here("output", "funnel.png"), res=300,width=3200,height=1300)
par(mfrow=c(1,2))
funnel(fit.macroinv)
text(-1.3, 0.02, "Macroinvertebrates")
funnel(fit.fish)
text(-1.3, 0.02, "Fishes")
dev.off()

###OLD####
###Testing linear relationship betwen citations and Fishers'z

data_cit <- read.table(here ("data","processed","citations.txt"), h=T) ##Citation data
str (data_cit)
colnames(data_cit)

correl_cit <- cor(data_cit[,c("Cit_SCOPUS","yi")]) 

test_citations <- lm (Cit_SCOPUS ~ yi-1,data = data_cit)
summary(test_citations)
plot(data_cit$yi~data_cit$Cit_SCOPUS)


