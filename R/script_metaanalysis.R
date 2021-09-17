#####Charge packages#####
library (metafor)
library(here)
library (robumeta)
library (vegan)
library(tidyverse)
library(patchwork)

####Meta-analysis and Meta-regression 
####MACROINVERTEBRATES####
data_macroinv <- read.table(here ("data","processed","macroinvertebrates.txt"), h=T) ##All data
colnames(data_macroinv)
str(data_macroinv)

##Modelling de Hedge's H

dat.macroinv<- escalc(measure="ZCOR", ri=Effect_size, ni=N_sample, vtype = "LS",
                     data=data_macroinv, append=TRUE) #método para valores de correlação, porém, extraindo o valor de z de fischer, já que valor do r (vr = ((1-r^2)^2)/(n-1))


##Calculating the cummulative effects of each study
cum_effect <- robu(yi ~1,
                   data = dat.macroinv, 
                   modelweights = "CORR", #Jean: Use os pesos correlacionados ao inves dos hierarquicos
                   studynum = ID,  
                   var.eff.size = vi , small = T)
cum_effect

#Excluding relationships with small number
dim(dat.macroinv)
length(unique(dat.macroinv$ID)) #Number of studies
dat.macroinv_filt <- dat.macroinv[dat.macroinv$Transition!="Openagr",]
dat.macroinv_filt <- dat.macroinv_filt[dat.macroinv_filt$Transition!="Openurb",]
dim(dat.macroinv_filt)
length(unique(dat.macroinv_filt$ID)) #We lost 5 studies

#Testing colinearity of climate variables
correl_climate_mac <- cor(dat.macroinv_filt[,c(6:9)]) 
correl_climate_mac

#Generating PCA axes
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c(6:9)]) #criando um dataframe com as var climáticas
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #fazendo a PCA
pca_result_mac
plot(pca_result_mac$x, pca_result_mac$y)
PCA1_mac<- as.data.frame(pca_result_mac$x) #Gerando uma matriz com os PC
Clim_mac <- PCA1_mac[,1] #Gerando uma matriz apenas com o PCA1
dat.macroinv_filt1 <- cbind(dat.macroinv_filt, Clim_mac) #adicionando à matriz de análise

#Metaregression
cor(dat.macroinv_filt1[,c(10:13,17)], use = "na.or.complete") #Range and position have >0.7
colnames(dat.macroinv_filt1)

model_macroinv_ecol<- robu(yi ~ Transition+Clim_mac+Topography+
                                            Historic_land,
                                            data = dat.macroinv_filt1, 
                                            modelweights = "CORR", studynum = ID,  
                                            var.eff.size = vi , small = T)
model_macroinv_ecol

model_macroinv_method<- robu(yi ~ Range,
                             data = dat.macroinv_filt1,
                             modelweights = "CORR", studynum = ID,  
                             var.eff.size = vi , small = T)
model_macroinv_method

#Safe-number
fsn (yi=yi, vi=vi, data= dat.macroinv_filt1,type="Orwin")

#Olthers
write.table(dat.macroinv_filt1, file= (here("output", "macroinv_filt.csv")),sep = ",", quote = TRUE,  row.names=F)

####FISHES####

data_fish <- read.table(here ("data","processed","fishes.txt"), h=T) ##Charge data
colnames(data_fish)
str(data_fish)

##Modelling de Hedge's H

dat.fish<- escalc(measure="ZCOR", ri=Effect_size, ni=N_sample, vtype = "LS",
                      data=data_fish, append=TRUE) #método para valores de correlação, porém, extraindo o valor de z de fischer, já que valor do r (vr = ((1-r^2)^2)/(n-1))


##Calculating the cummulative effects of each study
cum_effect <- robu(yi ~1,
                   data = dat.fish, 
                   modelweights = "CORR", #Jean: Use os pesos correlacionados ao inves dos hierarquicos
                   studynum = ID,  
                   var.eff.size = vi , small = T)
cum_effect

#Excluding relationships with small number
dim(dat.fish)
length(unique(dat.fish$ID)) #Number of studies
dat.fish_filt <- dat.fish[dat.fish$Transition!="openagr",]
dim(dat.fish_filt)
length(unique(dat.fish_filt$ID)) #We lost 3 studies

#Testing colinearity of climate variables
correl_climate_fish <- cor(dat.fish_filt[,c(6:9)])
correl_climate_fish
dat.fish_filt <- dat.fish_filt[, -8] # Excluding temperature sazonality (r= -0.86)
dat.fish_filt <- dat.fish_filt[, -7] # Excluding temperature sazonality (r= -0.86)
colnames (dat.fish_filt)

#Generating PCA axes
dat.fish_frame <- as.data.frame(dat.fish_filt[,c(6:7)]) #criando um dataframe com as var climáticas
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #fazendo a PCA
pca_result_fish
plot(pca_result_fish$x, pca_result_fish$y)
PCA1_fish<- as.data.frame(pca_result_fish$x) #Gerando uma matriz com os PC
Clim_fish <- PCA1_fish[,1] #Gerando uma matriz apenas com o PCA1
dat.fish_filt1 <- cbind(dat.fish_filt, Clim_fish) #adicionando à matriz de análise


#Metaregression
cor(dat.fish_filt1[,c(8:11,15)], use = "na.or.complete") #Range and position have >0.7
colnames(dat.macroinv_filt1)

model_fish_ecol<- robu(yi ~ Transition+Clim_fish+
                      Topography+Historic_land,
                      data = dat.fish_filt1, 
                      modelweights = "CORR", studynum = ID,  
                      var.eff.size = vi , small = T)

model_fish_ecol


model_fish_method<- robu(yi ~ Range,
                         data = dat.fish_filt1,
                         modelweights = "CORR", studynum = ID,
                         var.eff.size = vi , small = T)

model_fish_method


=#Safe-number
fsn (yi=yi, vi=vi, data= dat.fish_filt1,type="Orwin")

#Others
write.table(dat.fish_filt1, file= (here("output", "dados_plot_fish.csv")),sep = ",", quote = TRUE)

#####PLOTS####
###Effect Size

plot_macroinv <- dat.macroinv_filt1
plot_fish <- dat.fish_filt1

vetor_macro <- c(0:82) #vetor para eixo X
vetor_fish <- c(0:41) #vetor para eixo X

macro <- ggplot(plot_macroinv,aes(x=yi, y=vetor_macro, fill=Transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = 0.46, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = 0.10, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Relations Ordered by Effect Size")+
  xlab("Effect Size (Hedges' g)")+
  xlim(-1.5,1.0)+
  coord_flip()+
  annotate("text",x=0.94,y=10,label="Macroinvertebrates",family="Noto Sans",size=4.5,fontface="bold")+
  scale_fill_manual(name = "Transition Category",
                    values = c("#ff0000","#333333"),
                    labels = c("Forest to Agriculture","Forest to Urban"))+
  theme(axis.title.y = element_text(colour = "black",size=12,
                                    face="bold",family="Noto Sans"))+
  theme(axis.text.y = element_text(angle=90,size=9,
                                   hjust=0.2,family="Noto Sans"))+
  theme(axis.title.x = element_text(colour = "black",size=12,
                                    face="bold",family="Noto Sans"))+
  theme(axis.text.x = element_text(hjust=0.2,size=9,family="Noto Sans"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=11,family="Noto Sans"))+
  theme(legend.title = element_text(size=12,family="Noto Sans", face="bold"))

fish <- ggplot(plot_fish,aes(x=yi, y=vetor_fish, fill=Transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = -0.06, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.36, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Relations Ordered by Effect Size")+
  xlab("Effect Size (Hedges' g)")+
  xlim(-1.5,1.0)+
  coord_flip()+
  annotate("text",x=0.94,y=1,label="Fishes",family="Noto Sans",size=4.5,fontface="bold")+
  scale_fill_manual(name = "Transition Category",
                    values = c("#ff0000","#333333"),
                    labels = c("Forest to Agriculture","Forest to Urban"))+
  theme(axis.title.y = element_text(colour = "black",size=12,
                                    face="bold",family="Noto Sans"))+
  theme(axis.text.y = element_text(angle=90,size=9,
                                   hjust=0.2,family="Noto Sans"))+
  theme(axis.title.x = element_text(colour = "black",size=12,
                                    face="bold",family="Noto Sans"))+
  theme(axis.text.x = element_text(hjust=0.2,size=9,family="Noto Sans"))+
  theme(legend.position = "bottom")+
  theme(legend.text = element_text(size=11,family="Noto Sans"))+
  theme(legend.title = element_text(size=12,family="Noto Sans", face="bold"))

(macro|fish)+ 
  patchwork::plot_layout(guides="collect") & ggplot2::theme(legend.position = "bottom")

#Exporting plot
ggsave (here("output", "metaanal.png"), width = 12, height = 5)

###Funnel plot
fit.macroinv <- rma (dat.macroinv, yi, vi)
fit.fish <-  rma (dat.fish, yi, vi)

png(here("output", "funnel.png"), res=300,width=3200,height=1300)
par(mfrow=c(1,2))
funnel(fit.macroinv)
text(-1.15, 0.02, "Macroinvertebrates")
funnel(fit.fish)
text(-1.3, 0.02, "Fishes")

dev.off()




