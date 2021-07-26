#####Charge packages#####
library (metafor)
library(here)
library (robumeta)
library (vegan)
library(tidyverse)
library(patchwork)

####Meta-analysis and Meta-regression#### 
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
dat.macroinv_filt <- dat.macroinv[dat.macroinv$transition!="openagr",]
dat.macroinv_filt <- dat.macroinv_filt[dat.macroinv_filt$transition!="openurb",]
dim(dat.macroinv_filt)
length(unique(dat.macroinv_filt$ID)) #We lost 5 studies

#Testing colinearity of climate variables
correl_climate_mac <- cor(dat.macroinv_filt[,c(24:27)]) 

#Generating PCA axes
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c(24:27, 29)]) #criando um dataframe com as var climáticas
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #fazendo a PCA
PCA1_mac<- as.data.frame(pca_result_mac$x) #Gerando uma matriz com os PC
Clim_mac <- PCA1_mac[,1] #Gerando uma matriz apenas com o PCA1
dat.macroinv_filt1 <- cbind(dat.macroinv_filt, Clim_mac) #adicionando à matriz de análise

#Metaregression

model_macroinv_redu_covars<- robu(yi ~ transition+Clim_mac,
                                  data = dat.macroinv_filt1, 
                                  modelweights = "CORR", studynum = ID,  
                                  var.eff.size = vi , small = T)
model_macroinv_redu_covars
#Safe-number
fsn (yi=yi, vi=vi, data= dat.macroinv_filt1,type="Orwin")

#Olthers
write.table(dat.macroinv_filt1, file= (here("data","processed", "dados_plot.csv")),sep = ",", quote = TRUE)

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
dat.fish_filt <- dat.fish[dat.fish$transition!="openagr",]
dim(dat.fish_filt)
length(unique(dat.fish_filt$ID)) #We lost 1 study

#Testing colinearity of climate variables
correl_climate_fish <- cor(dat.fish_filt[,c(24:27)])
correl_climate_fish
dat.fish <- dat.fish[, -25] # Excluding temperature sazonality (r= -0.86)

#Generating PCA axes
dat.fish_frame <- as.data.frame(dat.fish_filt[,c(24:26, 29)]) #criando um dataframe com as var climáticas
pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #fazendo a PCA
PCA1_fish<- as.data.frame(pca_result_fish$x) #Gerando uma matriz com os PC
Clim_fish <- PCA1_fish[,1] #Gerando uma matriz apenas com o PCA1
dat.fish_filt1 <- cbind(dat.fish_filt, Clim_fish) #adicionando à matriz de análise

#Metaregression

model_fish_redu_covars<- robu(yi ~ transition+Clim_fish,
                                  data = dat.fish_filt1, 
                                  modelweights = "CORR", studynum = ID,  
                                  var.eff.size = vi , small = T)
model_fish_redu_covars

#Safe-number
fsn (yi=yi, vi=vi, data= dat.fish_filt1,type="Orwin")

#Others
write.table(dat.fish_filt1, file= (here("data","processed", "dados_plot_fish.csv")),sep = ",", quote = TRUE)
x11()
forest (fit.fish)

#####PLOTS####
##Effect Size

plot_macroinv <- dat.macroinv_filt1
plot_fish <- dat.fish_filt1

vetor_macro <- c(0:82) #vetor para eixo X
vetor_fish <- c(0:41) #vetor para eixo X

macro <- ggplot(plot_macroinv,aes(x=yi, y=vetor_macro, fill=transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = -0.23, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.39, alpha = 0.8, colour = "#333333", linetype = "dashed")+
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

fish <- ggplot(plot_fish,aes(x=yi, y=vetor_fish, fill=transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = 0.06, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.25, alpha = 0.8, colour = "#333333", linetype = "dashed")+
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
ggsave (here("output", "metaanal.png"), width = 12, height = 5 )

##Funnel plot
fit.macroinv <- rma (dat.macroinv, yi, vi)
fit.fish <-  rma (dat.fish, yi, vi)



png(here("output", "funnel.png"), res=300,width=3200,height=1300)
par(mfrow=c(1,2))
funnel(fit.macroinv)
text(-1.15, 0.02, "Macroinvertebrates")
funnel(fit.fish)
text(-1.3, 0.02, "Fishes")

dev.off()






