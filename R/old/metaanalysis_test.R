#####Metaanalysis test#####
library (metafor)
library(here)
library (robumeta)
library (vegan)
library(tidyverse)
library(patchwork)

####Meta-analysis and Meta-regression#### 
####MMACROINVERTEBRATES####
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
dat.macroinv_filt1 <- cbind(dat.macroinv_filt, PCA1) #adicionando à matriz de análise

#Metaregression

model_macroinv_redu_covars<- robu(yi ~ transition+Clim_mac-27,
                                  data = dat.macroinv_filt1, 
                                  modelweights = "CORR", studynum = ID,  
                                  var.eff.size = vi , small = T)
model_macroinv_redu_covars

##Funnel Plot
fit.macroinv <- rma(yi, vi,
                    data= dat.macroinv_filt1, method = "DL")

funnel(fit.macroinv)

##Polygon plot
transition<-as.numeric(as.factor(dat.macroinv_filt1$transition))
y <- c(0:76)
length (y)
x11()
plot( dat.macroinv_filt1$yi, main = "Macroinvertebrates", ylab = "Effect Size (Hedge's H)", xlab = "Studies Ordered by Effect Size", type = 'p', 
      frame = F,pch =18, col = transition, cex = 2, ylim = c(-1.5, 0.5 ), xlim = c(0,80))
abline (0, 0)
abline(-0.27, 0, lty=2, col =2)
abline (-0.43, 0, lty=2)
text(70,-0.21 ,lab= "Forest->Agriculture", col = 2, font=2)
text(6,-0.46 ,lab= "Forest->Urban", col = 1, font =2)

#Olthers
write.table(dat.macroinv_filt1, file= (here("data","processed", "dados_plot.csv")),sep = ",", quote = TRUE)

####FISHES####
#All values in one model#
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

model_fish_redu_covars<- robu(yi ~ transition+Clim_fish-1,
                                  data = dat.fish_filt1, 
                                  modelweights = "CORR", studynum = ID,  
                                  var.eff.size = vi , small = T)
model_fish_redu_covars

##Funnel Plot
fit.fish <- rma(yi, vi,
                    data= dat.fish_filt1, method = "DL")

funnel(fit.fish)

#Olthers
write.table(dat.fish_filt1, file= (here("data","processed", "dados_plot_fish.csv")),sep = ",", quote = TRUE)
x11()
forest (fit.fish)

#####PLOTS####
##Effect Size

plot_macroinv <- dat.macroinv_filt1
plot_fish <- dat.fish_filt1

vetor_macro <- c(0:76) #vetor para eixo X
vetor_fish <- c(0:36) #vetor para eixo X

macro <- ggplot(plot_macroinv,aes(x=yi, y=vetor_macro, fill=transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = -0.27, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.43, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Studies Ordered by Effect Size")+
  xlab("Effect Size (Hedge's H)")+
  xlim(-1.5,1.0)+
  coord_flip()+
  annotate("text",x=0.94,y=10,label="Macroinvertebrate",family="Noto Sans",size=4.5,fontface="bold")+
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
  geom_vline(xintercept = 0.0578, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.2208, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Studies Ordered by Effect Size")+
  xlab("Effect Size (Hedge's H)")+
  xlim(-1.5,1.0)+
  coord_flip()+
  annotate("text",x=0.94,y=0,label="Fish",family="Noto Sans",size=4.5,fontface="bold")+
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
png(here("output", "funnel.png"), res=300,width=3200,height=1300)
par(mfrow=c(1,2))
funnel(fit.macroinv)
funnel(fit.fish)
dev.off()





