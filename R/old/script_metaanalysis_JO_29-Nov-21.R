#####Charge packages#####
library (metafor)
library(here)
library (robumeta)
library (vegan)
library(tidyverse)
library(patchwork)

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

#Generating PCA axes
dat.macroinv_frame <- as.data.frame(dat.macroinv_filt[,c("Temp_media","Precipitacao","Temp_sazonal","Prec_sazonal")]) #climatic data to a data frame
pca_result_mac <- prcomp(dat.macroinv_frame, scale=TRUE) #PCA with climatic data
pca_result_mac
plot(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"])
text(pca_result_mac$x[,"PC1"], pca_result_mac$x[,"PC2"], dat.macroinv_filt[,"ID"])

PCA1_mac<- as.data.frame(pca_result_mac$x[,c("PC1", "PC2")]) #PC1 data frame
Clim_mac <- PCA1_mac #Climatic matrix
dat.macroinv_filt1 <- cbind(dat.macroinv_filt, Clim_mac) #Climatic matrix added to macroinv matrix

#Metaregression
colnames(dat.macroinv_filt1)
round(as.dist(cor(dat.macroinv_filt1[,c("Range", "Position", "Historic_land", "Topography", "PC1", "PC2")], use = "na.or.complete")),3) #Range and position have >0.7

boxplot(dat.macroinv_filt1$yi~dat.macroinv_filt1$Transition)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$Topography)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$Historic_land)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$PC1)
plot(dat.macroinv_filt1$yi~dat.macroinv_filt1$PC2)
#There were no evident non-linearities
       
model_macroinv_ecol<- robu(yi ~ Transition+Topography+
                                Historic_land+PC1+PC2,
                                data = dat.macroinv_filt1, 
                                modelweights = "CORR", studynum = ID,  
                                var.eff.size = vi , small = T)

model_macroinv_ecol

plot (yi ~ Historic_land, data = dat.macroinv_filt1)
#Ecological predictors

#Jean: Alternative graph showing the weights of each effect-size:
wi_fig <- 1/dat.macroinv_filt1$vi #Approximate weights (the true weight considered for analysis include the T2, mean variance by study and the correlation)
wi_fig_stan <- (wi_fig - min(wi_fig))/(max(wi_fig)-min(wi_fig))
wi_fig_stan <- (wi_fig_stan*3)+.5 #To make the range be between 0.5 and 3.5

par(las=1,mar=c(5,4.5,1,1)+.1)
plot (yi ~ Historic_land, data = dat.macroinv_filt1,
      cex=wi_fig_stan, ylab=expression(paste("Fisher's ",italic(z))))
#Mateus, this figure is interesting because it informs the weight of each effect-size on the meta-regression.
#Note that all effects with largest weights are more or less spread at a horizontal line, which explain a non-significant relationship of effect with historic land

#####
model_macroinv_method<- robu(yi ~ Range,
                                  data = dat.macroinv_filt1,
                                  modelweights = "CORR", studynum = ID,  
                                  var.eff.size = vi , small = T)

model_macroinv_method #Methodological predictor

plot (yi ~ Range, data = dat.macroinv_filt1)

summary(wi_fig)
summary(wi_fig_stan)

par(las=1,mar=c(5,4.5,1,1)+.1)
plot (yi ~ Range, data = dat.macroinv_filt1,cex=wi_fig_stan,
      ylab=expression(paste("Fisher's ",italic(z))))
legend(legend = c("5.0","52.0","682.0"),bty="n",
       x=0,y=-1,pt.cex = c(.5,.7082,3.5),pch=1)
text("Weigths",x = 10,y=-1)

#Fail Safe-number
fsn (yi=yi, vi=vi, data= dat.macroinv_filt1,type="Orwin")

#Jean: best to make publication bias analysis with mean z and Vz by study to avoid bias due to multiplicity:
mean_yi_macro <- as.numeric(unlist(by(dat.macroinv_filt1$yi,INDICES = dat.macroinv_filt1$ID,FUN = mean)))
mean_vi_macro <- as.numeric(unlist(by(dat.macroinv_filt1$vi,INDICES = dat.macroinv_filt1$ID,FUN = mean)))

fsn (yi=mean_yi_macro, vi=mean_vi_macro,type="Orwin")

#Olthers
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
dat.fish_filt1 <- dat.fish[dat.fish$Transition!="Openagr",]
dim(dat.fish_filt1)
length(unique(dat.fish_filt$ID)) #We lost 3 studies

#Testing colinearity of climate variables
colnames(data_fish)
correl_climate_fish <- cor(dat.fish_filt[,c("Temp_media","Precipitacao","Temp_sazonal","Prec_sazonal")]) 
as.dist(correl_climate_fish) # Precipitacao and Temp_sazonal will be excluded from the analisys
dat.fish_filt <- subset(dat.fish_filt, select = -c(Precipitacao,Temp_sazonal))
colnames(dat.fish_filt)

#Generating PCA axes

#dat.fish_frame <- as.data.frame(dat.fish_filt[,c("Temp_media","Precipitacao","Temp_sazonal","Prec_sazonal")]) #climatic data to a data frame
#pca_result_fish <- prcomp(dat.fish_frame, scale=TRUE) #PCA with climatic data
#pca_result_mac
#plot(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"])
#text(pca_result_fish$x[,"PC1"], pca_result_fish$x[,"PC2"], dat.fish_filt[,"ID"])
#PCA1_fish<- as.data.frame(pca_result_fish$x[,c("PC1", "PC2")]) #PC1 data frame
#Clim_fish <- PCA1_fish #Climatic matrix
#dat.fish_filt1 <- cbind(dat.fish_filt, Clim_fish) #Climatic matrix added to macroinv matrix

boxplot(dat.fish_filt1$yi~dat.fish_filt1$Transition)
plot(dat.fish_filt1$yi~dat.fish_filt1$Topography)
plot(dat.fish_filt1$yi~dat.fish_filt1$Historic_land)
plot(dat.fish_filt1$yi~dat.fish_filt1$Temp_media)
plot(dat.fish_filt1$yi~dat.fish_filt1$Prec_sazonal)
#There were no evident non-linearities

#Metaregression
colnames(dat.fish_filt1)
as.dist(cor(dat.macroinv_filt[,c("Range", "Position", "Historic_land", 
                                 "Topography", "Temp_media", "Prec_sazonal")], use = "na.or.complete")) #Range and position have >0.7


model_fish_ecol<- robu(yi ~ Transition+ Topography+Historic_land+
                            Temp_media+Prec_sazonal, #Mateus, keep the intercept in this model     
                            data = dat.fish_filt1, 
                            modelweights = "CORR", studynum = ID,  
                            var.eff.size = vi , small = T)

model_fish_ecol #Ecological predictors

plot (yi ~ Prec_sazonal, data = dat.fish_filt)

model_fish_method<- robu(yi ~ Range,
                              data = dat.fish_filt,
                              modelweights = "CORR", studynum = ID,
                              var.eff.size = vi , small = T)

model_fish_method #Methodological predictors

plot (yi ~ Range, data = dat.fish_filt)

#Safe-number
fsn (yi=yi, vi=vi, data= dat.fish_filt,type="Orwin")

#Jean: best to make publication bias analysis with mean z and Vz by study to avoid bias due to multiplicity:
mean_yi_fish <- as.numeric(unlist(by(dat.fish_filt$yi,INDICES = dat.fish_filt$ID,FUN = mean)))
mean_vi_fish <- as.numeric(unlist(by(dat.fish_filt$vi,INDICES = dat.fish_filt$ID,FUN = mean)))

fsn (yi=mean_yi_fish, vi=mean_vi_fish,type="Orwin")


#Others
write.table(dat.fish_filt, file= (here("output", "dados_plot_fish.csv")),sep = ",", quote = TRUE)

#####PLOTS####
###Effect Size

plot_macroinv <- dat.macroinv_filt1
plot_fish <- dat.fish_filt

vetor_macro <- c(0:80) #vetor para eixo X
vetor_fish <- c(0:35) #vetor para eixo X

macro <- ggplot(plot_macroinv,aes(x=yi, y=vetor_macro, fill=Transition))+
  theme_bw()+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  geom_vline(xintercept = 0,alpha = 0.8)+
  geom_vline(xintercept = -0.47, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  geom_vline(xintercept = -0.61, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Relations Ordered by Effect Size")+
  xlab("Effect Size (Fisher's z)")+
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
  #geom_vline(xintercept = 0.43, alpha = 0.8, colour = "#ff0000", linetype = "dashed")+
  #geom_vline(xintercept = 0.25, alpha = 0.8, colour = "#333333", linetype = "dashed")+
  geom_point(size=4, shape = 21,colour="black",alpha=0.8)+
  ylab("Relations Ordered by Effect Size")+
  xlab("Effect Size (Fisher's z)")+
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
text(-1.3, 0.02, "Macroinvertebrates")
funnel(fit.fish)
text(-1.3, 0.02, "Fishes")
dev.off()

#Jean: best to make publication bias analysis with mean z and Vz by study to avoid bias due to multiplicity:
fit.macroinv <- rma (mean_yi_macro, mean_vi_macro)
fit.fish <-  rma (mean_yi_fish,mean_vi_fish)

png(here("output", "funnel.png"), res=300,width=3200,height=1300)
par(mfrow=c(1,2))
funnel(fit.macroinv)
text(-1.3, 0.02, "Macroinvertebrates")
funnel(fit.fish)
text(-1.3, 0.02, "Fishes")
dev.off()


