#####Metaanalysis test#####
library (metafor)
library(here)
library(robumeta)

####First, Ill make a more simple model####
data_test <- read.table(here ("data","processed","data_test.txt"), h=T) ##All data
colnames(data_test)
str(data_test)
dat.test <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                   data=data_test, append=TRUE) #método para valores de correlação

res.fe1 <- rma(yi, vi, slab = paste(N),
               data= dat.test, method = "DL")
x11()
forest(res.fe1, atransf = exp)
text(-3.5, 29, "Article ID", pos = 4, font = 2)
text(3, 29, "Effect Size [95% CI]", pos = 2, font = 2)

fsn (yi, vi, data=dat.test, type="Orwin", target= -0.132) ## Para obter quantos estudos ainda precisariam para deixar de ser significativo

x11()
funnel (res.fe1)

###Now, Ill try make a logic order to forest plot based in a collum with cummulative effect
res.fe2 <- rma(yi, vi,
               data= dat.test, method = "DL", slab = paste(ID, sep=","))
rcu <- cumul(res.fe2, order = order(dat.test$year_interval))
x11()
forest(rcu, atransf = exp)
text(-3, 29, "ID", pos = 4, font = 2)
text(3, 29, "Effect Size [95% CI]", pos = 2, font = 2)


#Test for AgrXUrb## Não deu certo
test_use <- read.table(here::here ("data","processed","test_use.txt"), h=T)
colnames(test_use)

use.test<- escalc(measure="SMD", m1i = eff_agr,n1i = N_agr, sd1i = sd_agr,
                   m2i = eff_urb, n2i = N_urb, sd2i = sd_urb,
                   vtype = "LS",
                   data=test_use, append=TRUE)
res.fe3 <- rma(yi, vi,
               data= use.test, method = "DL")

#Test native veg and use as moderators 
data_test <- read.table(here ("data","processed","data_test.txt"), h=T) ##All data
colnames(data_test)
str(data_test)

dat.test <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                   data=data_test, append=TRUE) #método para valores de correlação
res.fe4 <- rma(yi, vi,
               data= dat.test, method = "DL", slab = paste(veget, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
res.fe4.1 <- rma(yi, vi,
               data= dat.test, method = "DL",mods = cbind(veget)-1, slab = paste(veget, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.fe4.1, atransf = exp)
text(-3.5, 29, "Native Vegetation", pos = 4, font = 2)
text(3, 29, "Effect Size [95% CI]", pos = 2, font = 2)



res.fe5 <- rma(yi, vi,
               data= dat.test, method = "DL",mods = cbind(land_use)-1, slab = paste(land_use, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.fe5, atransf = exp)
text(-3.5, 29, "Land Use", pos = 4, font = 2)
text(3, 29, "Effect Size [95% CI]", pos = 2, font = 2)


fsn (yi, vi, data=dat.test, type="Orwin", target= -0.132) ## Para obter quantos estudos ainda precisariam para deixar de ser significativo

x11()
funnel (res.fe1)


####Trying make a single robumeta model#### 
#Fishes
data_fish <- read.table(here ("data","processed","fishes.txt"), h=T) ##All data
colnames(data_fish)
str(data_fish)

dat.fish<- escalc(measure="COR", ri=Effect_size, ni=N_sample, vtype = "LS",
                     data=data_fish, append=TRUE) #método para valores de correlação

model_fish<- robu(yi ~ transition-1,
                     data = dat.fish, 
                     modelweights = "HIER", studynum = ID,  
                     var.eff.size = vi , small = T)
print (model_fish)

fit.fish <- rma(yi, vi, slab = paste(Original_cover,Land_use_class, sep = "->"),
               data= dat.fish, method = "DL")

x11()
forest(fit.fish, atransf = exp)
text(-3.5, 43, "Transition", pos = 4, font = 2)
text(3, 43, "Effect Size [95% CI]", pos = 2, font = 2)

x11()
funnel (fit.fish)


###Macroinvertebrates

data_macroinv <- read.table(here ("data","processed","macroinvertebrates.txt"), h=T) ##All data
colnames(data_macroinv)
str(data_macroinv)

dat.macroinv<- escalc(measure="COR", ri=Effect_size, ni=N_sample, vtype = "LS",
                     data=data_macroinv, append=TRUE) #método para valores de correlação

model_macroinv<- robu(yi ~ transition-1,
                     data = dat.macroinv, 
                     modelweights = "HIER", studynum = ID,  
                     var.eff.size = vi , small = T)
print (model_macroinv)

res.macroinv <- rma(yi, vi, slab = paste(Original_cover,Land_use_class, sep = "->"),
                 data= dat.macroinv, method = "DL")
x11()
forest(res.macroinv, atransf = exp)
text(-3.5, 87, "Transition", pos = 4, font = 2)
text(3, 87, "Effect Size [95% CI]", pos = 2, font = 2)

