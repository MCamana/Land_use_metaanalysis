#####Metaanalysis test#####
library (metafor)
library(here)

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

####Now, Ill use the data of RICHNESS as a test. First, Ill make the analisys with all data, just as exploratory test####
####Fishes#### 
##All
data_fishall <- read.table(here ("data","processed","fish_all.txt"), h=T) ##All data
colnames(data_fishall)
str(data_fishall)

dat.fishall<- escalc(measure="COR", ri=Effect_size, ni=N_sample, vtype = "LS",
                   data=data_fishall, append=TRUE) #método para valores de correlação

res.fe1 <- rma(yi, vi, slab = paste(ID),
               data= dat.fishall, method = "DL")
x11()
forest(res.fe1, atransf = exp)
text(-3.5, 34, "Article ID", pos = 4, font = 2)
text(3, 34, "Effect Size [95% CI]", pos = 2, font = 2)

fsn (yi, vi, data=dat.test, type="Orwin", target= -0.132) ## Para obter quantos estudos ainda precisariam para deixar de ser significativo

x11()
funnel (res.fe1)

#####Native VEgetation
fish_veg <- read.table(here ("data","processed","fish_veg.txt"), h=T) ##All data
colnames(fish_veg)
str(fish_veg)

fish.veg <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                   data=fish_veg, append=TRUE) #método para valores de correlação
res.fishveg <- rma(yi, vi,
                 data= fish.veg, method = "DL",mods = cbind(Original_cover)-1, 
                 slab = paste(Original_cover, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.fe4.1, atransf = exp)
text(-3.5, 34, "Native Vegetation", pos = 4, font = 2)
text(3, 34, "Effect Size [95% CI]", pos = 2, font = 2)

#####Land Use

fish_land <- read.table(here ("data","processed","fish_land.txt"), h=T) ##All data
colnames(fish_land)
str(fish_land)

fish.land <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                   data=fish_land, append=TRUE) #método para valores de correlação
res.land <- rma(yi, vi,
                 data= fish.land, method = "DL",mods = cbind(Land_use_class)-1, 
                 slab = paste(Land_use_class, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.land, atransf = exp)
text(-3.5, 34, "Land Use", pos = 4, font = 2)
text(3, 34, "Effect Size [95% CI]", pos = 2, font = 2)

##Temperature

fish_temp <- read.table(here ("data","processed","fish_temp.txt"), h=T) ##All data
colnames(fish_temp)
str(fish_temp)

fish.temp <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                    data=fish_temp, append=TRUE) #método para valores de correlação
res.temp <- rma(yi, vi,
                data= fish.temp, method = "DL",mods = cbind(mean_temp)-1, 
                slab = paste(legend, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot

summary(lm(lm(yi~mean_temp, data = fish.temp))) ##Seems different. Perharps a metaregression will be more effective

x11()
forest(res.temp, atransf = exp)
text(-3.5, 33, "Temperature [°C]", pos = 4, font = 2)
text(3, 33, "Effect Size [95% CI]", pos = 2, font = 2)

##Time
fish_time <- read.table(here ("data","processed","fish_time.txt"), h=T) ##All data
colnames(fish_time)
str(fish_time)

fish.time <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                    data=fish_time, append=TRUE) #método para valores de correlação
res.time <- rma(yi, vi,
                data= fish.time, method = "DL",mods = cbind(Year_diff)-1, 
                slab = paste(Year_diff, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot

x11()
forest(res.time, atransf = exp)
text(-3.5, 32, "Years", pos = 4, font = 2)
text(3, 32, "Effect Size [95% CI]", pos = 2, font = 2)

####Macroinvertebrates####
##All
data_macall <- read.table(here ("data","processed","mac_all.txt"), h=T) ##All data
colnames(data_macall)
str(data_macall)

dat.macall<- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                     data=data_macall, append=TRUE) #método para valores de correlação

res.macall <- rma(yi, vi, slab = paste(ID),
               data= dat.macall, method = "DL")
x11()
forest(res.macall, atransf = exp)
text(-3.5, 168, "Article ID", pos = 4, font = 2)
text(3, 168, "Effect Size [95% CI]", pos = 2, font = 2)

#####Native VEgetation
mac_veg <- read.table(here ("data","processed","mac_veg.txt"), h=T) ##All data
colnames(mac_veg)
str(mac_veg)

mac.veg <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                   data=mac_veg, append=TRUE) #método para valores de correlação
res.macveg <- rma(yi, vi,
                   data= mac.veg, method = "DL",mods = cbind(Original_cover)-1, 
                   slab = paste(Original_cover, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.macveg, atransf = exp)
text(-3.5, 77, "Native Vegetation", pos = 4, font = 2)
text(3, 77, "Effect Size [95% CI]", pos = 2, font = 2)

#####Land Use
mac_land <- read.table(here ("data","processed","mac_land.txt"), h=T) ##All data
colnames(mac_land)
str(mac_land)

mac.land <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                    data=mac_land, append=TRUE) #método para valores de correlação
res.macland <- rma(yi, vi,
                data= mac.land, method = "DL",mods = cbind(Land_use_class)-1, 
                slab = paste(Land_use_class, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot
x11()
forest(res.macland, atransf = exp)
text(-3.5, 84, "Land Use", pos = 4, font = 2)
text(3, 84, "Effect Size [95% CI]", pos = 2, font = 2)

##Temperature

fish_temp <- read.table(here ("data","processed","fish_temp.txt"), h=T) ##All data
colnames(fish_temp)
str(fish_temp)

fish.temp <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                    data=fish_temp, append=TRUE) #método para valores de correlação
res.temp <- rma(yi, vi,
                data= fish.temp, method = "DL",mods = cbind(mean_temp)-1, 
                slab = paste(legend, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot

summary(lm(lm(yi~mean_temp, data = fish.temp))) ##Seems different. Perharps a metaregression will be more effective

x11()
forest(res.temp, atransf = exp)
text(-3.5, 33, "Temperature [°C]", pos = 4, font = 2)
text(3, 33, "Effect Size [95% CI]", pos = 2, font = 2)

##Time
fish_time <- read.table(here ("data","processed","fish_time.txt"), h=T) ##All data
colnames(fish_time)
str(fish_time)

fish.time <- escalc(measure="COR", ri=effect_size, ni=N, vtype = "LS",
                    data=fish_time, append=TRUE) #método para valores de correlação
res.time <- rma(yi, vi,
                data= fish.time, method = "DL",mods = cbind(Year_diff)-1, 
                slab = paste(Year_diff, sep=",")) ##I paste the information of native vegetation that will be used in the forest plot

x11()
forest(res.time, atransf = exp)
text(-3.5, 32, "Years", pos = 4, font = 2)
text(3, 32, "Effect Size [95% CI]", pos = 2, font = 2)


