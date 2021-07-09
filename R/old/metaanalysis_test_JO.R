#####Charge packages#####
library (metafor)
library(here)
library (robumeta)

####Meta-analysis and Meta-regression#### 
#Macroinvertebrates#

data_macroinv <- read.table(here ("data","processed","macroinvertebrates.txt"), h=T) ##All data
colnames(data_macroinv)
str(data_macroinv)

#Jean: Mateus use o z de Fisher (measure="ZCOR"). O "r" como tamanho de efeito apresenta um vies em sua variancia pois ela eh dependente do proprio valor do r (vr = ((1-r^2)^2)/(n-1))
dat.macroinv<- escalc(measure="ZCOR", ri=Effect_size, ni=N_sample, vtype = "LS",
                     data=data_macroinv, append=TRUE) #método para valores de correlação

#Jean: Antes do modelo abaixo, sugiro voce calcular um tamanho de efeito acumulado
#entre todos os estudos:
cum_effect <- robu(yi ~1,
                   data = dat.macroinv, 
                   modelweights = "CORR", #Jean: Use os pesos correlacionados ao inves dos hierarquicos
                   studynum = ID,  
                   var.eff.size = vi , small = T)
cum_effect
#Esse artigo aquitem uma explicacao "palpavel" sobre a diferenca entre pesos correlacionados
#e pesos hierarquicos: Tanner-Smith & Tipton (2013). Res. Synth. Methods. https://doi.org/10.1002/jrsm.1091
#(veja o exemplo 1 e exemplo 2 na pg. 14)

#Jean: Nesse seu modelo abaixo voce retirou o intercepto, assim o teste dos coeficientes de cada nivel
#de "transition" eh se o coeficiente difere de zero. Se significativo,
#indica que o coeficiente eh diferente de zero (i.e. nao esta testando
#se os coeficientes diferem entre os niveis de transition)
model_macroinv<- robu(yi ~ transition-1,
                     data = dat.macroinv, 
                     modelweights = "HIER", studynum = ID,  
                     var.eff.size = vi , small = T)

model_macroinv<- robu(yi ~ transition, #Aqui a hipotese testada eh se cada nivel difere do nivel referencia (Intercept que no caso eh "foragr")
                      data = dat.macroinv, 
                      modelweights = "CORR", studynum = ID,  
                      var.eff.size = vi , small = T)
print (model_macroinv) #Mateus, note que na tabela dos coeficientes os niveis "openagr" e "openurb" apresentam df < 4 (sempre o robumeta envia um aviso no rodape "Note: If df < 4, do not trust the results"; muito poucos estudos com esses niveis para uma inferencia segura)

#Jean: devido ao pequeno k (nº de estudos) desses dois niveis, recomendo que voce os exclua da meta-regressao (mas, os inclua no tamanho de efeito acumulado acima):
summary(dat.macroinv$Temp)
summary(dat.macroinv$Prec)
summary(dat.macroinv$S_temp)
summary(dat.macroinv$S_prec) #Ok todas essas variaveis nao tem NA

#Filtre os niveis com poucos estudos:
dim(dat.macroinv)
length(unique(dat.macroinv$ID)) #Seu k

dat.macroinv_filt <- dat.macroinv[dat.macroinv$transition!="openagr",]
dat.macroinv_filt <- dat.macroinv_filt[dat.macroinv_filt$transition!="openurb",]
dim(dat.macroinv_filt)
length(unique(dat.macroinv_filt$ID)) #Excluindo ambos os niveis, voce perde 5 estudos

model_macroinv_redu<- robu(yi ~ transition,
                      data = dat.macroinv_filt, 
                      modelweights = "CORR", studynum = ID,  
                      var.eff.size = vi , small = T)
model_macroinv_redu#"foragr" apresenta um tamanho de efeito negativo e "forurb" eh mais negativo ainda 

#Inclua as covariaveis nesse modelo reduzido:
model_macroinv_redu_covars<- robu(yi ~ transition+Temp+Prec+S_temp+S_prec-1,
                           data = dat.macroinv_filt, 
                           modelweights = "CORR", studynum = ID,  
                           var.eff.size = vi , small = T)
model_macroinv_redu_covars #Mesmo controlando as variaveis climaticas, o que importa pro
#efeito mesmo eh o tipo de transicao ("forurb" tem efeito mais negativo que "foragr")
#Mateus, esse modelo "model_macroinv_redu_covars" eh um modelo "global". Eh legal manter as
#inferencias somente nele, ja que voce tem expectativas teoricas para esse modelo

#For the forestplot

res.macroinv <- rma(yi, vi, slab = paste(Original_cover,Land_use_class, sep = "->"),
                 data= dat.macroinv, method = "DL")
x11()
forest(res.macroinv, atransf = exp)
text(-3.5, 87, "Transition", pos = 4, font = 2)
text(3, 87, "Effect Size [95% CI]", pos = 2, font = 2)

#Mateus, eu acho que o forest plot padrao do metafor ficam estranhos pra gente da ecologia porque temos muitos estudos com varios tamanhos de efeitos...
#Da para montar um manualmente utilizando o plot, yi, vi, funcao "polygon" para desenhar o diamante (tamanho de efeito acumulado)
#e omitindo os nomes dos estudos e os valores numericos dos tamanhos de efeitos e seus CI95%