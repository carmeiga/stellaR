setwd("C:/Users/Carlos M Garcia/Google Drive/astronomia")

library(tidyverse)
library(magick)
im <- image_read("Mejide_Garcia_Carlos.gif")

im_proc <- im %>%
  image_channel("blue") %>%
  image_threshold("white", "75%") 

  im_proc2=image_threshold(im_proc,"white", "10%") 
  
  imfin=image_compare(im_proc,im_proc2, metric = "", fuzz = 0) %>% 
  image_threshold("white", "30%")  %>%
image_threshold("black", "80%") 


im_proc3 <- imfin %>%
  image_negate()
imfin

dat <- image_data(im_proc3)[1,,] %>%
  as.data.frame() %>%
  mutate(Row = 1:nrow(.)) %>%
  select(Row, everything()) %>%
  mutate_all(as.character) %>%
  gather(key = Column, value = value, 2:ncol(.)) %>%
  mutate(Column = as.numeric(gsub("V", "", Column)),
         Row = as.numeric(Row),
         value = ifelse(value == "00", NA, 1)) %>%
  filter(!is.na(value))


agg = aggregate(dat,
       
               by = list(dat$Row),
                FUN = mean)[,c(2,3)]

n=nrow(agg)


names(agg)=c('lambda','int')
agg$lambda=seq(360,480, length.out = n)

m=(max(agg$int)-min(agg$int))/2

agg$trans=2*m-agg$int

B=max(agg$trans)
A=min(agg$trans)

D=1.35e-10
C=5.5e-11

agg$fin=((agg$trans-A)*(D-C)/(B-A))+C

paso=agg$lambda[2]-agg$lambda[1]



agg=agg[,c('lambda','fin')]

#ata aqui recupera o espectro da imaxe

########

library(np)
#bw.subset <- npregbw(fin~lambda, bwmethod = "cv.aic",data=agg,regtype = "ll")
#res.np <- npreg(bws = bw.subset)
#np=predict(res.np,newdata=data.frame(lambda=agg$lambda))


####### AXUSTE PARAMETRICO

kb=8.617333262e-11 #MeV/K

hc=2*pi*197*1e-6 #MeV*nm

k=hc/kb

yl=agg$fin
yl=yl/yl[355]
xl=agg$lambda

plancknorm=function(x) {
  lambda0=430.9182 #elemento 355 de xl
  t=30919
  f=((exp(k/(lambda0*t))-1)*lambda0^5)/((exp(k/(x*t))-1)*x^5)
  
  
  return (f)
}



mod.np=npreg(yl ~ xl, bws = 3)
npfit=fitted(mod.np)





plot(xl,yl,type='l',xlab='Lonxitude de onda (nm)',ylab='Intensidade normalizada (-)')
lines(xl,npfit,col='purple')

lambda0=430.9182

y=plancknorm(xl)


m<-nls(npfit~((exp(k/(lambda0*t))-1)*lambda0^5)/((exp(k/(xl*t))-1)*xl^5),start=list(t=20000))


summary(m)

plot(xl,npfit,type='l',col='purple',xlab='Lonxitude de onda (nm)',ylab='Intensidade normalizada (-)')
lines(xl,y,col='darkorange',lty='dashed',lwd='3')


############## PARTE ESQUERDA


y1=yl[1:200]
x1=xl[1:200]



err=y1-npfit[1:200]


plot(x1q,y1q,type='l',xlab='Lonxitude de onda (nm)',ylab='Residuos negativos (-)')

#plot(agg$lambda,err)
plot(x1,y1,type='l')

quitar=(err<0)



points(y1[quitar]~x1[quitar],col='green',pch=17)
lines(x1,npfit[1:200],col="red",lwd=3)


y1q=err[quitar]
y11=y1[quitar]
x1q=x1[quitar]


require(pracma)


pea=findpeaks(-y1q, npeaks=12, threshold=0, sortstr=TRUE)[,2]

plot(x1,y1,type='l',xlab='Lonxitude de onda (nm)',ylab='Intensidade normalizada (-)')

points(x1q[pea],y11[pea],col='blue',pch=16)

linhas1=x1q[pea]


###### PARTE DEREITA


y1=yl[200:length(xl)]
x1=xl[200:length(xl)]



err=y1-npfit[200:length(xl)]



#plot(agg$lambda,err)
plot(x1,y1,type='l')

quitar=(err<0)

points(y1[quitar]~x1[quitar],col='green',pch=17)
lines(x1,npfit[200:length(xl)],col="red",lwd=3)


y1q=err[quitar]
y11=y1[quitar]
x1q=x1[quitar]

plot(x1q,y1q,type='l',xlab='Lonxitude de onda (nm)',ylab='Residuos negativos (-)')


require(pracma)


pea=findpeaks(-y1q, npeaks=14, threshold=0, sortstr=TRUE)[,2]

plot(x1,y1,type='l',xlab='Lonxitude de onda (nm)',ylab='Intensidade normalizada (-)')

points(x1q[pea],y11[pea],col='blue',pch=16)

linhas2=x1q[pea]

round(c(linhas1,linhas2)*10)


